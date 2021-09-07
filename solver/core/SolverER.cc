#include "core/Solver.h"
#include "mtl/Sort.h"

#define MICROSEC_PER_SEC (1000000)

// Timers
static struct rusage g_timer_start, g_timer_end;
inline void timerStart() {
    getrusage(RUSAGE_SELF, &g_timer_start);
}
inline void timerStop(struct rusage& timer) {
    getrusage(RUSAGE_SELF, &g_timer_end);
    
    // Semantics: timer += end_time - start_time

    // Add to total overhead
    timer.ru_utime.tv_sec  += g_timer_end.ru_utime.tv_sec - g_timer_start.ru_utime.tv_sec;
    timer.ru_utime.tv_usec += g_timer_end.ru_utime.tv_usec;

    // Check if subtracting the initial time would result in underflow
    if (g_timer_start.ru_utime.tv_usec > timer.ru_utime.tv_usec) {
        timer.ru_utime.tv_usec += MICROSEC_PER_SEC;
        timer.ru_utime.tv_sec  -= 1;
    }
    timer.ru_utime.tv_usec -= g_timer_start.ru_utime.tv_usec;

    // Check if we carry over to the next second
    if (timer.ru_utime.tv_usec >= MICROSEC_PER_SEC) {
        timer.ru_utime.tv_usec -= MICROSEC_PER_SEC;
        timer.ru_utime.tv_sec  += 1;
    }
}

namespace Minisat {

static inline std::pair<Lit, Lit> mkLitPair(Lit a, Lit b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

static inline void removeLits(std::set<Lit>& set, const std::vector<Lit>& toRemove) {
    for (std::vector<Lit>::const_iterator it = toRemove.begin(); it != toRemove.end(); it++) set.erase(*it);
}

static inline bool contains (std::map< Lit, std::map<Lit, Lit> > extVarMap, Lit a, Lit b) {
    std::map< Lit, std::map<Lit, Lit> >::const_iterator it1 = extVarMap.find(a);
    if (it1 == extVarMap.end()) return false;
    std::map<Lit, Lit>::const_iterator it2 = it1->second.find(b);
    return it2 != it1->second.end();
}

inline void Solver::er_substitute(vec<Lit>& clause, std::map<Lit, std::map<Lit, Lit> >& extVarDefs) {
    // Get set of all literals in clause
    std::set<Lit> learntLits;
    for (int i = 1; i < clause.size(); i++) learntLits.insert(clause[i]);

    // Possible future investigation: should we ignore really long clauses in order to save time?
    // This would need a command-line option to toggle ignoring the clauses

    // Check for extension variables over literal pairs
    for (int i = 1; i < clause.size(); i++) {
        // Check if any extension variables are defined over this literal
        std::map< Lit, std::map<Lit, Lit> >::const_iterator tmp1 = extVarDefs.find(clause[i]);
        if (tmp1 == extVarDefs.end()) continue;

        // TODO: would it be better to iterate over the extension definitions here?
        // Also consider sorting or doing a set intersection here
        const std::map<Lit, Lit>& possibleDefs = tmp1->second;
        for (int j = i + 1; j < clause.size(); j++) {
            // Check if any extension variables are defined over this literal pair
            std::map<Lit, Lit>::const_iterator tmp2 = possibleDefs.find(clause[j]);
            if (tmp2 == possibleDefs.end()) continue;

            // Replace disjunction with intersection
            learntLits.erase(tmp1->first);
            learntLits.erase(tmp2->first);
            learntLits.insert(tmp2->second);
            // printf("(%d v %d) = %d\n", var(tmp1->first) * (sign(tmp1->first) ? -1 : 1), var(tmp2->first) * (sign(tmp2->first) ? -1 : 1), var(tmp2->second) * (sign(tmp2->second) ? -1 : 1));
            goto ER_SUBSTITUTE_DOUBLE_BREAK;
        }
    }
    ER_SUBSTITUTE_DOUBLE_BREAK:;

    // Generate reduced learnt clause
    std::set<Lit>::iterator it = learntLits.begin();
    unsigned int i;
    for (i = 1; i <= learntLits.size(); i++) {
        clause[i] = *it;
        it++;
    }
    clause.shrink(clause.size() - i);
}

void Solver::substituteExt(vec<Lit>& out_learnt) {
    timerStart();
    er_substitute(out_learnt, extVarDefs);
    timerStop(ext_sub_overhead);
}

// Prioritize branching on a given set of literals
inline void Solver::er_prioritize(const std::vector<Var>& toPrioritize) {
    const double desiredActivity = activity[order_heap[0]] * 1.5;
    for (std::vector<Var>::const_iterator i = toPrioritize.begin(); i != toPrioritize.end(); i++) {
        Var v = *i;
        // Prioritize branching on our extension variables
        activity[v] = desiredActivity;
#if EXTENSION_FORCE_BRANCHING
        // This forces branching because of how branching is implemented when ANTI_EXPLORATION is turned on
        // FIXME: this only forces branching on the last extension variable we add here - maybe add a queue for force branch variables?
        canceled[v] = conflicts;
#endif
        if (order_heap.inHeap(v)) order_heap.decrease(v);
    }
}

static inline void insert (std::map< Lit, std::map<Lit, Lit> >& defMap, Lit x, Lit a, Lit b) {
    // Insert for tuple <a, b>
    std::map< Lit, std::map<Lit, Lit> >::iterator it1 = defMap.find(a);
    if (it1 == defMap.end()) {
        std::map<Lit, Lit> submap;
        submap.insert(std::make_pair(b, x));
        defMap.insert(std::make_pair(a, submap));
    } else {
        it1->second.insert(std::make_pair(b, x));
    }

    // Insert for tuple <b, a>
    std::map< Lit, std::map<Lit, Lit> >::iterator it2 = defMap.find(b);
    if (it2 == defMap.end()) {
        std::map<Lit, Lit> submap;
        submap.insert(std::make_pair(a, x));
        defMap.insert(std::make_pair(b, submap));
    } else {
        it2->second.insert(std::make_pair(a, x));
    }
}

inline std::vector<Var> Solver::er_add(
    vec<CRef>& er_def_db,
    std::map< Lit, std::map<Lit, Lit> >& er_def_map,
    const std::map< Var, std::pair<Lit, Lit> >& newDefMap
) {
    // Add extension variables
    // TODO: verify that we do not already have an extension variable for this literal pair before adding clauses
    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v a)
    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v -a)
    std::vector<Var> new_variables;
    for (unsigned int i = 0; i < newDefMap.size(); i++) {
        new_variables.push_back(newVar());
    }

    // Add extension clauses
    for (std::map< Var, std::pair<Lit, Lit> >::const_iterator i = newDefMap.begin(); i != newDefMap.end(); i++) {
        // Get literals
        Lit x = mkLit(i->first);
        Lit a = i->second.first;
        Lit b = i->second.second;
        assert(var(x) > var(a) && var(x) > var(b));

        // Create extension clauses and add them to the extension definition database
        addClauseToDB(er_def_db, ~x, a, b);
        addClauseToDB(er_def_db, x, ~a);
        addClauseToDB(er_def_db, x, ~b);

        // Save definition
        insert(er_def_map, x, a, b);
    }

    return new_variables;
}

// Add extension variables to our data structures and prioritize branching on them.
// This calls a heuristic function which is responsible for identifying extension variable
// definitions and adding the appropriate clauses and variables.
void Solver::addExtVars(
    std::vector<CRef>(*er_select_heuristic)(Solver&, unsigned int),
    std::map< Var, std::pair<Lit, Lit> >(*er_add_heuristic)(Solver&, std::vector<CRef>&, unsigned int),
    unsigned int numClausesToConsider,
    unsigned int maxNumNewVars
) {
    // Get extension clauses according to heuristics
    timerStart();
    std::vector<CRef> candidateClauses = er_select_heuristic(*this, numClausesToConsider);
    const std::map< Var, std::pair<Lit, Lit> > newDefMap = er_add_heuristic(*this, candidateClauses, maxNumNewVars);
    timerStop(ext_add_heuristic_overhead);

    // Add the extension variables to our data structures
    timerStart();
    const std::vector<Var> new_variables = er_add(extDefs, extVarDefs, newDefMap);
    er_prioritize(new_variables); 
    timerStop(ext_add_ds_overhead);
}

void Solver::delExtVars(std::vector<Var>(*er_delete_heuristic)(Solver&)) {
    extTimerStart();

    // Get variables to delete
    std::vector<Var> varsToDelete = er_delete_heuristic(*this);
    std::set<Var> varsToDeleteSet(varsToDelete.begin(), varsToDelete.end());

    // Delete variables

    // TODO: add an option to switch between these two modes

    // option 1: delete all clauses containing the extension variables
    delExtVars(extLearnts, varsToDeleteSet);    
    delExtVars(extDefs, varsToDeleteSet);    

    // option 2: substitute extension variable with definition
    // TODO

    extTimerStop(ext_delV_overhead);
}

}