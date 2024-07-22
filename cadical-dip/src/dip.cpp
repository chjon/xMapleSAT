#include "internal.hpp"

namespace CaDiCaL {

  string Internal::graphLit2String (int l) {
    string tmp = to_string(abs(l));
    if (l < 0) return "-" + tmp;
    else return tmp;
  }
  
  void Internal::writeFinalInfoGraph (ofstream& out, set<int>& current, set<int>& previous, int UIP, bool writePrevious) {
    for (auto lit : current)
      if (lit != UIP) 
	out << "\"" << graphLit2String(lit) << "\"[ fillcolor = lightblue, style = filled ];" << endl;  
    
    if (writePrevious)
      for (auto lit : previous)
	out << "\"" << graphLit2String(lit) << "\"[ fillcolor = lightgrey, style = filled ];" << endl;
    
    out << "\"" << graphLit2String(UIP) << "\"[ fillcolor = red, style = filled ];" << endl;
    out << "}" << endl;
  }
  
  void Internal::writeEdgeInGraph (ofstream& out, const int& orig, const int& dest, bool colored) {
    if (abs(dest) == max_var + 1) out << "\"" << graphLit2String(orig) << "\" -> CONFLICT";
    else out << "\"" << graphLit2String(orig) << "\" -> \"" << graphLit2String(dest) << "\" ";
    if (colored)  out << "[color=red];";
    out << endl;
  }

  void Internal::writeDIPComputationInfo (TwoVertexBottlenecks& info, DIPGraphEncoder& encoder, const vector<int>& predecessors, const vector<int>& predecessorsLits, const vector<int>& predIndex, const vector<int>& literalsInAnalysis, bool foundDIP, const vector<int>& pathA, const vector<int>& pathB) {
    cout << endl;
    cout << endl;
    cout << string(60,'=') << endl;
    cout << "Conflict " << stats.conflicts << " at level " << level << endl;
    cout << "N " << encoder.numVertices() << endl;
    cout << "predecessors: ";
    for (int x : predecessors) cout << x << " ";
    cout << endl;
    cout << "predIndex: ";
    for (int x : predIndex) cout << x << " ";
    cout << endl;
    
    cout << "Literals in analysis: ";
    for (auto lit : literalsInAnalysis) cout << lit << " ";
    cout << endl;
    cout << "Encoding: " << endl;
    for (auto lit : predecessorsLits)
      cout << "int " << lit << " is DIP algorithm vertex " << encoder.Solver2Sam(lit) << endl;

    cout << "Found DIP: " << (foundDIP?"true":"false") << endl;
    
    vector<TwoVertexBottlenecks::VertPairInfo> listA = info.GetVertListA();
    vector<TwoVertexBottlenecks::VertPairInfo> listB = info.GetVertListB();
    
    cout << string(30,'-') << endl;
    cout << "pathA: ";
    for (int x : pathA) cout << x << " ";
    cout << endl;
    cout << "pathB: ";
    for (int x : pathB) cout << x << " ";
    cout << endl << endl;
    
    cout << "List A: " << endl;
    for (uint i = 0; i < listA.size(); ++i){
      cout << "Vertex[" << i << "]: " << listA[i].vertNum << " --> lit " << encoder.Sam2Solver(listA[i].vertNum) << endl;
      cout << "minPair: " << listA[i].minPair << endl; 
      cout << "maxPair: " << listA[i].maxPair << endl;
      cout << "minAncestor: " << listA[i].minAncestor << endl;
    }
    
    cout << endl;
    cout << "List B: " << endl;
    for (uint i = 0; i < listB.size(); ++i){
      cout << "Vertex[" << i << "]: " << listB[i].vertNum << " --> lit " << encoder.Sam2Solver(listB[i].vertNum) << endl;
      cout << "minPair: " << listB[i].minPair << endl; 
      cout << "maxPair: " << listB[i].maxPair << endl;
      cout << "minAncestor: " << listB[i].minAncestor << endl;
    }
  }
  
  
  Internal::DIPGraphEncoder::DIPGraphEncoder (int nVars ):
    nextVertex(0), _solver2sam(nVars+2,DIP_ENCODER_UNDEF) { }

  inline int Internal::DIPGraphEncoder::numVertices ( )      const { return nextVertex; }
  inline int Internal::DIPGraphEncoder::Sam2Solver (int var) const { return _sam2solver[var]; };
  inline int Internal::DIPGraphEncoder::Solver2Sam (int lit)       {
    int v = abs(lit);
    if (_solver2sam[v] == DIP_ENCODER_UNDEF) {
      _solver2sam[v] = nextVertex;
      _sam2solver.push_back(lit);
      ++nextVertex;
      //      cout << "Lit " << lit << " gets " << nextVertex - 1 << std::endl;
      return nextVertex-1;
    }
    else return _solver2sam[v];
  };


  inline void Internal::analyze_literal_dip (int lit, int &open) {
    assert (lit);
    Var &v = var (lit);
    Flags &f = flags (lit);
    
    if (!v.level) return;
    if (f.seen) return;

    f.seen = true;
    analyzed.push_back (lit);
    
    assert (val (lit) < 0);
    assert (v.level <= level);
    assert (v.reason != external_reason);
    Level &l = control[v.level];
    if (v.level < level)
      clause.push_back (lit);
    if (!l.seen.count++) { // I keep this in order to get the LBD
      LOG ("found new level %d contributing to conflict", v.level);
      levels.push_back (v.level);
    }
    if (v.trail < l.seen.trail) // this is needed for lemma minimization only (I think)
      l.seen.trail = v.trail;
    LOG ("analyzed literal %d assigned at level %d", lit, v.level);
    if (v.level == level)
      open++;
  }
  
  inline void Internal::analyze_reason_dip (int lit, Clause *reason, int &open) {
    assert (reason);
    assert (reason != external_reason);    
    // cout << "Explore reason(" << lit;
    // if (lit) cout << "-h" << var(lit).trail;
    // cout << "): ";
    // for (int x : *reason) cout << x << " (h" << var(x).trail << ") ";
    // cout << endl;
    for (const auto &other : *reason)
      if (other != lit)
	analyze_literal_dip (other, open);
  }
  
  // Finds the 1UIP. This function has no side effects: no backtracking, no statistics,
  // no activity update, etc.
  
  bool Internal:: try_dip_analysis ( vector<int>& dip2conflict, vector<int>& uip2dip, vector<int>& analyzed_lits){
    //cout << "Call to find1UIP at DL " << level << " and control size is " << control.size() << endl;
    assert(clause.empty ());
    assert(check_all_literals_cleared());
    assert(check_all_levels_cleared());


    
    // Stuff for writing the conflict graph (useful for debugging/understanding)
    // Currently the writing is interleaved with the 1UIP computation
    // It might be cleaner to do it separately (but I think that is not a problem for efficiency)
    string filename1 = "graph-"+to_string(stats.conflicts) + ".dot"; //conflict graph
    string filename2 = "graph-current-"+to_string(stats.conflicts) + ".dot"; // conflict graph with only last-DL lits
    ofstream outAll;
    ofstream outCurrent;
    set<int> currentNodes;
    set<int> previousNodes;
    
    // if (stats.conflicts == 6730534) exit(1);
    // bool write = stats.conflicts == 6730533; // To write conflict graph only of a concrete conflict
#define write 0// quicker for release mode
    if (write) {
      outAll.open(filename1.c_str(),fstream::out);
      outCurrent.open(filename2.c_str(),fstream::out);
      outAll << "digraph D {" << endl;
      outCurrent << "digraph D {" << endl;
    }

    // Data structure for generating correct input for DIP-detection algorithm
    DIPGraphEncoder encoder(max_var);
    vector<int> predecessorsLits; // only needed for DIP-detection algorithm
    vector<int> predIndex; // only needed for DIP-detection algorithm
    vector<int> literalsInAnalysis; // only needed for DIP-detection algorithm
    // The latter stores the literals as they have appeared in the analysis (inverse topological order)
    // This will allow us to assign numbers respecting the order (the algorithm needs this)
    // Remember that conflict node is the sink (node 0), and the 1UIP is the source (node with largest N)
    // Once the analysis is done, we only have to convert the predecessorLits into a vector<int>
    // to get the correct input to DIP-detection

    // Start by adding the edges from the negation of the lits  the conflicting clause to the "CONFLICT" node
    int conflictLit = max_var + 1;
    literalsInAnalysis.push_back(conflictLit);// Fake literal corresponding to conflict
    predIndex.push_back(predecessorsLits.size());      
    Clause* confClause = conflict;
    for (int lit : *confClause) {
      if (write) writeEdgeInGraph(outAll,-lit,conflictLit, confClause->redundant);
      if (var(lit).level == level) {
	predecessorsLits.push_back(-lit);
	if (write) {
	  currentNodes.insert(-lit);
	  writeEdgeInGraph(outCurrent,-lit,conflictLit,confClause->redundant);
	}
      }
      else if (write) previousNodes.insert(-lit);
    }

    Clause *reason = conflict;
    LOG (reason, "analyzing conflict");
    
    assert (clause.empty ());
    assert (lrat_chain.empty ());
    clause.push_back(0); // space for uip
    const auto &t = &trail;
    int i = t->size ();      // Start at end-of-trail.
    int open = 0;            // Seen but not processed on this level.
    int uip = 0;             // The first UIP literal.

    for (;;) {
      analyze_reason_dip (uip, reason, open);
      
      uip = 0;
      while (!uip) {
	assert (i > 0);
	const int lit = (*t)[--i];
	if (!flags (lit).seen)  continue;
	if (var (lit).level == level) uip = lit;
      }
      if (!--open) {
	literalsInAnalysis.push_back(uip);
	predIndex.push_back(predecessorsLits.size());
	break;
      }
      literalsInAnalysis.push_back(uip);
      predIndex.push_back(predecessorsLits.size());	    
      reason = var (uip).reason;
      assert (reason != external_reason);
      LOG (reason, "analyzing %d reason", uip);

      for (int lit : *reason) {
	if (lit != uip) {
	  if (write) writeEdgeInGraph(outAll,-lit,uip,reason->redundant);
	  if (var(lit).level == level) {
	    predecessorsLits.push_back(-lit);
	    if (write) {
	      currentNodes.insert(-lit);
	      writeEdgeInGraph(outCurrent,-lit,uip,reason->redundant);
	    }
	  }
	  else if (write) previousNodes.insert(-lit);
	}
      }      
    }

    clause[0] = -uip;
    // cout << "Lemma 1UIP (conflict " << stats.conflicts << "): ";
    // for (int x : clause) cout << x << "(v " << int(val(x)) << ", lev " << var(x).level << ") ";
    // cout << endl;

    minimize_clause();

    // cout << "Lemma minimized: ";
    // for (int x : clause) cout << x << " ";
    // cout << endl;

    // This is super slow!!!!!
    // We place "uip" at the position 0
    int pos = -1;
    for (uint i = 0; pos == -1 and i < clause.size(); ++i) {
      if (clause[i] == -uip) pos = i;
    }
    assert(pos != -1);
    swap(clause[0],clause[pos]);
    
    // cout << "UIP " << uip << endl;
    if (write) {
      cout << "After placing it properly: ";
      for (int x : clause) cout << x << " ";
      cout << endl;
    }

    
    vector<int> UIP_clause = clause;
    analyzed_lits = analyzed;
    
    if (write) {
      writeFinalInfoGraph(outAll,currentNodes,previousNodes,uip,true);
      writeFinalInfoGraph(outCurrent,currentNodes,previousNodes,uip,false);
      outAll.close();
      outCurrent.close();
    }
    
    clear_analyzed_levels ();
    clear_analyzed_literals ();
    clause.clear(); // Here we should probably minimize this clause and store it somewhere
                    // We need to clear it for DIP clause computation that go next
    assert(check_all_levels_cleared());
    
    LOG ("first UIP %d", uip);
    
    vector<int> predecessors;
    for (auto lit : literalsInAnalysis) encoder.Solver2Sam(lit); // encode in topological order
    for (auto lit : predecessorsLits) predecessors.push_back(encoder.Solver2Sam(lit)); // map vector<int> to vec<int>

    TwoVertexBottlenecks dip;
    vector<int> pathA, pathB; // the two disjoint paths
    int res = dip.CalcBottlenecks(predecessors,predIndex,pathA,pathB);
    
    bool foundDIP =  (res > 0);
    if (foundDIP) ++stats.dip_exists;
      
    //if (write) writeDIPComputationInfo(dip,encoder,predecessors,predecessorsLits,predIndex,literalsInAnalysis,foundDIP,pathA,pathB);

    if (foundDIP) {
      res -=4;
      int b = res%2;
      res/=2;
      int a = res;     
      
      vector<int> lits_to_bump;
      bool ok = computeDIPClauses(a,b,conflict,dip,encoder,dip2conflict,uip2dip,uip,pathA, pathB, lits_to_bump);
      assert(check_all_literals_cleared());
      assert(check_all_levels_cleared());
      if (ok) {
	analyzed_lits = lits_to_bump;
	return true;
      }
      else {
	dip2conflict = UIP_clause;
	return false;
      }
    }
    

    assert(check_all_literals_cleared());
    assert(check_all_levels_cleared());
    dip2conflict = UIP_clause;
    return false;

  }


  bool Internal::computeClosestDIPToConflict (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, int& x, int& y) {
    const vector<TwoVertexBottlenecks::VertPairInfo>& listA = info.GetVertListA();
    
    // a == 1 if the first element of VertPairListA is an immediate ancestor of the conflict node
    // b == 1 if the first element of VertPairListB is an immediate ancestor of the conflict node
    
    // Closest to the conflict
    int idxA;
    if (a == 1) {
      if (listA.size() == 1) return false;
      else idxA = 1;
    }
    else idxA = 0;
    
    const vector<TwoVertexBottlenecks::VertPairInfo>& listB = info.GetVertListB();    
    int idxB = -1;
    for (uint i = 0; i < listB.size(); ++i) {
      if (b == 1 and i == 0) continue;
      if (listA[idxA].minPair <= listB[i].vertNum and
	  listB[i].vertNum <= listA[idxA].maxPair and
	  listA[idxA].minAncestor > listB[i].vertNum) { // avoid ancestors but I believe is not necessary
	idxB = i;
      }
    }
    
    if (idxB == -1) return false;
    
    x = encoder.Sam2Solver(listA[idxA].vertNum);
    y = encoder.Sam2Solver(listB[idxB].vertNum);
    
    return true;
  }

  bool Internal::computeRandomDIP (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, int& x, int& y) {
    // a == 1 if the first element of VertPairListA is an immediate ancestor of the conflict node
    // b == 1 if the first element of VertPairListB is an immediate ancestor of the conflict node
    
    vector<TwoVertexBottlenecks::VertPairInfo> listA = info.GetVertListA();
    int idxA;
    if (a == 1) {
      if (listA.size() == 1) return false;
      else idxA = rand()%(listA.size()-1) + 1;
    }
    else idxA = rand()%listA.size();
    
    x = encoder.Sam2Solver(listA[idxA].vertNum);
    
    vector<TwoVertexBottlenecks::VertPairInfo> listB = info.GetVertListB();  
    vector<int> candidatesY;
    for (uint i = 0; i < listB.size(); ++i){
      if (b == 1 and i == 0) continue;
      if (listA[idxA].minPair <= listB[i].vertNum and
	  listB[i].vertNum <= listA[idxA].maxPair)
	{
	  if (listA[idxA].minAncestor <= listB[i].vertNum) { // skip ancestors (I do not think this is necessary)
	    //cout << "Skip candidate " << encoder.Sam2Solver(listB[i].vertNum) << " because it is an ancestor of " << x << endl;
	  }
	  else candidatesY.push_back(encoder.Sam2Solver(listB[i].vertNum));
	}
    }
    
    if (candidatesY.size() == 0) return false;
    y = candidatesY[rand()%candidatesY.size()];
    
    return true;
  }
  
  bool Internal::computeBestMiddleDIP (const TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, const vector<int>& pathA, const vector<int>& pathB, int& x, int& y) {
    
    if (pathA.size() <= 3 and pathB.size() <= 3) return false; // does not seem interesting
  
    
    const vector<TwoVertexBottlenecks::VertPairInfo>& listA = info.GetVertListA();
    const vector<TwoVertexBottlenecks::VertPairInfo>& listB = info.GetVertListB();
    
    // The idea is:
    // 1) Assume wlog that pathA is the longest path between pathA, pathB
    // 2) Among all DIPs in pathA, choose the one (DIP_A) that is closest to the middle position of pathA
    //    (i.e. ideally half the vertex in pathA)
    // 3) Having chosen DIP_A, choose DIP_B a vertex in pathB that, among the ones that form a DIP with DIP_A,
    //    is the one which is closes to the middle position inside pathB
    
    // Determine who is longest/shortest
    const vector<int>& longestPath = (pathA.size() > pathB.size() ? pathA : pathB);
    const vector<int>& shortestPath = (pathA.size() > pathB.size() ? pathB : pathA);
    const vector<TwoVertexBottlenecks::VertPairInfo>& longestList = (pathA.size() > pathB.size() ? listA : listB);
    const vector<TwoVertexBottlenecks::VertPairInfo>& shortestList = (pathA.size() > pathB.size() ? listB : listA);
    
    // First we compute the position in the path of each vertex of the lists (which are the DIP candidates)
    // For every k, posInLongestPath[k] = p  iff vertex longestList[k].vertNum appears in longestPath[p]
    // similarly for shortest
    vector<int> posInLongestPath(longestList.size());
    vector<int> posInShortestPath(shortestList.size());
    
    // Remember that pathA, pathB both start with the 1UIP (largest vertex) and finish with the conflict (vertex 0).
    // We know that both paths are strictly ordered from large to small
    // However, note that longest/shortest List are ordered strictly from smaller to larger
    for (uint i = 0, j = longestPath.size() - 1; i < longestList.size(); ++i) {
      while (longestPath[j] != longestList[i].vertNum) --j;
      posInLongestPath[i] = j;
      //cout << "Pos of " << longestList[i].vertNum << " in longestPath is " << j << endl;
    }
    
    for (uint i = 0, j = shortestPath.size() - 1; i < shortestList.size(); ++i) {
      while (shortestPath[j] != shortestList[i].vertNum) --j;
      posInShortestPath[i] = j;
      //cout << "Pos of " << shortestList[i].vertNum << " in shortestPath is " << j << endl;
    }
    
    int medianLongest = longestPath.size()/2;
    int medianShortest = shortestPath.size()/2;
    
    // We know that the longest path has size >= 3 (we have returned false otherwise)  
    int idxLongest = 0;
    int longestDistToMedian  = abs(medianLongest - posInLongestPath[0]);
    // Look for the one closest to the middle (DIP_A in the description)
    uint k = 1;
    while (k < posInLongestPath.size() and longestDistToMedian > 0) {
      int dist = abs(medianLongest - posInLongestPath[k]);
      if (dist < longestDistToMedian) {
	longestDistToMedian = dist;
	idxLongest = k;
      }
      ++k;
    }
    
    //cout << "For longest we choose DIP " << longestList[idxLongest].vertNum << " at position " << posInLongestPath[idxLongest] << endl;
    
    // minV and maxV tell that any number inside [minV, maxV] is DIP-compatible with DIP_A
    int minV = longestList[idxLongest].minPair;
    int maxV = longestList[idxLongest].maxPair;
    
    int idxShortest = -1;
    int shortestDistToMedian = 1e9;
    
    k = 0;
    while (k < posInShortestPath.size() and shortestDistToMedian > 0) {
      if (shortestList[k].vertNum < minV or shortestList[k].vertNum > maxV) {++k; continue;}
      int dist = abs(medianShortest - posInShortestPath[k]);
      if (dist < shortestDistToMedian) {
	shortestDistToMedian = dist;
	idxShortest = k;
      }
      ++k;
    }
    
    assert(idxShortest >= 0);
    // if (idxShortest < 0) {  // Just in case. If assertion fails sometime, add this code
    //   return false; 
    // }
    
    //cout << "For shortest we choose DIP " << shortestList[idxShortest].vertNum << " at position " << posInShortestPath[idxShortest] << endl;
    
    
    if (posInShortestPath[idxShortest] == int(shortestPath.size()) - 2 and
	posInLongestPath[idxLongest] == int(longestPath.size()) - 2)
      {
	//cout << "TWO IMMEDIATE ANCESTORS" << endl;
	return false;}
    
    
    x = encoder.Sam2Solver(longestList[idxLongest].vertNum);
    y = encoder.Sam2Solver(shortestList[idxShortest].vertNum);
        
    return true;
    
  }

  bool Internal::computeDIPClauses (int a, int b, Clause* confl, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, vector<int>& clause_to_learn, vector<int>& clause_to_learn2, int uip, vector<int>& pathA, vector<int>& pathB, vector<int>& lits_to_bump) {
    assert(clause.empty());
    assert(check_all_levels_cleared());
    assert(check_all_literals_cleared());

    //bool write = stats.conflicts == 6730533; // To write conflict graph only of a concrete conflict
    //    Clause* originalConflict = confl; 
    
    int x, y; // {x,y} are a DIP
    // 3 DIP computations: random, closest to conflict, the one in the middle

    bool dip_found = false;
    static const int MIDDLE_DIP = 1;
    static const int CLOSEST_TO_CONFLICT = 2;
    static const int RANDOM_DIP = 3;

    int dip_type = MIDDLE_DIP;
    //int dip_type = CLOSEST_TO_CONFLICT;
    // int dip_type = RANDOM_DIP;
  
    if (dip_type == MIDDLE_DIP) dip_found = computeBestMiddleDIP(info,encoder,pathA,pathB,x,y);
    else if (dip_type == CLOSEST_TO_CONFLICT) dip_found = computeClosestDIPToConflict(a,b,info,encoder,x,y);
    else if (dip_type == RANDOM_DIP) dip_found = computeRandomDIP(a,b,info,encoder,x,y);
    else assert(false);

    if (not dip_found) return false;
    
    // cout << endl << endl;
    if (write) {
      cout << "DIP " << x << " height " << var(x).trail << " and " << y << " height " << var(y).trail << " at conflict " << stats.conflicts << endl;
    }

    // Decidir si es prou bo.
    notify_DIP(x,y);

    //cout << "DIP: " << x << " " << y << endl;
    
    if (occurrences_DIP(x,y) < 20) return false;

    int ext_var = er_manager.find_definition(-x,-y);
    if (ext_var != 0) {
      if (flags(ext_var).fixed()) { 
	assert(val(ext_var) > 0); // in this case we should be able to remove the DIP from the clause!
	++stats.dip_dangerous;
	return false;
      }

      if (int(val(ext_var)) != 0 and var(ext_var).level == 0) {	++stats.dip_dangerous;return false;}

      if (not(val(ext_var) >= 0 or var(ext_var).level == level)) {	++stats.dip_dangerous;return false;} // due to chronological BT
      //assert(val(z) >= 0 or var(z).level == level); // If false it has to be at the same DL as x and y
      
      if (val(ext_var) > 0) {
	//cout << "CAREFUL: extended variable is true but should be false/undef" << endl;
	++stats.dip_dangerous;
	return false;
      }

    }
    
    // if (dip_pair_threshold != -1) { // common pair used
    //   solver.branchingHeuristicManager.notifyDIPCandidateCommonPair(~x,~y); 
    //   if (not solver.branchingHeuristicManager.isPairCommon(~x,~y,dip_pair_threshold)) return false;
    // }
    // else if (dip_window_size != -1) { // use sliding window
    //   solver.branchingHeuristicManager.notifyDIPCandidateWindow(~x,~y);
    //   if (not solver.branchingHeuristicManager.isPairBetterThanAverage(~x,~y)) return false;
    // }

    // 1) Compute first clause: DIP ^ after -> conflict
    // This is done starting from conflict, doing standard 1UIP reasoning but stopping as soon as we hit the
    // two elements of the DIP

    int dip_reached = 0;
    const auto &t = &trail;
    int i = t->size();
    int p = 0;
    int open = 0; // do not use here, I believe
    Clause *reason = confl;
    clause.push_back(0); // space for DIP variable
    assert(clause.size() == 1);

    for ( ; dip_reached != 2; ) {
      analyze_reason_dip(p,reason,open);
      
      while (dip_reached != 2) {
	p = 0;
	const int lit = (*t)[--i];
	if (!flags(lit).seen) continue;
	if (var(lit).level == level) {
	  p = lit;
	  if (p != x and p != y) break;
	  else ++dip_reached;
	}
      }
      reason = var(p).reason;
      assert (reason != external_reason);
    }

    // Pensar que passa amb els analitzats perque s'ha d'incrementar heuristica
    // for (Lit z : toIncreaseActivity)
    //   branchingHeuristicManager.handleEventLitInConflictGraph(z, solver.conflicts);

    clause_to_learn = clause; // do it before add_extended_var since this function modifies clause
    lits_to_bump = analyzed;
    clear_analyzed_levels();
    clear_analyzed_literals();
    clause.clear();
    
    // z <--> -x v -y
        
    int z = add_extended_var(-x,-y);

    if (write) 
      cout << "Def " << z << " for " << -x << " v " << -y << endl;

    // Check whether z appears (repeated) or -z appears (tautology)
    for (uint i = 1; i < clause_to_learn.size(); ++i) {
      if (abs(clause_to_learn[i]) == z) {
	if (clause_to_learn[i] == -z) {	++stats.dip_dangerous;return false;}
	else { // Remove repeaded literal
	  cout << "REPEATED LITERAL" << endl;
	  clause_to_learn[i] = clause_to_learn.back();
	  clause_to_learn.pop_back();
	}
      }
    }
    // cout << "DL " << level << endl;
    // cout << "Extended definition: " << z << "(val " << int(val( z)) << ", lev " << var(z).level << ")  <--> ";
    // cout                           << -x << "(val " << int(val(-x)) << ", lev " << var(x).level << ") v ";
    // cout                           << -y << "(val " << int(val(-y)) << ", lev " << var(y).level << ")" << endl;
    assert(val(x) > 0); // x is true
    assert(var(x).level == level);
    assert(val(y) > 0); // y is true
    assert(var(y).level == level);

    clause_to_learn[0] = z;
    if (write) {
      cout << "DIP -> conflict clause is: ";
      for (int x : clause_to_learn) cout << x << "(val " << int(val(x)) << ",lev " << var(x).level << ") ";
      cout << endl;
    }

    assert(clause.empty());

    // if (newDef and solver.produce_proof) { // Add definition to proof
    //   static vec<Lit> ERClause;
    //   ERClause.clear();
    //   ERClause.push(~extLit);     ERClause.push(~x);     ERClause.push(~y); 
    //   solver.proofLogger.addClause(ERClause);

    //   ERClause.clear();
    //   ERClause.push(extLit);    ERClause.push(x);     
    //   solver.proofLogger.addClause(ERClause);
    
    //   ERClause.clear();
    //   ERClause.push(extLit);    ERClause.push(y);     
    //   solver.proofLogger.addClause(ERClause);
    // }

    return true;
    // 3) UIP ^ before --> DIP
    // This is done as standard 1UIP but starting with the two DIP lits being marked
    // Also, we start having a look at the trail at the max height of the DIPs lits  

    int heightX = var(x).trail;
    int heightY = var(y).trail;
    // cout << x << " at height " << heightX << endl;
    // cout << y << " at height " << heightY << endl;
    assert((*t)[heightX] == x and (*t)[heightY] == y);

    assert(clause.empty());
    clause.push_back(0); // space for UIP
    clause.push_back(0); // space for DIP
	
    p = (heightX > heightY ? x : y);
    i = max(heightX,heightY);
    reason = var(p).reason;
    analyzed.push_back(x); flags(x).seen = true; // think whether this is necessary
    analyzed.push_back(y); flags(y).seen = true;
    open = 1;
    for(;;){
      analyze_reason_dip (p, reason, open);
      p = 0;
      while (!p) {
	assert (i > 0);
	const int lit = (*t)[--i];
	//	cout << "Trobo " << lit << " amb seen " << flags(lit).seen << endl;
	if (!flags (lit).seen)  continue;
	if (var (lit).level == level) p = lit;
      }
      if (!--open) break;
      reason = var(p).reason;
      assert (reason != external_reason);      
    }

    clause_to_learn2 = clause;
    clause_to_learn2[0] = -uip; // UIP
    clause_to_learn2[1] = -z; // DIP
#ifndef NDEBUG
    for (uint i = 2; i < clause_to_learn2.size(); ++i)
      assert(var(clause_to_learn2[i]).level < level);
#endif
    if (write) {
      cout << "UIP --> DIP clause is: "; // Falta posar UIP i DIP a la clausula
      for (int x : clause_to_learn2) cout << x << "(val " << int(val(x)) << ",lev " << var(x).level << ") ";
      cout << endl;
    }
    clear_analyzed_levels();
    clear_analyzed_literals();
    clause.clear();    

    // clause_to_learn2.push(~p);
    // clause_to_learn2.push(~extLit); // Two first lits are of highest DL
    // // First is UIP
    // assert(assignmentTrail.level(var(clause_to_learn2[0])) == assignmentTrail.decisionLevel());  
    // for (auto l : beforeLits) {
    //   clause_to_learn2.push(l);
    //   assert(assignmentTrail.level(var(l)) < assignmentTrail.decisionLevel());
    // }

    // // Clear marks
    // for (int i = 0; i < clause_to_learn2.size(); ++i) seen3[var(clause_to_learn2[i])] = false;
    // seen3[var(p)] = seen3[var(x)] = seen3[var(y)] = false;
  
    // assert(checkSeen3());
  
    return true;
  }

  void Internal::substitute_definitions_in_clause ( ){
    // PRE: UIP is in the first position
    
    for (uint i = 1; i < clause.size(); ++i) {
      if (not er_manager.part_of_definition(clause[i])) continue;
      for (uint j = i + 1; j < clause.size(); ++j) {
	if (not er_manager.part_of_definition(clause[j])) continue;
	int def = er_manager.find_definition(clause[i], clause[j]);
	if (def != 0){
	  assert(val(clause[i]) < 0);
	  assert(val(clause[j]) < 0);
	  if (val(def) == 0) continue;
	  if (var(def).level == level) continue;
	  // def is defined at some previous level
	  
	   // cout << "Found at level " << level << " -- " << clause[i] << "(pos " << i << ") and " << clause[j] << "(pos " << j << ") that are part of definition " << er_manager.find_definition(clause[i],clause[j]) << endl;
	   // cout << "Clause is: ";
	   // for (int x : clause) 
	   //   cout << x << "(val " << int(val(x)) << ", lev " << var(x).level << ") ";
	   // cout << endl;


	   bool tautology = false;
	   int  def_found = -1;
	   for (uint k = 0; k < clause.size(); ++k){
	     int x = clause[k];
	     if (x == -def) tautology = true;
	     if (x == def) def_found = k;
	   }
	   
	   if (tautology) continue;
	   
	   vector<int> orig_clause = clause;
      
	     
	   if (def_found != -1) { // Remove positions i and j
	     clause[i] = clause.back(); clause.pop_back();
	     clause[j] = clause.back(); clause.pop_back();
	   }
	   else {
	     clause[i] = def;	 
	     clause[j] = clause.back();
	     clause.pop_back();
	   }
	   
#ifndef NDEBUG
	   set<int> lits;
	   for (int x : clause) {
	     if (lits.count(x) != 0) {
	       cout << "Repetit " << x << ": ";
	       for (int y : clause) cout << y << " ";
	       cout << endl;
	       assert(false);
	     }
	     lits.insert(x);
	   }
#endif
	   // cout << "After substitution (level " << level << "): ";
	   // for (int x : clause) cout << x << "(val " << int(val(x)) << ", lev " << var(x).level << ") ";
	   // cout << endl << endl;
	   
	   // if (val(def) >= 0) {
	   //   cout << "We are at DL " << level << endl;
	   //   for (auto c : clauses){
	   //     for (auto it = c->begin(); it != c->end(); ++it)
	   // 	cout << *it << " (val " << int(val(*it)) << ", lev " << var(*it).level << ") ";
	   //     cout << endl;
	   //   }
	   // }
	   
	   //clause = orig_clause;
	   //return;
	}
      }
      
    }
  }
  
} // namespace CaDiCaL



