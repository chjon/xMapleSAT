#include "internal.hpp"

namespace CaDiCaL {
  // return z such that z <--> x v y
  int Internal::add_extended_var (int x, int y) {
    // First add extended variable
    int z = er_manager.add_extended_var(x,y,max_var+1);
    if (z == max_var + 1) {
      // This is a new variable and we should add clauses
      //      cout << "Create new definition for " << z << endl;
      external->init_er(max_var + 1);
      mark_active(z);
      mark_var_actively_deleted(z,false);
      
      assert (flags(z).status == Flags::ACTIVE);
      
      assert(val(x) < 0);
      assert(val(y) < 0);

      //cout << -z << " " << x << " " << y << endl;
      
      clause = {-z,x,y};  new_extended_definition_clause();

      // assert(val(z) <= 0);
      // if (val(z) == 0)
      //   search_assign(-z,new_clause); // This is not necessary in DIP extended variables because we will immediately BT
      //                                 // and hence z will be unassigned

      //cout << z << " " << -x << endl;
      clause = {z,-x};   new_extended_definition_clause();
      //cout << z << " " << -y << endl;
      clause = {z,-y};   new_extended_definition_clause();
      clause.clear();
      
      assert(wtab.size() > 0);
      //cout << "Add " << z << " <-> " << x << " v " << y << " with clauses " << c1 << " " << c2 << " " << c3 << endl;
	    
    }
    else { // Variable has been created
      if (is_var_actively_deleted(z) or find_definition_clauses(z).size() == 0) { // Should create clauses

	
	cout << "Var " << z << " already existed but clauses removed " << endl;	
	cout << z << "(val " << int(val(z)) << ", lev " << var(z).level << ") <--> ";
	cout << x << "(val " << int(val(x)) << ", lev " << var(x).level << ")  OR ";
	cout << y << "(val " << int(val(y)) << ", lev " << var(y).level << ")" << endl;

	// if (max_var >= 493) {
	//   Clause* tmp = find_binary_extended_clause_in_DB(493,209);	
	//   if (tmp == NULL)
	//     cout << "Clause not found before crashing" << endl; 
	//   else cout << "Found before crashing (garbage " << tmp->garbage << ")" << endl;
	// }

	cout << endl << endl;
	for (auto c : clauses){
	  for (int l : *c) cout << l << " ";
	  cout << endl;
	}
	exit(1);
	if (!flags(z).active() and flags(z).unused()) mark_active(z);

	clause = {-z,x,y}; new_extended_definition_clause();
	clause = {z,-x};   new_extended_definition_clause();
	clause = {z,-y};   new_extended_definition_clause();
	clause.clear();

	mark_var_actively_deleted(z,false);
      }
      
    }
    return z;
  }

  bool Internal::is_extended_var (int v) {
    return abs(v) > orig_max_var and abs(v) <= max_var;
  }
  
  void Internal::notify_DIP (int x, int y) {
    er_manager.order_pair(x,y);
    ++dip2times[{x,y}];
  }

  int Internal::occurrences_DIP(int x, int y) {
    er_manager.order_pair(x,y);
    return dip2times[{x,y}];
  }

  bool Internal::removing_extended ( ) {
    return level == 0 and (max_var - orig_max_var - num_extended_vars_actively_deleted > lim.delete_ext);
    //return level == 0 and stats.conflicts >= lim.delete_ext;
  }

  void Internal::unwatch_and_mark_for_deletion (Clause *c) {
    if (not c->garbage) mark_garbage(c); // TODO: have a look at the treatment of redundant vs irredundant
    unwatch_clause(c);
  }
  
  void Internal::remove_extended ( ) {
    cout << "REMOVE EXTENDED (should improve to not remove the last ones, and keep a reasonable number of extended vars" << endl;
    cout << "Currently active should be " << max_var - orig_max_var - num_extended_vars_actively_deleted << endl;
    assert(level == 0);
    cout << "Use scores: " << use_scores() << endl;
    //lim.delete_ext += 300;
    lim.delete_ext += 2000;
    //lim.delete_ext += 100000;
    vector<pair<double,int>> candidates;
    for (int v = orig_max_var + 1; v <= max_var-30; ++v)  // do not remove last 30
      if (er_manager.is_deletable(v) and not is_var_actively_deleted(v)) candidates.push_back({score(v),v});
    cout << "Candidates " << candidates.size() << endl;
  // for (auto x : candidates) cout << x.second << " ";
    // cout << endl;
    
    sort(candidates.begin(),candidates.end());
    int first_kept = 0.3*candidates.size(); // remove candidates[0...first_kept-1]
    //if (candidates.size()) cout << "Largest score " << candidates.back().first << endl;

    // TO DO: if we are not using VSDIS but QUEUE then the scores are not very up to date.....
    
    int rem = 0;
    int skip = 0;
    for (int k = 0; k < first_kept; ++k) {
      int v = candidates[k].second; // Delete v
      vector<Clause*> cls = find_definition_clauses(v); // This is slow, but so far the only way to do it
                                                        // Since not many rounds of ext_var_delete are performed
                                                        // this might be a big problem
      if (cls.size() == 0) {
	//cout << "Definition for " << v << " not deleted because not found" << endl;
	mark_var_actively_deleted(v,true);
	++skip; continue;       
      } // do not remove if cannot find clauses (TO DO BE IMPROVED, BECAUSE THIS CAUSES THAT SOME VARS WILL NEVER BE DELETED, EVEN THOUGH THEY MIGHT BE USELESS)
      
      
      //      cout << "Remove definition for " << v << endl;
      
      unwatch_and_mark_for_deletion(cls[0]);
      unwatch_and_mark_for_deletion(cls[1]);
      unwatch_and_mark_for_deletion(cls[2]);
      er_manager.delete_definition(v);
      mark_var_actively_deleted(v,true);
      //cout << "Delete " << v << endl;
      ++rem;
    }

    cout << "Removed " << rem << " out of " << rem+skip << "(" << 100*double(rem)/(rem+skip) << " %)" << endl;
    cout << "Currently active is now " << max_var - orig_max_var - num_extended_vars_actively_deleted << endl;
    //if (rem > 0) exit(1);
  }

    Clause* Internal::find_binary_extended_clause (int x, int y) {
      //      cout << x << " and idx " << vlit(x) << " and size " << wtab.size() << endl;
      const auto end = watches(x).end();
      auto i = watches(x).begin();
      while (i != end) {
	if (i->binary() and i->blit == y and not i->clause->garbage) return i->clause;
	++i;
      }
      return NULL;
    }

    Clause* Internal::find_binary_extended_clause_in_DB (int x, int y) {
      //      cout << x << " and idx " << vlit(x) << " and size " << wtab.size() << endl;
      for (auto c : clauses) {
	if (c->size == 2) {
	  auto it = c->begin();
	  int lit1 = *it; ++it;
	  int lit2 = *it;
	  if ((lit1 == x and lit2 == y) or
	      (lit2 == x and lit1 == y)) return c;
	}	
      }
      return NULL;
    }

  
    Clause* Internal::find_ternary_extended_clause (int x, int y, int z) { // TO BE IMPROVED
      set<int> expected = {x,y,z};
      // cout << "Expected: ";
      // for (auto m : expected) cout << m << " ";
      // cout << endl;
      const auto endx = watches(x).end();      
      auto i = watches(x).begin();
      while (i != endx) {
	if (i->size == 3 and i->clause->extended and not i->clause->garbage) {
	  auto it = i->clause->begin();
	  set<int> found = {*it,*(it+1),*(it+2)};
	  if (found == expected) return i->clause;
	}
	++i;
      }

      //      cout << "First" << endl;
      const auto endy = watches(y).end();      
      i = watches(y).begin();
      while (i != endy) {
	if (i->size == 3 and i->clause->extended and not i->clause->garbage) {
	  auto it = i->clause->begin();
	  set<int> found = {*it,*(it+1),*(it+2)};
	  if (found == expected) return i->clause;
	}
	++i;
      }

      //      cout << "Second" << endl;
      const auto endz = watches(z).end();      
      i = watches(z).begin();
      while (i != endz) {
	if (i->size == 3 and i->clause->extended and not i->clause->garbage) {
	  auto it = i->clause->begin();
	  set<int> found = {*it,*(it+1),*(it+2)};
	  if (found == expected) return i->clause;
	}
	++i;
      }

      //      cout << "Final" << endl;
      // for (auto cla : clauses) {
      // 	for (auto it = cla->begin(); it != cla->end(); ++it)
      // 	  cout << *it << " ";
      // 	cout << endl;
      // }
      // exit(1);
      return NULL;
    }

  void Internal::write_all_clauses ( ) {
    for (auto cl : clauses){
      for (auto l : *cl) cout << l << " ";
      cout << "(deleted " << cl->garbage << " extended " << cl->extended << ")" << endl;
    }
  }
    
    vector<Clause*> Internal::find_definition_clauses(int z) {
      assert(z > 0);
      assert(is_extended_var(z));
      pair<int,int> def = er_manager.expand_definition(z);
      int x = def.first; int y = def.second; // z <--> x v y



      //cout << "Find definition " << z << "(val " << int(val(z)) << ") <--> " << x << "(val " << int(val(x)) << ") v " << y << "(val " << int(val(y)) << ")" << endl;

            
      // -x v z

      //cout << "Look for " << z << " " << -x << endl;	    
      Clause* c1 = find_binary_extended_clause(z,-x);
      if (c1 == NULL)  return {};
      //assert(c1 != NULL);
      // for (auto it = c1->begin(); it != c1->end(); ++it)
      // 	cout << *it << " ";
      // cout << endl;
      
      // -y v z
      //cout << "Look for " << z << " " << -y << endl;	    
      Clause* c2 = find_binary_extended_clause(z,-y);
      if (c2 == NULL) return {};
      //assert(c2 != NULL);
      // for (auto it = c2->begin(); it != c2->end(); ++it)
      // 	cout << *it << " ";
      // cout << endl;

      
      // -z v x v y
      //cout << "Look for " << -z << " " << x << " " << y << endl;
      Clause* c3 = find_ternary_extended_clause(-z,x,y);
      if (c3 == NULL) return {};
      //cout << c3 << endl;

      
      return {c1,c2,c3};
    }

  void Internal::mark_var_actively_deleted (int v, bool b){
    assert(v >= 0 and v <= max_var);
    if (b and not is_extended_actively_deleted[v]) ++num_extended_vars_actively_deleted;
    else if (not b and is_extended_actively_deleted[v]) --num_extended_vars_actively_deleted;    
    is_extended_actively_deleted[v] = b;
  }
  
  bool Internal::is_var_actively_deleted (int v) {
    assert(v >= 0 and v <= max_var);
    return is_extended_actively_deleted[v];
  }
}
