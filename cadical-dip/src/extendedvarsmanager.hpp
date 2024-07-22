#ifndef _extendedvarsmanager_hpp_INCLUDED
#define _extendedvarsmanager_hpp_INCLUDED

#include <map>

// We cannot delete an extended variable x if it used in another definition
// For that purpose, we can construct a graph such that for every definition x <-> y v z
// we add two edges x->y and x->z. We will be able to delete a variable if it is extended
// and has indegree zero

namespace CaDiCaL {

  struct ExtendedVarsManager {

    map<pair<int,int>,int> pair2lit;
    vector<pair<int,int>> var2pair;
    
    vector<int> indegree; // indegree[vlit(l)] is the indegree of literal l
    
    int max_ext_var;

    unsigned lit2idx (int lit) { return (lit < 0) + 2u * (unsigned) abs (lit); }
    unsigned int max_lit_idx ( ) { return lit2idx(-max_ext_var);}

    bool part_of_definition (int l) {
      if (lit2idx(l) >= indegree.size()) return false;
      //      assert(lit2idx(l) < indegree.size());
      return indegree[lit2idx(l)] > 0;
    }
    
    void order_pair (int& x, int& y) {
      assert(abs(x) != abs(y));
      if (abs(x) > abs(y)) swap(x,y);
    }

    bool exists_extended_var (int x, int y) {
      order_pair(x,y);
      auto it = pair2lit.find({x,y});
      return it != pair2lit.end();
    }


    int find_definition (int x, int y) {
      if (abs(x) == abs(y)) {cout << "ERROR: " << x << " " << y << endl; exit(1);}
      order_pair(x,y);
      auto it = pair2lit.find({x,y});
      if (it == pair2lit.end()) return 0;
      else return it->second;
    }
    
    pair<int,int> expand_definition (int v) {
      assert( v > 0);
      assert(v < int(var2pair.size()));
      return var2pair[v];
    }

    int add_extended_var(int x, int y, int z) {
      order_pair(x,y);
      auto it = pair2lit.find({x,y});
      if (it != pair2lit.end()) return it->second;
      else {
	max_ext_var = max(max_ext_var, z);
	if (indegree.size() <= max_lit_idx()) {
	  indegree.resize( max_lit_idx() + 1, 0);
	  var2pair.resize( max_ext_var + 1, {0,0});
	}
	pair2lit.insert({{x,y},z});
	var2pair[z] = {x,y};
	
	// We are adding z <-> x v y. That is, edges z->x and z->y
	assert(lit2idx(x) < indegree.size());
	assert(lit2idx(y) < indegree.size());
	++indegree[lit2idx(x)];
	++indegree[lit2idx(y)];
	//	cout << "Adding " << z << " <--> " << x << " v " << y << endl;
	// cout << "in(" << x << ") = " << indegree[lit2idx(x)] << endl;
	// cout << "in(" << y << ") = " << indegree[lit2idx(y)] << endl;
	// cout << z << " is extended " << is_extended[abs(z)] << endl;
	// cout << endl;
	return z;
      }      
    }

    bool is_deletable (int v) {
      assert(v > 0);
      assert(lit2idx(v) < indegree.size());  assert(lit2idx(-v) < indegree.size());
      return indegree[lit2idx(v)] == 0 and indegree[lit2idx(-v)] == 0;
    }

    void delete_definition (int v) {
      //      cout << "Delete definition " << v << endl;
      assert(v > 0);
      assert(v < int(var2pair.size()));
      pair<int,int> def = var2pair[v];
      assert(def.first and def.second);
      var2pair[v] = {0,0};
      pair2lit.erase(def);
      int x = def.first, y = def.second;
      assert(indegree[lit2idx(x)] > 0); assert(indegree[lit2idx(y)] > 0);
      assert(lit2idx(x) < indegree.size());
      assert(lit2idx(y) < indegree.size());
      --indegree[lit2idx(x)]; --indegree[lit2idx(y)];
      
      // we do not modify max_ext_var. This is not too precise, but since we are not reusing vars this is ok
    }
    
    ExtendedVarsManager ( ) : max_ext_var(0) { }

  };
  
} // namespace CaDiCaL

#endif
