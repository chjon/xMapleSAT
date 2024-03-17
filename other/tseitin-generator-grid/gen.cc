#include <iostream>
#include <vector>
#include <map>
#include <assert.h>

using namespace std;

int n;
map<pair<pair<int,int>,pair<int,int>>, int> edge2var;
vector<pair<pair<int,int>,pair<int,int>>> var2edge;
int nextVar = 1;
vector<vector<int>> clauses;

vector<pair<int,int>> dirs = {{0,1}, {0,-1}, {1,0}, {-1,0}};

int encodeEdge ( pair<int,int> v1, pair<int,int> v2) {

  //cout << "Orig era (" << v1.first << "," << v1.second << ") -- (" << v2.first << "," << v2.second << ")" << endl;
  
  if (v1.first > v2.first) swap(v1,v2);
  else if (v1.first == v2.first and v1.second > v2.second) swap(v1,v2);

  // cout << "Despres es (" << v1.first << "," << v1.second << ") -- (" << v2.first << "," << v2.second << ")" << endl;
  // cout << "Encode edge (" << v1.first << "," << v1.second << ") -- (" << v2.first << "," << v2.second << ")" << endl;
  if (edge2var.count({v1,v2}) == 1) {
    //    cout << "Retorno " << edge2var[{v1,v2}] << endl << endl;
    return edge2var[{v1,v2}];
  }
  else {
    edge2var[{v1,v2}] = nextVar;
    var2edge.push_back({v1,v2});      
    ++nextVar;
    //cout << "Retorno " << nextVar - 1 << endl << endl;
    return nextVar-1;
  }
}


vector<pair<int,int>> neighbours (int i, int j) {
  vector<pair<int,int>> res;
  for (auto d : dirs) {
    int n_i = i + d.first;
    int n_j = j + d.second;
    if (n_i < 0) n_i += n;
    if (n_j < 0) n_j += n;
    if (n_i >= n) n_i -= n;
    if (n_j >= n) n_j -= n;
    res.push_back({n_i,n_j});
  }
  return res;
}

vector<vector<int>> forbidEven = {
  {0,0,0,0},
  {1,1,0,0},
  {1,0,1,0},
  {1,0,0,1},
  {0,1,1,0},
  {0,1,0,1},
  {0,0,1,1},
  {1,1,1,1}
};

vector<vector<int>> forbidOdd = {
  {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {0,0,0,1},
  {0,1,1,1},
  {1,0,1,1},
  {1,1,0,1},
  {1,1,1,0}
};




void encode (const pair<int,int>& v, const vector<pair<int,int>>& neigh, int c) {
  assert(neigh.size() == 4);

  // cout << "Encode vertex (" << v.first << "," << v.second << ") with neighbours ";
  // for (auto p :  neigh) cout << "(" << p.first << "," << p.second << ") ";
  // cout << " and charge " << c << endl;
  // cout << "Edges are: ";
  // for (int k = 0; k < 4; ++k) {
  //   cout << "(" << v.first << "," << v.second << ") -- (" << neigh[k].first << "," << neigh[k].second << ") ===> " << encodeEdge(v,neigh[k]) << endl;
  // }
  
  vector<vector<int>>& pattern = (c % 2 == 1 ? forbidEven : forbidOdd);
  for (auto p : pattern) {
    clauses.push_back({});
    for (int k = 0; k < 4; ++k){
      if (p[k] == 0)
	clauses.back().push_back(-encodeEdge( v ,  neigh[k] ));
      else
	clauses.back().push_back( encodeEdge( v ,  neigh[k] ));				   
    }
    // for (auto x : clauses.back()) cout << x << " ";
    // cout << 0 << endl;
  }
}

int main ( ){
  cin >> n;

  int nVertices = n*n; // verticesvariables
  vector<int> charges(nVertices,0);
  charges[0] = 1;
  var2edge.push_back({{0,0},{0,0}}); // fake
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      vector<pair<int,int>> neigh = neighbours(i,j);
      encode({i,j},neigh,charges[i*n+j]);
    }    
  }

  cout << "p cnf " << nextVar - 1 << " " << clauses.size() << endl;
  for (auto c : clauses) {
    for (auto l : c) cout << l << " ";
    cout << 0 << endl;
  }

  
  for (int v = 0; v < var2edge.size(); ++v) {
    cout << "{" << var2edge[v].first.first << "," << var2edge[v].first.second << ",";
    cout << var2edge[v].second.first << "," << var2edge[v].second.second << "}," << endl;
  }

}
