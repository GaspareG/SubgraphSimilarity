#include <bits/stdc++.h>
using namespace std;

map<int, set<int>> G;

int N, M;

int nextInt() {
  int tmp;
  scanf("%d", &tmp);
  return tmp;
}

void dfs(int node, set<int> *visited) {
  if (visited->find(node) == visited->end()) {
    visited->insert(node);
    for (int j : G[node]) dfs(j, visited);
  }
}

set<int> component(int node, set<int> *visited) {
  if (visited == NULL) visited = new set<int>();
  dfs(node, visited);
  return *visited;
}

bool connected() {
  if (G.size() == 0) return false;
  int first = G.begin()->first;
  return component(first, NULL).size() == G.size();
}

set<set<int>> components() {
//  printf("COMPONENT\n");

  set<set<int>> ret;

//  printf("COMPONENT\n");
  set<int> key;
  for (auto i : G) key.insert(i.first);
//  printf("COMPONENT 2\n");
  while (key.size()) {
  //  printf("COMPONENT 3\n");
    set<int> cc = component(*key.begin(), NULL);
  //    printf("COMPONENT 4\n");
      for( int ccc : cc ){ key.erase(ccc);
  //    printf("COMPONENT 5\n");
    ret.insert(cc); }
  //    printf("COMPONENT 6\n");
  }

  return ret;
}

void slashBurn() {
  int iteration = 0;
  int removedHubs = 0;

  printf("OK\n");
  vector<int> hubs;
  map<int, set<int>> replications;

    printf("OK\n");
  vector<set<int>> coms;

    printf("OK\n");
  vector<int> degOrd;
  for (int i = 0; i < N; i++) degOrd.push_back(i);
  sort(degOrd.begin(), degOrd.end(),
       [](int x, int y) { return G[x].size() < G[y].size(); });

         printf("OK\n");
  bool lastCommunity = false;

  while (degOrd.size()) {
    iteration++;
    printf("IT = %d || %d\n", iteration, degOrd.size());

    while (connected()) {
      if (degOrd.size() <= 2) {
        lastCommunity = true;
        break;
      }

      int maxN = degOrd[degOrd.size() - 1];
      degOrd.pop_back();

      removedHubs++;
      hubs.push_back(maxN);

      replications[maxN] = set<int>();

      for (int j : G[maxN]) {
        G[j].erase(maxN);
        replications[maxN].insert(j);
      }

      G.erase(maxN);

      sort(degOrd.begin(), degOrd.end(),
           [](int x, int y) { return G[x].size() < G[y].size(); });
    }
//    printf("DISCO\n");
    set<set<int>> ccs = components();
//    printf("OK\n");
    size_t maxSize = 0;
//    printf("OK\n");
    set<int> maxc;
//    printf("OK\n");
    for (auto cc : ccs) {
      if (cc.size() > maxSize) {
        maxSize = cc.size();
        maxc = cc;
      }
    }
//    printf("OK\n");
    ccs.erase(maxc);


//      printf("CCS OK\n");
    if (lastCommunity) {
      maxc.insert(hubs.begin(), hubs.end());
      coms.push_back(maxc);
    }

    for (auto cc : ccs) {
      set<int> icc = set<int>(cc.begin(), cc.end());

      for (int i : cc) {
        if (replications.find(i) != replications.end()) {
          icc.insert(replications[i].begin(), replications[i].end());
        }
      }
      coms.push_back(icc);
    }

    for (int i : degOrd)
      if (maxc.find(i) == maxc.end()) G.erase(i);

    for (pair<int, set<int>> g : G)
      for (int j : maxc) G[g.first].erase(j);

    vector<int> nDegOrd;
    for (int i : degOrd)
      if (maxc.find(i) == maxc.end()) nDegOrd.push_back(i);

    degOrd = nDegOrd;

    sort(degOrd.begin(), degOrd.end(),
         [](int x, int y) { return G[x].size() < G[y].size(); });
  }

//  coms.push_back(hubs);

}

int main(int argc, char **argv) {
  N = nextInt();
  M = nextInt();

  for (int i = 0; i < N; i++) G[i] = set<int>();

  for (int i = 0; i < M; i++) {
    int a = nextInt();
    int b = nextInt();
    G[a].insert(b);
    G[b].insert(a);
  }

  slashBurn();
}
