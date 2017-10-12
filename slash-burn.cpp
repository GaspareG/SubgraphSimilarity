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
  set<set<int>> ret;

  set<int> key;
  for (auto i : G) key.insert(i.first);

  while (key.size()) {
    set<int> cc = component(*key.begin(), NULL);
    key.erase(cc.begin(), cc.end());
    ret.insert(cc);
  }

  return ret;
}

void slashBurn() {
  int iteration = 0;
  int removedHubs = 0;

  vector<int> hubs;
  map<int, set<int>> replications;

  vector<set<int>> coms;

  vector<int> degOrd;
  for (int i = 0; i < N; i++) degOrd.push_back(i);
  sort(degOrd.begin(), degOrd.end(),
       [](int x, int y) { return G[x].size() < G[y].size(); });

  bool lastCommunity = false;

  while (degOrd.size()) {
    iteration++;
    printf("IT = %d\n", iteration);

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

    set<set<int>> ccs = components();
    size_t maxSize = 0;
    set<int> maxc;
    for (auto cc : ccs) {
      if (cc.size() > maxSize) {
        maxSize = cc.size();
        maxc = cc;
      }
    }
    ccs.erase(maxc);

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
