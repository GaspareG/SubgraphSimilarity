/*
  Author: Gaspare Ferraro
  Count simple k-path in a graph using
  the divide-and-color technique
*/
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

int N, M, k, kp;
set<int> *G;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

set<pair<int, int> > list_k_path(unordered_set<int> Gp, int k) {
  printf("%d [%zu]\n", k, Gp.size());
  set<pair<int, int> > ret;

  // Caso base
  if (k == 1) {
    for (int v : Gp) ret.insert(make_pair(v, v));
  } else {
    // Ripeto per 3*2^k come consigliato nell'articolo originale
    for (int t = 0; t < 3 * (1 << k); t++) {
      // Partiziono G nei due sottografi G1 e G2
      unordered_set<int> Gl[2];
      for (int v : Gp) Gl[rand() % 2].insert(v);

      // Ricorro nei due sottografi cercando path lunghi k/2
      auto L1 = list_k_path(Gl[0], k - (int)(k / 2));
      auto L2 = list_k_path(Gl[1], (int)(k / 2));

      /*
      Versione O(n^4)
      for(int u : Gp)
      for(int v : Gp)
      for(int w : Gp)
      for(int x : Gp)
      {
        if( L1.find(make_pair(u,v)) == L1.end() ) continue;
        if( L2.find(make_pair(w,x)) == L2.end() ) continue;
        if( G[v].find(w) == G[v].end() ) continue;
        ret.insert( make_pair(u, x) );
      }
      */

      // Versione O( |L1| * |L2| )
      for (auto l1 : L1)  // l1 = (u,v)
      {
        int u = l1.first;
        int v = l1.second;
        for (auto l2 : L2)  // l2 = (w,x)
        {
          // Controllo se (u,w) € E
          int w = l2.first;
          int x = l2.second;
          if (G[v].find(w) != G[v].end()) ret.insert(make_pair(u, x));
        }
      }
    }
  }
  return ret;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Usage: %s k\n", argv[0]);
    return 1;
  }
  srand(42);
  N = nextInt();
  M = nextInt();
  k = atol(argv[1]);

  G = new set<int>[N];

  for (int i = 0; i < M; i++) {
    int a = nextInt();
    int b = nextInt();
    G[a].insert(b);
    G[b].insert(a);
  }

  unordered_set<int> Gp;
  for (int i = 0; i < N; i++) Gp.insert(i);

  auto result = list_k_path(Gp, k);

  printf("%zu\n", result.size());
  //    for(auto r : result)
  //      printf("%d - %d\n", r.first, r.second);
}
