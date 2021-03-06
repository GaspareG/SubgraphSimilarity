/*
  Autore: Gaspare Ferraro
  Conta i k-path indotti in un grafo non orientato
  con una dfs
*/
#include <bits/stdc++.h>

using namespace std;

int N, M, k;
set<int> *G;

int nextInt() {
  int tmp;
  scanf("%d", &tmp);
  return tmp;
}

bool *in;
vector<int> path;

long long cont = 0;
void dfs(int a) {
  bool induced = true;
  for (int i = 0; i < (int)path.size() - 1 && induced; i++) {
    if (G[path[i]].find(a) != G[path[i]].end()) induced = false;
  }
  if (!induced) return;
  path.push_back(a);
  in[a] = true;
  if (path.size() == (size_t)k) {
    cont++;
    // for(int i : path)
    //   printf("%d ", i);
    // printf("\n");
  } else {
    for (int b : G[a]) {
      if (in[b]) continue;
      dfs(b);
    }
  }
  in[a] = false;
  path.pop_back();
}

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Usage: ./%s k\n", argv[0]);
  }

  k = atol(argv[1]);
  N = nextInt();
  M = nextInt();

  G = new set<int>[N];
  in = new bool[N];

  for (int i = 0; i < M; i++) {
    int a = nextInt();
    int b = nextInt();
    G[a].insert(b);
    G[b].insert(a);
  }

  for (int i = 0; i < N; i++) dfs(i);
  printf("%llu\n", cont);
  return 0;
}
