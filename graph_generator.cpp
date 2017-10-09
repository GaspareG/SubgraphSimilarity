#include <bits/stdc++.h>
using namespace std;

set<int> *G;

int main(int argc, char **argv) {
  srand(42);
  if (argc < 3) {
    printf("Usage: %s N M\n", argv[0]);
    return 1;
  }
  int N = atol(argv[1]);
  int M = atol(argv[2]);

  printf("%d %d\n", N, M);

  G = new set<int>[N];

  while (M) {
    int i = rand() % N;
    int j = rand() % N;
    if (i > j) swap(i, j);
    if (i != j && G[i].find(j) == G[i].end()) {
      G[i].insert(j);
      M--;
    }
  }

  for (int i = 0; i < N; i++)
    for (auto j : G[i]) printf("%d %d\n", i, j);

  return 0;
}
