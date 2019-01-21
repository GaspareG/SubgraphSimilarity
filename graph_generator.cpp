#include <bits/stdc++.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

using namespace std;

set<int> *G;
int *labels;

int main(int argc, char **argv) {
  int seed = 42;
  if (argc < 3) {
    printf("Usage: %s N M output label A B seed\n", argv[0]);
    printf("\tN, number of nodes\n");
    printf("\tM, number of edges\n");
    printf("\toutput, filename (optional, default: stdout)\n");
    printf("\tA, size of set A (jaccard)\n");
    printf("\tB, size of set B (jaccard)\n");
    printf(
        "\tlabel, labeled graph in [0, label-1] (optional, 0 if not "
        "labeled)\n");
    printf("\tseed, (optional)\n");
    return 1;
  }
  bool jaccard = false;
  int N = atol(argv[1]);
  int M = atol(argv[2]);
  int Mc = M;
  int fd = -1;
  int label = 0;
  int sampleS[2] = {0, 0};
  int *sample[2];

  if (argc >= 4 && strcmp(argv[3],"stdout") != 0) {
    fd = open(argv[3], O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
    if (fd == -1) {
      perror("Opening file");
      return -1;
    }
  }
  if (argc >= 5) label = atoi(argv[4]);

  if( argc >= 7 )
  {
    jaccard = true;
    sampleS[0] = atoi(argv[5]);
    sampleS[1] = atoi(argv[6]);
  }

  if (argc >= 8) seed = atoi(argv[7]);

  assert(seed >= 0);
  assert(N > 0);
  assert(M > 0);
  assert(M >= N);
  assert(sampleS[0] >= 0);
  assert(sampleS[0] < N);
  assert(sampleS[1] >= 0);
  assert(sampleS[1] < N);
  assert(label >= 0);
  assert(label < 255);
  assert((long long)M <= (long long)N * (N - 1) / 2);

  srand(seed);

  G = new set<int>[N];
  labels = new int[N];

  for (int i = 0; i < N - 1; i++) {
    G[i].insert(i + 1);
    G[i + 1].insert(i);
    M--;
  }


  if (label > 0)
    for (int i = 0; i < N; i++) labels[i] = rand() % label;

  int *toWrite = new int[2 + Mc * 2 + (jaccard ? (2+sampleS[0]+sampleS[1]) : 0 )];
  int idx = 0;
  do {
    int i = rand() % N;
    if (G[i].size() == (size_t)N - 1) continue;
    while (1) {
      int j = rand() % N;
      if (j == i || G[i].find(j) != G[i].end()) continue;
      G[i].insert(j);
      G[j].insert(i);
      break;
    }
    M--;
  } while (M);


  toWrite[idx++] = N;
  toWrite[idx++] = Mc;
  for (int i = 0; i < N; i++)
    for (auto j : G[i])
      if (i < j) {
        toWrite[idx++] = i;
        toWrite[idx++] = j;
      }

  if(jaccard)
  {
    sample[0] = new int[sampleS[0]];
    sample[1] = new int[sampleS[1]];
    vector<int> V;
    for(int i=0; i<N; i++) V.push_back(i);
    random_shuffle(V.begin(), V.end());
    for(int i=0; i<sampleS[0]; i++)
      sample[0][i] = V[i];
    for(int i=0; i<sampleS[0]; i++)
      sample[1][i] = V[N-1-i];
  }

  if (fd != -1) {
    write(fd, (void *)toWrite, 2 * sizeof(int));
    if (label > 0) write(fd, (void *)labels, N * sizeof(int));
    write(fd, (void *)&toWrite[2], 2 * Mc * sizeof(int));
    if(jaccard)
    {
      write(fd, (void *) sampleS, 2*sizeof(int));
      write(fd, (void *) sample[0], sampleS[0]*sizeof(int));
      write(fd, (void *) sample[1], sampleS[1]*sizeof(int));
    }

  } else {
    printf("%d %d\n", toWrite[0], toWrite[1]);
    if (label > 0) {
      for (int i = 0; i < N; i++) printf("%d ", labels[i]);
      printf("\n");
    }
    for (int i = 1; i < Mc + 1; i++)
      printf("%d %d\n", toWrite[2 * i], toWrite[2 * i + 1]);

    if( jaccard )
    {
      printf("%d %d\n", sampleS[0], sampleS[1]);
      for(int i=0; i<sampleS[0]; i++) printf("%d ", sample[0][i]);
      printf("\n");
      for(int i=0; i<sampleS[1]; i++) printf("%d ", sample[1][i]);
      printf("\n");
    }
  }
  close(fd);
  return 0;
}
