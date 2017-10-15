#include <bits/stdc++.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

using namespace std;

set<int> *G;

int main(int argc, char **argv) {
  int seed = 42;
  if (argc < 3) {
    printf("Usage: %s N M output seed\n", argv[0]);
    printf("\tN, number of nodes\n");
    printf("\tM, number of edges\n");
    printf("\toutput, filename (optional, default: stdout)\n");
    printf("\tseed, (optional)\n");
    return 1;
  }
  int N = atol(argv[1]);
  int M = atol(argv[2]);
  int Mc = M;
  int fd = -1;
  if( argc >= 4 )
  {
    fd = open(argv[3], O_CREAT | O_WRONLY, S_IRUSR | S_IWUSR);
    if( fd == -1 )
    {
      perror("Opening file");
      return -1;
    }
  }
  if( argc >= 5 
)
    seed = atoi(argv[4]);

  assert(seed >= 0);
  assert(N > 0);
  assert(M > 0);
  assert(M >= N);
  assert((long long) M <= (long long) N*(N-1)/2 );

  srand(seed);

  G = new set<int>[N];
  
  for(int i=0; i<N-1; i++)
  {
    G[i].insert(i+1);
    G[i+1].insert(i);
    M--;
  }
  
  int *toWrite = new int[2+Mc*2];
  int idx = 0;
  do
  {
    int i = rand() % N;
    if( G[i].size() == (size_t)N-1 ) continue;
    while(1)
    {
      int j = rand() % N ;
      if( j == i || G[i].find(j) != G[i].end() ) continue;
      G[i].insert(j);
      G[j].insert(i);
      break;
    }
    M--;
  }
  while( M );

  toWrite[idx++] = N;
  toWrite[idx++] = Mc;
  for (int i = 0; i < N; i++)
    for (auto j : G[i])
      if( i < j )
      {
        toWrite[idx++] = i;
        toWrite[idx++] = j;
      }

  if( fd != -1) 
    write(fd, (void*) toWrite, (2+Mc*2)*sizeof(int) );
  else
    for(int i=0; i<Mc+1; i++)
      printf("%d %d\n",toWrite[2*i], toWrite[2*i+1]);
  return 0;
}
