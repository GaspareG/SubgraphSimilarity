/*
  Author: Gaspare Ferraro
  Count the frequency of labels sequency of length k
  start/finish in every node using color-coding technique (parallel version)
  (by Alessio Conte)
*/

#include <vector>
#include <set>
#include <map>
#include <string>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <omp.h>
#include <getopt.h>
#include <unistd.h>

#ifdef K_8
#define MAXK 8
#define COLORSET uint8_t
#elif K_16
#define MAXK 16
#define COLORSET uint16_t
#elif K_64
#define MAXK 64
#define COLORSET uint64_t
#else
#define MAXK 32
#define COLORSET uint32_t
#endif

using namespace std;
typedef long long ll;

unsigned int N, M;
unsigned int k = 0, kp = 0;
int x = -1;
unsigned int thread_count = 0;
static int verbose_flag, help_flag;

ll cont = 0;
char *labels;
int *color;
vector<int> *G;

set<string> W;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Get pos-th bit in n
bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }

// Set pos-th bit in n
COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }

// Reset pos-th bit in n
COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }

// Random coloring graph using kp color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = rand() % kp;
}

// Dynamic programming processing
map< pair<COLORSET, string>, ll > *DP[MAXK+1];

void processDP() {

    if( verbose_flag ) printf("K = %u\n", 1);
    for(unsigned int u=0; u<N; u++)
      DP[1][u][ make_pair( setBit(0, color[u] ), string(&labels[u], 1) ) ] = 1;

    for(unsigned int i=2; i <= k; i++)
    {
      if( verbose_flag ) printf("K = %u\n", i);
      #pragma omp parallel for
      for(unsigned int u=0; u<N; u++)
      {
        for(int v : G[u])
        {
          for(auto d : DP[i-1][v])
          {
            COLORSET s = d.first.first;
            string l = d.first.second;
            ll f = d.second;
            if( getBit(s, color[u]) ) continue;

            COLORSET su = setBit(s, color[u]);
            string lu = string(l);
            lu += labels[u];
            ll fp = DP[i][u][ make_pair( setBit(s, color[u]), lu ) ];

            auto it = W.lower_bound(lu);
            if( it == W.end() ) continue;
            if( mismatch(lu.begin(), lu.end(), (*it).begin()).first != lu.end() ) continue;

            DP[i][u][make_pair(su, lu)] = f+fp;
          }
        }
      }
    }
}


void print_usage(char *filename) {
  printf(
      "Usage: ./%s -w filename -k length -K number -g filename -p threadcount -n index"
      "--help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-w, --sample filename\n");
  printf("\tInput file with sampled paths (default stdin)\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph in .nme.bin format (default stdin)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("-n, --node index\n");
  printf("\tIndex of the node to print sequence and frequency\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\tPrint status messages.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;

bool sampled_path_flag = false;
char *sampled_path = NULL;

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'k'},
      {"color", required_argument, 0, 'K'},
      {"input", required_argument, 0, 'g'},
      {"parallel", required_argument, 0, 'p'},
      {"sample", required_argument, 0, 'w'},
      {"node", required_argument, 0, 'n'},
      {"help", no_argument, &help_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  while (1) {
    c = getopt_long(argc, argv, "k:K:g:p:n:w:", long_options, &option_index);

    if (c == -1) break;

    switch (c) {
      case 'w':
        sampled_path_flag = true;
        if (optarg != NULL) sampled_path = optarg;
        break;
      case 'k':
        if (optarg != NULL) k = atoi(optarg);
        break;
      case 'K':
        if (optarg != NULL) kp = atoi(optarg);
        break;
      case 'n':
        if (optarg != NULL) x = atoi(optarg);
        break;
      case 'g':
        input_graph_flag = true;
        if (optarg != NULL) input_graph = optarg;
        break;
      case 'p':
        if (optarg != NULL) thread_count = atoi(optarg);
        break;
    }
  }

  if (help_flag || argc == 1) {
    print_usage(argv[0]);
    return 0;
  }

  if (k == 0) {
    printf("Invalid or missing path length value.\n");
    return 1;
  }

  if (k > MAXK || kp > MAXK) {
    printf("k or kp to high! (max value: %d\n", MAXK);
    return 1;
  }

  if (kp == 0) kp = k;

  if (thread_count > 0 && (int)thread_count < omp_get_max_threads()) {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
  }

  if (verbose_flag) {
    printf("Options:\n");
    printf("k = %d\n", k);
    printf("kp = %d\n", kp);
    printf("thread = %d\n", thread_count);
    printf("input_graph = %s\n", input_graph != NULL ? input_graph : "stdin");
  }
  if (verbose_flag) printf("Reading graph...\n");

  if (input_graph_flag) {
    if (input_graph == NULL) {
      printf("Input file name missing!\n");
      return 1;
    }

    int input_fd = open(input_graph, O_RDONLY, 0);
    if (input_fd == -1) {
      perror("Error opening input file");
      return 1;
    }
    read(input_fd, &N, sizeof(int));
    read(input_fd, &M, sizeof(int));

    color = new int[N];
    labels = new char[N];
    G = new vector<int>[N];

    read(input_fd, labels, N*sizeof(char));

    for(unsigned int i=0; i<N; i++)
      labels[i] += 'A';

    int ab[2];
    for (unsigned int i = 0; i < M; i++) {
      read(input_fd, ab, 2 * sizeof(int));
      G[ab[0]].push_back(ab[1]);
      G[ab[1]].push_back(ab[0]);
    }

  } else {
    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();

    color = new int[N];
    labels = new char[N];
    G = new vector<int>[N];

    for(unsigned int i = 0; i<N; i++)
      labels[i] = 'A'+nextInt();

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  if (verbose_flag) printf("Reading sampled paths (empty line to stop)\n");

  FILE *fd_w = stdin;
  if( sampled_path_flag && sampled_path != NULL )
  {
    if( (fd_w=fopen(sampled_path,"r")) == NULL)
    {
      perror("Error opening sampled paths file");
      return 1;
    }
  }

  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  while ( !feof(fd_w) && (read = getline(&line, &len, fd_w)) > 1 ) {
    line[read-1] = '\0';
    W.insert( string(line) );
  }

  if( verbose_flag ) printf("%zu string in W!\n", W.size());

  if( x >= (int) N )
  {
    x = -1;
    if( verbose_flag ) printf("Invalid node index, it will be ignored!\n");
  }

  // Create DP Table
  for (unsigned int i = 0; i <= k; i++)
    DP[i] = new map< pair<COLORSET, string>, ll >[N];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Processing DP
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("Finish processing DP table...\n");

  if( x != -1 )
  {
    ll cont = 0;
    printf("Tuple<ColorSet, String, Frequency> from node %d:\n", x);
    for(auto v : DP[k][x] )
    {
      cont += v.second;
      printf("\tC=%d S=%s F=%llu\n", v.first.first, v.first.second.c_str(), v.second );
    }
    if( verbose_flag ) printf("%llu k-labeled colorful-path", cont);
  }
  return 0;
}
