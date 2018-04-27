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
#include <algorithm>
#include <iterator>
#include <random>
#include <limits>
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
unsigned int k = 0, kp = 0, r = 0;
unsigned int thread_count = 0;
static int verbose_flag, help_flag;

ll cont = 0;
char *labels;
int *color;
vector<int> *G;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Random generator
random_device rd;
mt19937_64 eng(rd());
uniform_int_distribution<unsigned long long> distr;

// Get pos-th bit in n
bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }

// Set pos-th bit in n
COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }

// Reset pos-th bit in n
COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }

// Complementary set of a COLORSET
COLORSET getCompl(COLORSET n) { return ((1 << kp) - 1) & (~n); }

// Random coloring graph using kp color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = rand() % kp;
}

// Dynamic programming processing
map<COLORSET, ll> *DP[MAXK + 1];

void processDP() {
  if (verbose_flag) printf("K = %u\n", 1);
  for (unsigned int u = 0; u < N; u++) DP[1][u][setBit(0, color[u])] = 1;

  for (unsigned int i = 2; i <= k; i++) {
    if (verbose_flag) printf("K = %u\n", i);
#pragma omp parallel for
    for (unsigned int u = 0; u < N; u++) {
      for (int v : G[u]) {
        for (auto d : DP[i - 1][v]) {
          COLORSET s = d.first;
          ll f = d.second;

          if (getBit(s, color[u])) continue;

          ll fp = DP[i][u][setBit(s, color[u])];
          DP[i][u][setBit(s, color[u])] = f + fp;
        }
      }
    }
  }
}

// Sample r(*2) path from DP matrix
set<string> sample() {
  set<string> out;
  while (r) {
    string pi;
    COLORSET D = 0;
    COLORSET kD = getCompl(D);
    int u = N;

    for (unsigned int i = k; i > 0; i--) {
      vector<int> A;
      vector<ll> freqA;
      ll last = 0ll;
      ll sum = 0ll;
      for (int v : G[u]) {
        auto it = DP[i][v].find(kD);
        if (it != DP[i][v].end()) {
          A.push_back(v);
          freqA.push_back(last);
          last = (*it).second;
          sum += last;
        }
      }

      if (A.size() == 0) break;
      ll rndIdx = distr(eng) % sum;
      int v = A[distance(freqA.begin(), upper_bound(freqA.begin(), freqA.end(), rndIdx) - 1)];

      u = v;
      pi += labels[v];
      D = setBit(D, color[v]);
      kD = getCompl(D);
    }

    if (pi.size() < k) continue;

    if (out.find(pi) == out.end()) {
      out.insert(pi);
      reverse(pi.begin(), pi.end());
      out.insert(pi);
      r--;
    }
  }
  return out;
}

void print_usage(char *filename) {
  printf(
      "Usage: ./%s -k length -K number -g filename -p threadcount -r number -R "
      "filename --help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph in .nme.bin format (default stdin)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("-r, --paths number\n");
  printf("\tNumber of paths to be sampled\n");

  printf("-R, --sample filename\n");
  printf("\tOutputfile for sampled paths (default stdout)\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\tPrint status messages.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;

bool output_sample_flag = false;
char *output_sample = NULL;

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'k'},
      {"color", required_argument, 0, 'K'},
      {"input", required_argument, 0, 'g'},
      {"parallel", required_argument, 0, 'p'},
      {"paths", required_argument, 0, 'r'},
      {"sample", required_argument, 0, 'R'},
      {"help", no_argument, &help_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  while (1) {
    c = getopt_long(argc, argv, "k:K:g:p:n:r:R:", long_options, &option_index);

    if (c == -1) break;

    switch (c) {
      case 'k':
        if (optarg != NULL) k = atoi(optarg);
        break;
      case 'K':
        if (optarg != NULL) kp = atoi(optarg);
        break;
      case 'g':
        input_graph_flag = true;
        if (optarg != NULL) input_graph = optarg;
        break;
      case 'r':
        if (optarg != NULL) r = atoi(optarg);
        break;
      case 'R':
        output_sample_flag = true;
        if (optarg != NULL) output_sample = optarg;
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
    G = new vector<int>[N + 1];

    read(input_fd, labels, N * sizeof(char));

    for (unsigned int i = 0; i < N; i++) labels[i] += 'A';

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
    G = new vector<int>[N + 1];

    for (unsigned int i = 0; i < N; i++) labels[i] = 'A' + nextInt();

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  for (unsigned i = 0; i < N; i++) G[N].push_back(i);

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  // Create DP Table
  for (unsigned int i = 0; i <= k; i++) DP[i] = new map<COLORSET, ll>[N];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Processing DP
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("Finish processing DP table...\n");

  if (r > 0) {
    if (verbose_flag) printf("Sampling paths...\n");
    set<string> sampledPath = sample();
    if (verbose_flag) printf("Finish sampling path...\n");

    if (output_sample_flag && output_sample != NULL) {
      FILE *fdout;
      if ((fdout = fopen(output_sample, "w")) == NULL) {
        perror("Error opening output sample file");
        return 1;
      }
      for (string s : sampledPath) fprintf(fdout, "%s\n", s.c_str());
    } else
      for (string s : sampledPath) printf("%s\n", s.c_str());
  }
  return 0;
}
