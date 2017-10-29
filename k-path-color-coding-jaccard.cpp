/*
  Author: Gaspare Ferraro
  Count and find simple k colorful-path in a graph
  using the color-coding technique (parallel version)
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
unsigned k = 0, kp = 0;
unsigned thread_count = 0;
static int verbose_flag, help_flag;

ll cont = 0;
int *color;
char *label;
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

// Path label
string L(vector<int> P)
{
    string l;
    l.reserve(P.size());
    for(size_t i=0; i<P.size(); i++) l[i] = label[P[i]];
    return l;
}

// Link
map<pair<int, COLORSET>, vector<int>> links;

inline void addLink(int x, COLORSET C, int j) {
  auto key = make_pair(x, C);
  if (links.find(key) == links.end()) links[key] = vector<int>();
  links[key].push_back(j);
}

// Oracle
vector<int> H(int x, COLORSET C) { return links[make_pair(x, C)]; }

void list_k_path(FILE *out, vector<int> ps, COLORSET cs, int x) {
  vector<int> oracle = H(x, cs);

  if (ps.size() + 1 == k) {
    cont++;
    for (int j : ps) fprintf(out, "%d ", j);
    fprintf(out, "\n");
  } else
    for (int v : oracle) {
      ps.push_back(v);
      list_k_path(out, ps, setBit(cs, color[v]), v);
      ps.pop_back();
    }
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

// Random sampling
vector<int> randomColorfulSample(vector<int> X)
{
  set<string> W;
  set< vector<int> > R;
  vector<int> AB;
  vector<ll> freqAB;
  for(int a : A) AB.push_back(a);
  for(int b : B) AB.push_back(b);
  sort(AB.begin(), AB.end() );
  AB.erase( unique( AB.begin(), AB.end() ), AB.end() );
  ll last = 0ll;
  ll sum = 0ll;
  for(int x : AB)
  {
    ll freq = DP[k][x][getCompl(0)];
    ll alpha = 0;
    if( A.find(x) != A.end() ) alpha++;
    if( B.find(x) != B.end() ) alpha++;
    freqAB.push_back(last);
    last = freq*alpha;
    sum += last;
  }

  while( r )
  {
    ll rndIdx = distr(eng) % sum;
    int u = AB[distance(freqAB.begin(), upper_bound(freqAB.begin(), freqAB.end(), rndIdx) - 1)];

    vector<int> P = randomPathFrom(u);
    if( R.find(P) == R.end() )
    {
      R.insert(P);
      r--;
    }
  }

  for(vector<int> r : R)  W.insert( L(r) );

  return W;
}


set<string> BCSampler(set<int> A, set<int> B, int r)
{
  set<string> W;
  set< vector<int> > R;
  vector<int> AB;
  vector<ll> freqAB;
  for(int a : A) AB.push_back(a);
  for(int b : B) AB.push_back(b);
  sort(AB.begin(), AB.end() );
  AB.erase( unique( AB.begin(), AB.end() ), AB.end() );
  ll last = 0ll;
  ll sum = 0ll;
  for(int x : AB)
  {
    ll freq = DP[k][x][getCompl(0)];
    ll alpha = 0;
    if( A.find(x) != A.end() ) alpha++;
    if( B.find(x) != B.end() ) alpha++;
    freqAB.push_back(last);
    last = freq*alpha;
    sum += last;
  }

  while( r )
  {
    ll rndIdx = distr(eng) % sum;
    int u = AB[distance(freqAB.begin(), upper_bound(freqAB.begin(), freqAB.end(), rndIdx) - 1)];

    vector<int> P = randomPathFrom(u);
    if( R.find(P) == R.end() )
    {
      R.insert(P);
      r--;
    }
  }

  for(vector<int> r : R)  W.insert( L(r) );

  return W;
}

void print_usage(char *filename) {
  //  printf("Usage: ./%s -k length -K number -g filename -f format -t filename
  //  -T filename -p threadcount -h -v\n",filename);
  printf(
      "Usage: ./%s -k length -K number -g filename -p threadcount "
      "--help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph (default stdin)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\nPrint status messages.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;

bool list_path_flag = false;
char *list_path = NULL;

bool format_name_flag = false;
char *format_name = NULL;

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'k'},
      {"color", required_argument, 0, 'K'},
      {"input", required_argument, 0, 'g'},
      {"format", required_argument, 0, 'f'},
      {"list", required_argument, 0, 'l'},
      {"parallel", required_argument, 0, 'p'},
      {"help", no_argument, &help_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  while (1) {
    c = getopt_long(argc, argv, "k:K:g:f:l:p:", long_options, &option_index);

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
      case 'f':
        format_name_flag = true;
        if (optarg != NULL) format_name = optarg;
        break;
      case 'l':
        list_path_flag = true;
        if (optarg != NULL) list_path = optarg;
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
    printf("format_name = %s\n", format_name != NULL ? format_name : "stdin");
    printf("list_path   = %s\n", list_path != NULL ? list_path : "stdin");
  }
  if (verbose_flag) printf("Reading graph...\n");

  if (input_graph_flag) {
    if (input_graph == NULL) {
      printf("Input file name missing!\n");
      return 1;
    }
    if (!format_name_flag || format_name == NULL) {
      printf("Input format missing!\n");
      return 1;
    }
    if (strcmp(format_name, "snap") == 0) {
      FILE *input_fd = fopen(input_graph, "r");
      if (input_fd == NULL) {
        perror("Error opening input file");
        return 1;
      }
      set<pair<int, int> > edge;
      char *buffer = NULL;
      size_t n;
      int ab[2];
      N = 0;
      do {
        if( buffer != NULL ) free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if (line_length == 0) continue;
        if (buffer[0] == '#') continue;
        sscanf(buffer, "%d %d", ab, ab + 1);
        edge.insert(make_pair(ab[0], ab[1]));
        N = (unsigned)ab[0] > N ? ab[0] : N;
        N = (unsigned)ab[1] > N ? ab[1] : N;
      } while (!feof(input_fd));

      M = edge.size();
      color = new int[N + 1];
      G = new vector<int>[N + 1];
      for (auto e : edge) {
        G[e.first].push_back(e.second);
        G[e.second].push_back(e.first);
      }
    } else if (strcmp(format_name, "nde") == 0) {
      FILE *input_fd = fopen(input_graph, "r");
      if (input_fd == NULL) {
        perror("Error opening input file");
        return 1;
      }
      char *buffer;
      size_t n;

      do {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if (line_length == 0) continue;
        if (buffer[0] == '#') continue;
        break;
      } while (!feof(input_fd));
      sscanf(buffer, "%u", &N);
      int Nrim = N;

      M = 0;
      int ab[2];
      while (Nrim > 0 && !feof(input_fd)) {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if (line_length == 0) continue;
        if (buffer[0] == '#') continue;
        sscanf(buffer, "%d %d", ab, ab + 1);
        Nrim--;
        M += ab[1];
        N = (unsigned)ab[0] > N ? (unsigned)ab[0] : N;
      }

      color = new int[N + 1];
      G = new vector<int>[N + 1];

      int Mrim = M;
      while (Mrim > 0 && !feof(input_fd)) {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if (line_length == 0) continue;
        if (buffer[0] == '#') continue;
        sscanf(buffer, "%d %d", ab, ab + 1);
        G[ab[0]].push_back(ab[1]);
        G[ab[1]].push_back(ab[0]);
      }

    } else if (strcmp(format_name, "nme") == 0) {
      int input_fd = open(input_graph, O_RDONLY, 0);
      if (input_fd == -1) {
        perror("Error opening input file");
        return 1;
      }
      read(input_fd, &N, sizeof(int));
      read(input_fd, &M, sizeof(int));

      color = new int[N + 1];
      G = new vector<int>[N + 1];
      int ab[2];
      for (unsigned int i = 0; i < M; i++) {
        read(input_fd, ab, 2 * sizeof(int));
        G[ab[0]].push_back(ab[1]);
        G[ab[1]].push_back(ab[0]);
      }
    } else {
      printf("Wrong input format (only 'snap', 'nde' or 'nme'\n");
      return 1;
    }
  } else {
    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();

    color = new int[N + 1];
    G = new vector<int>[N + 1];

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  for (unsigned int i = 0; i < N; i++) {
    G[N].push_back(i);
    G[i].push_back(N);
  }

  // Create DP Table
  for (unsigned int i = 0; i <= k + 1; i++)
    DP[i] = new map<COLORSET, ll>[N + 1];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();
  color[N] = kp;
  k++;

  // Fill dynamic programming table
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("End processing DP table...\n");


  // Count ad list k-colorful path
  // if (list_path_flag) {
  //   FILE *list_fd = stdout;
  //   if (list_path != NULL) {
  //     list_fd = fopen(list_path, "w");
  //     if (list_fd == NULL) {
  //       perror("Error opening list file");
  //       return 1;
  //     }
  //   }
  //   if (verbose_flag) printf("Listing k-path...\n");
  //   list_k_path(list_fd, vector<int>(), setBit(0, color[N]), N);
  //   if (verbose_flag) printf("%llu k-path found!\n", cont);
  //   fclose(list_fd);
  // }

  // cont = 0;
  // for (unsigned int i = 0; i <= k; i++)
  //   for (unsigned int j = 0; j <= N; j++) cont += DP[i][j].size();
  // if (verbose_flag) printf("DP elements: %llu\n", cont);
  //
  // cont = 0;
  // for (auto l : links) cont += l.second.size();
  // if (verbose_flag) printf("Oracle links: %llu\n", cont);

  return 0;
}
