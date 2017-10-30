/*
  Author: Gaspare Ferraro
  Compute JaccardIndex on multiset for 2 subgraph A / B
*/
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
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
char *label;
vector<int> *G;
int Sa, Sb;
int *A, *B;
set<int> setA, setB;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Jaccard index computation
map<string, vector< vector<int> > > freqP; // path -> string
vector<bool> in;
vector<int> path;
string s;

void dfs(int u)
{
  if( in[u] ) return;

  path.push_back(u);
  s.push_back( label[u] );
  in[u] = true;

  if( s.size() == k ) freqP[s].push_back(path);
  else for(int v : G[u]) dfs(v);

  path.pop_back();
  s.pop_back();
  in[u] = false;
}

// J(A,B) = sum_{s}{ min(A_{s}, B_{s}) } / sum_{s}{ max(A_{s}, B_{s}) }
double jaccardIndex()
{
  return 0.;
//   ll numerator = 0;
//   ll denominator = dict.size();
//   for(string s : dict)
//   {
//     ll fa = freq[0][s];
//     ll fb = freq[1][s];
//     if( fa > 0 && fb > 0 )
//       numerator++;
//   }
//   printf("%lld %lld\n", numerator, denominator);
//   return (double) numerator / (double) denominator;
}

ll freqV(int a, string x)
{
  ll out = 0ll;
  for(auto p: freqP[x])
  {
    if( p[0] != a ) continue;
    out++;
  }
  return out;
}

ll freqV(set<int> A, string x)
{
  ll f = 0ll;
  for(int a : A) f += freqV(a, x);
  return f;
}

// BC(A,B) = 2*sum_{s}{ min(A_{s}, B_{s}) } / sum_{s}{ A_{s}+B_{s}}
double BrayCurtisIndex()
{
  ll numerator = 0;
  ll denominator = 0;
  for(auto sv : freqP)
  {
    string s = sv.first;
    reverse(s.begin(), s.end());
    ll fa = freqV(setA, s);
    ll fb = freqV(setB, s);
    numerator += min(fa, fb);
    denominator += fa+fb;
  }
  numerator *= 2;
  printf("%lld %lld\n", numerator, denominator);
  return (double) numerator / (double) denominator;
}

void print_usage(char *filename) {
  printf(
      "Usage: ./%s -k length -K number -g filename -p threadcount -s filename"
      "--help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of labeled graph in nmle format (default stdin)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\nPrint status messages.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'k'},
      {"color", required_argument, 0, 'K'},
      {"input", required_argument, 0, 'g'},
      {"parallel", required_argument, 0, 'p'},
      {"help", no_argument, &help_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  while (1) {
    // c = getopt_long (argc, argv, "k:K:g:f:t:T:l:p:", long_options,
    // &option_index);
    c = getopt_long(argc, argv, "k:K:g:p:", long_options, &option_index);

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

    label = new char[N + 1];
    int *intLabel = new int[N+1];

    read(input_fd, intLabel, N*sizeof(int));
    for(unsigned int i=0; i<N; i++)
      label[i] = 'A' + intLabel[i];

    G = new vector<int>[N + 1];
    int ab[2];
    for (unsigned int i = 0; i < M; i++) {
      read(input_fd, ab, 2 * sizeof(int));
      G[ab[0]].push_back(ab[1]);
      G[ab[1]].push_back(ab[0]);
    }

    read(input_fd, ab, 2 * sizeof(int));
    Sa = ab[0];
    Sb = ab[1];
    A = new int[Sa];
    B = new int[Sb];

    read(input_fd, A, Sa * sizeof(int));
    read(input_fd, B, Sb * sizeof(int));

  } else {

    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();

    label = new char[N + 1];
    G = new vector<int>[N + 1];

    for(unsigned int i=0; i<N; i++)
      label[i] = 'A' + nextInt();

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }

    Sa = nextInt();
    Sb = nextInt();

    for(int i=0; i<Sa; i++) A[i] = nextInt();
    for(int i=0; i<Sb; i++) B[i] = nextInt();

  }

  setA = set<int>(A, A+Sa);
  setB = set<int>(B, B+Sb);

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);
  if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

  for(unsigned int i=0; i<N; i++) in.push_back(false);
  s.reserve(k);

  if (verbose_flag) printf("Finding paths...\n");
  for(unsigned int i=0; i<N; i++) dfs(i);

  if (verbose_flag) printf("Count: %lld\n", freqP.size());

  if (verbose_flag) printf("Calculating Jaccard index...\n");
  double jaccard = jaccardIndex();
  printf("Jaccard index: %f\n", jaccard);

  if (verbose_flag) printf("Calculating Bray-Curtis index...\n");
  double brayCurtis = BrayCurtisIndex();
  printf("Bray-Curtis index: %f\n", brayCurtis);
  return 0;

}
