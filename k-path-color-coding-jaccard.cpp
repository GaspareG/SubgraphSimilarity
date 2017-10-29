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
int Sa, Sb;
int *A, *B;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Random generator
random_device rd;
mt19937_64 eng = mt19937_64(rd());
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
  for (unsigned int i = 0; i < N; i++) color[i] = eng() % kp;
}

// Path label
string L(vector<int> P)
{
    string l = "";
    for(size_t i=0; i<P.size(); i++) l += label[P[i]];
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
  // Base case
  for (unsigned int u = 0; u < N; u++) DP[1][u][setBit(0, color[u])] = 1ll;
  // Induction
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

void backProp() {
  for (int i = k - 2; i >= 0; i--) {
    if (verbose_flag) printf("K = %d\n", i);
    #pragma omp parallel for
    for (unsigned int x = 0; x <= N; x++) {
      vector<ll> toDel;
      for (auto cf : DP[i][x]) {
        COLORSET C = cf.first;
        bool find = false;

        for (int j : G[x]) {
          if (getBit(C, color[j])) continue;

          if (DP[i + 1][j].find(setBit(C, color[j])) != DP[i + 1][j].end()) {
            find = true;
            #pragma omp critical
            { addLink(x, C, j); }
          }
        }
        if (!find) toDel.push_back(C);
      }
      for (COLORSET tod : toDel){
        // printf("\t\tCANCELLO %llu da [%d][%d]\n",tod,i,x);
         DP[i][x].erase(tod);
       }
    }
  }
}

// Random sampling
vector<int> randomPathFrom(int u)
{
  vector<int> P;
  P.push_back(u);
  COLORSET D = color[u];
  for(int i=k-1; i>0; i--)
  {
    vector<ll> freq;
    for(int v : G[u]) freq.push_back( DP[i][v][getCompl(D)] );
    discrete_distribution<int> distribution(freq.begin(), freq.end());
    u = G[u][ distribution(eng) ];
    P.push_back(u);
    D = setBit(D, color[u]);
  }
  return P;
}

// Random sampling
set<string> randomColorfulSample(vector<int> X, int r)
{
  set<string> W;
  set< vector<int> > R;
  vector<ll> freqX;

  ll last = 0ll;
  ll sum = 0ll;
  for(int x : X)
  {
    ll freq = DP[k-1][x][getCompl(0)];
    freqX.push_back(last);
    last = freq;
    sum += last;
  }

  while( r )
  {
    ll rndIdx = distr(eng) % sum;
    int u = X[distance(freqX.begin(), upper_bound(freqX.begin(), freqX.end(), rndIdx) - 1)];

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
  set<vector<int>> R;

  vector<int> AB;
  vector<ll> freqAB;
  for(int a : A) AB.push_back(a);
  for(int b : B) AB.push_back(b);
  sort(AB.begin(), AB.end());
  AB.erase( unique(AB.begin(), AB.end()), AB.end());
  for(int v : AB)
  {
    ll freq = DP[k][v][getCompl(0ll)];
    ll alpha = 0ll;
    if( A.find(v) != A.end() ) alpha++;
    if( B.find(v) != B.end() ) alpha++;
    freqAB.push_back(alpha*freq);
    // printf("FREQ[%6d] = [%6lld][%6lld][%6lld]\n",v,freq, alpha, alpha*freq);
  }

  discrete_distribution<int> distribution(freqAB.begin(), freqAB.end());

  while(r)
  {
    int u = AB[distribution(eng)];
    printf("\t R = %d U = %d\n", r, u);
    vector<int> P = randomPathFrom(u);
    if( R.find(P) == R.end() )
    {
      R.insert(P);
      r--;
    }
  }

  for(auto r : R) W.insert(L(r));
  return W;
}

// set<string> BCSampler(set<int> A, set<int> B, int r)
// {
//   set<string> W;
//   set< vector<int> > R;
//   vector<int> AB;
//   vector<ll> freqAB;
//   for(int a : A) AB.push_back(a);
//   for(int b : B) AB.push_back(b);
//   sort(AB.begin(), AB.end() );
//   AB.erase( unique( AB.begin(), AB.end() ), AB.end() );
//   ll last = 0ll;
//   ll sum = 0ll;
//   for(int x : AB)
//   {
// //    printf("for[%d]\n", x);
//     ll freq = DP[k-1][x].begin()->second;
// //    for(auto y : DP[k-1][x] )
// //      printf("[%d][%d] = [%u][%lld]\n",k,x,y.first, y.second);
//     ll alpha = 0;
//     if( A.find(x) != A.end() ) alpha++;
//     if( B.find(x) != B.end() ) alpha++;
//     freqAB.push_back(sum);
//     // printf("[%6d] [%6lld] [%6lld] [%6lld] [%6lld]\n", x, freq, alpha, last, sum);
//     last = freq*alpha;
//     sum += last;
//   }
//   // printf("\tBCSampler (%lld) [%zu][%zu]\n", sum, AB.size(), freqAB.size());
//
//   for(size_t i=0; i<AB.size(); i++)
//   {
//     // printf("TB[%d][%lld]\n", AB[i], freqAB[i]);
//   }
//
//   while( r )
//   {
//     // printf("OK %d\n", r);
//     ll rndIdx = distr(eng) % sum;
//     // printf("\trndIdx %lld\n", rndIdx);
//     int u = AB[distance(freqAB.begin(), upper_bound(freqAB.begin(), freqAB.end(), rndIdx) - 1)];
//     // printf("\tu= %d\n", u);
//     vector<int> P = randomPathFrom(u);
//     // printf("FIND: [%zu] [%s]\n", P.size(), L(P).c_str() );
//     // TODO ??
//     if( R.find(P) == R.end() )
//     {
//       R.insert(P);
//       r--;
//     }
//   }
//   // printf("|R| = %zu\n", R.size());
//   for(vector<int> r : R)  W.insert( L(r) );
//   // printf("|W| = %zu\n", W.size());
//   return W;
// }

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

    label = new char[N+1];
    color = new int[N + 1];
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
    color = new int[N + 1];
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

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

  // for (unsigned int i = 0; i < N; i++) {
  //   G[N].push_back(i);
  //   G[i].push_back(N);
  // }

  // Create DP Table
  for (unsigned int i = 0; i <= k + 1; i++)
    DP[i] = new map<COLORSET, ll>[N + 1];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();
  // color[N] = kp;
  // label[N] = 'Z';
  // k++;

  // Fill dynamic programming table
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("End processing DP table...\n");


  // Fill dynamic programming table
  // if (verbose_flag) printf("Backward-propagation...\n");
  // backProp();
  // if (verbose_flag) printf("End Backward-propagation...\n");

  set<int> vA = set<int>(A, A+Sa);
  set<int> vB = set<int>(B, B+Sb);

  if (verbose_flag) printf("Sampling strings...\n");
  set<string> W = BCSampler(vA, vB, 10);

  if (verbose_flag) printf("Sampled strings:\n");
  for(string w : W) printf("%s\n", w.c_str());

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
