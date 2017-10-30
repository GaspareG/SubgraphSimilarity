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
COLORSET getCompl(COLORSET n) { return ((1 << k) - 1) & (~n); }

// Random coloring graph using kp color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = eng() % k;
}

// Path label
string L(vector<int> P) {
  string l = "";
  for (size_t i = 0; i < P.size(); i++) l += label[P[i]];
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

void list_k_path(vector<int> ps, COLORSET cs, int x) {
  vector<int> oracle = H(x, cs);
  if (ps.size() + 1 == k) {
    cont++;
    for (int j : ps) printf("[%6d] ", j);
    printf("\n");
    for (int j : ps) printf("[%6d] ", color[j]);
    printf("\n");
  } else
    for (int v : oracle) {
      // printf("LINK TO %d\n",v);
      ps.push_back(v);
      list_k_path(ps, setBit(cs, color[v]), v);
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

bool isPrefix(set<string> W, string x)
{
  auto it = W.lower_bound(x);
  if( it == W.end() ) return false;
  return mismatch(x.begin(), x.end(), (*it).begin()).first == x.end() ;
}

map<string, ll> processFrequency(set<string> W, multiset<int> X)
{
  set<string> WR;
  for(string w : W)
  {
    reverse(w.begin(), w.end());
    WR.insert(w);
  }
  multiset< tuple<int, string, COLORSET> > old;

  for(int x : X)
    if( isPrefix(WR, string(&label[x],1)) )
      old.insert(make_tuple(x, string(&label[x],1), setBit(0ll, color[x])));

  for(int i=k-1; i>0; i--)
  {
    multiset< tuple<int, string, COLORSET> > current;
    for(auto o : old)
    {
      int u = get<0>(o);
      string LP = get<1>(o);
      COLORSET CP = get<2>(o);
      for(int v : G[u])
      {
        if(getBit(CP, color[v])) continue;
        COLORSET CPv = setBit(CP, color[v]);
        string LPv = LP + label[v];
        if( !isPrefix(WR, LPv) ) continue;
        current.insert(make_tuple(v, LPv, CPv));
      }
    }
    old = current;
  }

  map<string, ll> frequency;
  for(auto c : old)
  {
    string s = get<1>(c);
    reverse(s.begin(), s.end());
    frequency[s]++;
  }
  return frequency;
}

void backProp() {
  // for (int i = k - 1; i > 0; i--) {
  //   if (verbose_flag) printf("K = %d\n", i);
  //   #pragma omp parallel for
  //   for (unsigned int x = 0; x < N; x++) {
  //     vector<COLORSET> toDelete;
  //     for (pair<COLORSET, ll> CF : DP[i][x]) {
  //       bool find = false;
  //       COLORSET C = CF.first;
  //       for (int j : G[x]) {
  //         if (color[j] == color[x]) continue;
  //         if (DP[i + 1][j].find(setBit(C, color[j])) != DP[i + 1][j].end()) {
  //           find = true;
  //           break;  // Don't stop if we need to construct the oracle H
  //           // addLink(x, C, j);
  //         }
  //       }
  //       if (!find) toDelete.push_back(C);
  //     }
  //     for (COLORSET C : toDelete) DP[i][x].erase(C);
  //   }
  // }
}

vector<int> randomPathTo(int u) {
  vector<int> P;
  P.push_back(u);
  COLORSET D = getCompl(setBit(0l, color[u]));
  for (int i = k - 1; i > 0; i--) {
    vector<ll> freq;
    for (int v : G[u]) freq.push_back(DP[i][v][D]);
    discrete_distribution<int> distribution(freq.begin(), freq.end());
    u = G[u][distribution(eng)];
    P.push_back(u);
    D = clearBit(D, color[u]);
  }
  reverse(P.begin(), P.end());
  return P;
}

set<string> randomColorfulSample(vector<int> X, int r) {
  set<string> W;
  set<vector<int>> R;
  vector<ll> freqX;
  for (int x : X) freqX.push_back(DP[k][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while (R.size() < (size_t)r) {
    int u = X[distribution(eng)];
    vector<int> P = randomPathTo(u);
    if (R.find(P) == R.end()) {
      // for (int p : P) printf("[%6d] ", p);
      // printf("\n");
      // for (int p : P) printf("[%6d] ", color[p]);
      // printf("\n");
      R.insert(P);
    }
  }
  for (auto r : R) W.insert(L(r));
  return W;
}

set<string> BCSampler(set<int> A, set<int> B, int r) {
  vector<int> X;
  for (int a : A) X.push_back(a);
  for (int b : B) X.push_back(b);
  return randomColorfulSample(X, r);
}

double FJW(set<string> W, set<int> A, set<int> B)
{
  multiset<int> AiB, AB;
  for(int a : A) AB.insert(a);
  for(int b : B) if( AB.find(b) == AB.end() )  AB.insert(b);
  for(int a : A) if(  B.find(a) !=  B.end() ) AiB.insert(a);
  long long num = 0ll;
  long long den = 0ll;
  map<string, ll> freqAiB = processFrequency(W, AiB);
  map<string, ll> freqAB = processFrequency(W, AB);
  for(string w : W)
  {
    num += freqAiB[w];
    den += freqAB[w];
  }
  return (double) num / (double) den;
}

double BCW(set<string> W, set<int> A, set<int> B)
{
  long long num = 0ll;
  long long den = 0ll;
  multiset<int> mA, mB;
  for(int a : A) mA.insert(a);
  for(int b : B) mB.insert(b);
  map<string, ll> freqA = processFrequency(W, mA);
  map<string, ll> freqB = processFrequency(W, mB);
  for(string w : W)
  {
    long long fax = freqA[w];
    long long fbx = freqB[w];
    num += 2 * min(fax, fbx);
    den += fax + fbx;
  }
  return (double) num / (double) den;
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
    color = new int[N + 1];
    int *intLabel = new int[N + 1];

    read(input_fd, intLabel, N * sizeof(int));
    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + intLabel[i];

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

    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + nextInt();

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }

    Sa = nextInt();
    Sb = nextInt();

    for (int i = 0; i < Sa; i++) A[i] = nextInt();
    for (int i = 0; i < Sb; i++) B[i] = nextInt();
  }

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

  // Create DP Table
  for (unsigned int i = 0; i <= k + 1; i++) DP[i] = new map<COLORSET, ll>[N + 1];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Add fake node N connected to all the nodes
  for (unsigned int i = 0; i < N; i++) {
    G[N].push_back(i);
    G[i].push_back(N);
  }
  color[N] = k;
  N++;
  k++;

  // Fill dynamic programming table
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("End processing DP table...\n");

  // Backward-propagation of DP
  if (verbose_flag) printf("Backward-propagation...\n");
  backProp();
  if (verbose_flag) printf("End Backward-propagation...\n");

  // list_k_path(vector<int>(), setBit(0ll, color[N-1]), N-1);
  N--;
  k--;
  for (unsigned int i = 0; i < N; i++) G[i].pop_back();

  set<int> vA = set<int>(A, A + Sa);
  set<int> vB = set<int>(B, B + Sb);

  double sum = 0.;

  for(int i=0; i<1000; i++)
  {
    if (verbose_flag) printf("Sampling strings...\n");
    set<string> W = BCSampler(vA, vB, 1);

    if (verbose_flag) printf("Sampled strings:\n");
    for (string w : W) printf("%s\n", w.c_str());

    if (verbose_flag) printf("Find frequency(A+B)\n");
    multiset<int> mAB = multiset<int>(A, A + Sa);
    mAB.insert(B, B + Sb);
    map<string, ll> freqAB = processFrequency(W, mAB);
    // if (verbose_flag) printf("Freq(A+B):\n");
    // for(auto f : freqAB)
    //   printf("[%10s] = [%6lld]\n", f.first.c_str(), f.second);
    double bcw = BCW(W,vA,vB);
    sum += bcw;
    if (verbose_flag) printf("BCW(W,A,B) = %.6f\n", bcw);
  }
  printf("E[BCW(W,A,B)] = %.6f\n", sum/1000.);
  return 0;
}
