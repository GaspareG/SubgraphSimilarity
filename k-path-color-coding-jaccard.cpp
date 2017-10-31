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
#include <time.h>
#include <sys/time.h>

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
unsigned q = 0;
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
COLORSET getCompl(COLORSET n) { return ((1 << q) - 1) & (~n); }

// Random coloring graph using k color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = eng() % q;
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
  if (ps.size() + 1 == q) {
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
  for (unsigned int i = 2; i <= q; i++)
  {
    if (verbose_flag) printf("K = %u\n", i);
    #pragma omp parallel for schedule(static, 1)
    for (unsigned int u = 0; u < N; u++)
    {
      for (int v : G[u])
      {
        for (auto d : DP[i - 1][v])
        {
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
  vector< tuple<int, string, COLORSET> > old;

  for(int x : X)
    if( isPrefix(WR, string(&label[x],1)) )
      old.push_back(make_tuple(x, string(&label[x],1), setBit(0ll, color[x])));

  for(int i=q-1; i>0; i--)
  {
    vector< tuple<int, string, COLORSET> > current;
    #pragma omp parallel for schedule(dynamic)
    for(int j=0; j<(int)old.size(); j++)
    {
      auto o = old[j];
      int u = get<0>(o);
      string LP = get<1>(o);
      COLORSET CP = get<2>(o);
      for(int v : G[u])
      {
        if(getBit(CP, color[v])) continue;
        COLORSET CPv = setBit(CP, color[v]);
        string LPv = LP+label[v];
        if( !isPrefix(WR, LPv) ) continue;
        #pragma omp critical
        {
          current.push_back(make_tuple(v, LPv, CPv));
        }
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

vector<int> randomPathTo(int u) {
  vector<int> P;
  P.push_back(u);
  COLORSET D = getCompl(setBit(0l, color[u]));
  for (int i = q - 1; i > 0; i--) {
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
  for (int x : X) freqX.push_back(DP[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while (R.size() < (size_t)r) {
    int u = X[distribution(eng)];
    vector<int> P = randomPathTo(u);
    if (R.find(P) == R.end()) R.insert(P);
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

// vector<int> naiveRandomPathTo(int u)
// {
//   vector<int> P;
//   // set<int> Ps;
//   // for(int i=q-1; i>0; i--)
//   // {
//   //   vector<int> Nu;
//   //   for(int j : G[]]);
//   // }
//   return P;
// }
//
// map<string, ll> baselineSampler(vector<int> X, int r)
// {
//   set<vector<int>> R;
//   while(R.size() < (size_t)r)
//   {
//     int u = X[eng()%X.size()];
//     vector<int> P = naiveRandomPathTo(u);
//     if( P.size() == q && R.find(P) == R.end() ) R.insert(P);
//     }
//   }
//   map<string, ll> fx;
//   for(auto P : R) fx[L(P)]++;
//   return fx;
// }

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
  ll num = 0ll;
  ll den = 0ll;
  multiset<int> mA, mB;
  for(int a : A) mA.insert(a);
  for(int b : B) mB.insert(b);
  map<string, ll> freqA = processFrequency(W, mA);
  map<string, ll> freqB = processFrequency(W, mB);
  vector<string> vW = vector<string>(W.begin(), W.end());
  #pragma omp parallel for schedule(static, 1) reduction(+:num), reduction(+: den)
  for(int i=0; i<(int)vW.size(); i++)
  {
    string w = vW[i];
    long long fax = freqA[w];
    long long fbx = freqB[w];
    num += 2 * min(fax, fbx);
    den += fax + fbx;
  }
  return (double) num / (double) den;
}

void print_usage(char *filename) {
  printf(
      "Usage: ./%s -q length -g filename -p threadcount -s filename"
      "--help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-q, --path length\n");
  printf("\tLength of the path.\n");

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

long long current_timestamp() {
    struct timeval te;
    gettimeofday(&te, NULL); // get current time
    long long milliseconds = te.tv_sec*1000LL + te.tv_usec/1000;
    return milliseconds;
}

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'q'},
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
    c = getopt_long(argc, argv, "q:g:p:", long_options, &option_index);

    if (c == -1) break;

    switch (c) {
      case 'q':
        if (optarg != NULL) q = atoi(optarg);
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

  if (q == 0) {
    printf("Invalid or missing path length value.\n");
    return 1;
  }

  if (q > MAXK) {
    printf("q to high! (max value: %d\n", MAXK);
    return 1;
  }

  if (thread_count > 0 && (int)thread_count < omp_get_max_threads()) {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
  }

  if (verbose_flag) {
    printf("Options:\n");
    printf("Q = %d\n", q);
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
    if (verbose_flag) printf("N = %d | M = %d\n", N, M);

    label = new char[N + 1];
    color = new int[N + 1];
    int *intLabel = new int[N + 1];

    if (verbose_flag) printf("Reading labels...\n");
    read(input_fd, intLabel, N * sizeof(int));
    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + intLabel[i];

    if (verbose_flag) printf("Reading edges...\n");
    G = new vector<int>[N + 1];
    int *ab = new int[2*M];
    read(input_fd, ab, 2 * M * sizeof(int));
    for (unsigned int i = 0; i < M; i++) {
      G[ab[2*i]].push_back(ab[2*i+1]);
      G[ab[2*i+1]].push_back(ab[2*i]);
    }
    free(ab);
    // read(input_fd, ab, 2 * sizeof(int));
    // Sa = ab[0];
    // Sb = ab[1];
    // A = new int[Sa];
    // B = new int[Sb];
    //
    // read(input_fd, A, Sa * sizeof(int));
    // read(input_fd, B, Sb * sizeof(int));

  } else {
    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();
    if (verbose_flag) printf("N = %d | M = %d\n", N, M);

    label = new char[N + 1];
    color = new int[N + 1];
    G = new vector<int>[N + 1];

    if (verbose_flag) printf("Reading labels...\n");
    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + nextInt();

    if (verbose_flag) printf("Reading edges...\n");
    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }

    // Sa = nextInt();
    // Sb = nextInt();
    //
    // for (int i = 0; i < Sa; i++) A[i] = nextInt();
    // for (int i = 0; i < Sb; i++) B[i] = nextInt();
  }

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

  // Create DP Table
  for (unsigned int i = 0; i <= q + 1; i++) DP[i] = new map<COLORSET, ll>[N + 1];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Add fake node N connected to all the nodes
  for (unsigned int i = 0; i < N; i++) {
    G[N].push_back(i);
    G[i].push_back(N);
  }
  color[N] = q;
  N++;
  q++;

  // Fill dynamic programming table
  if (verbose_flag) printf("Processing DP table...\n");
  ll time_a = current_timestamp();
  processDP();
  ll time_b = current_timestamp() - time_a;
  if (verbose_flag) printf("End processing DP table [%llu]ms\n", time_b);

  // list_k_path(vector<int>(), setBit(0ll, color[N-1]), N-1);
  N--;
  q--;
  for (unsigned int i = 0; i < N; i++) G[i].pop_back();

  vector<int> sampleV;
  for(unsigned int i=0; i<N; i++) sampleV.push_back(i);
  random_shuffle(sampleV.begin(), sampleV.end());

  vector<int> size;
  size.push_back(10);
  size.push_back(25);
  size.push_back(50);
  size.push_back(100);
  size.push_back(250);
  size.push_back(500);
  size.push_back(1000);
  size.push_back(2500);
  size.push_back(5000);
  size.push_back(10000);
  size.push_back(25000);
  size.push_back(50000);

  for(size_t i=0; i<size.size(); i++)
  {
    printf("|A| = |B| = %d\n",size[i]);

    set<int> A = set<int>(sampleV.begin(), sampleV.begin()+size[i]);
    set<int> B = set<int>(sampleV.end()-size[i]-1, sampleV.end());
    multiset<int> AB = multiset<int>(A.begin(), A.end());
    AB.insert(B.begin(), B.end());


    if (verbose_flag) printf("Sampling 1000 string...\n");
    time_a = current_timestamp();
    set<string> W = BCSampler(A, B, 1000);
    time_b = current_timestamp() - time_a;
    if (verbose_flag) printf("End sampling 1000 string [%llu]ms\n",time_b);

    if (verbose_flag) printf("Calculate BCW(A,B)...\n");
    time_a = current_timestamp();
    double bcw = BCW(W,A,B);
    time_b = current_timestamp() - time_a;
    if (verbose_flag) printf("BCW(A,B) = [%.6f]  [%llu]ms\n",bcw,time_b);

    printf("\n\n");
  }

  // set<int> vA = set<int>(A, A + Sa);
  // set<int> vB = set<int>(B, B + Sb);
  // multiset<int> mAB = multiset<int>(A, A + Sa);
  // mAB.insert(B, B + Sb);
  //
  // set<int> AB = set<int>(mAB.begin(), mAB.end());
  // if (verbose_flag) printf("Freq{AB}(W)...\n");
  // time_a = clock();
  // for(string w : W)
  // {
  //   set<string> ws;
  //   ws.insert(w);
  //   processFrequency(ws, mAB);
  // }
  // time_b = clock()-time_a;
  // if (verbose_flag) printf("End Freq{AB}(W) [%.6f]sec\n", (float)(time_b)/CLOCKS_PER_SEC);

  //
  // double sum = 0.;
  //
  // for(int i=0; i<1000; i++)
  // {
  //   if (verbose_flag) printf("Sampling strings...\n");
  //   set<string> W = BCSampler(vA, vB, 1);
  //
  //   if (verbose_flag) printf("Sampled strings:\n");
  //   for (string w : W) printf("%s\n", w.c_str());
  //
  //   if (verbose_flag) printf("Find frequency(A+B)\n");
  //   map<string, ll> freqAB = processFrequency(W, mAB);
  //   // if (verbose_flag) printf("Freq(A+B):\n");
  //   // for(auto f : freqAB)
  //   //   printf("[%10s] = [%6lld]\n", f.first.c_str(), f.second);
  //   double bcw = BCW(W,vA,vB);
  //   sum += bcw;
  //   if (verbose_flag) printf("BCW(W,A,B) = %.6f\n", bcw);
  // }
  // printf("E[BCW(W,A,B)] = %.6f\n", sum/1000.);
  return 0;
}
