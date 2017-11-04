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

#ifdef Q_8
#define MAXQ 8
#define COLORSET uint8_t
#elif Q_16
#define MAXQ 16
#define COLORSET uint16_t
#elif Q_64
#define MAXQ 64
#define COLORSET uint64_t
#else
#define MAXQ 32
#define COLORSET uint32_t
#endif

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif(defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) ||       \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
#include <stdio.h>

size_t getCurrentRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
  /* OSX ------------------------------------------------------ */
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                &infoCount) != KERN_SUCCESS)
    return (size_t)0L; /* Can't access? */
  return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE *fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL)
    return (size_t)0L; /* Can't open? */
  if (fscanf(fp, "%*s%ld", &rss) != 1) {
    fclose(fp);
    return (size_t)0L; /* Can't read? */
  }
  fclose(fp);
  return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
  return (size_t)0L;
#endif
}

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

using namespace std;
typedef long long ll;

unsigned int N, E;
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

// Random coloring graph using q color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = eng() % q;
}

// Path label
string L(vector<int> P) {
  string l = "";
  for (size_t i = 0; i < P.size(); i++) l += label[P[i]];
  return l;
}

// bruteforce
set<string> dict;
map<pair<int, string>, ll> freqBrute;

vector<int> P[30];
string Pstring[30];
set<int> Pset[30];

void dfs(int t, int u, int k) {
  if (Pset[t].find(u) != Pset[t].end()) return;

  Pset[t].insert(u);
  Pstring[t].push_back(label[u]);
  P[t].push_back(u);

  if (k == 0) {
    #pragma omp critical
    {
      dict.insert(Pstring[t]);
      freqBrute[make_pair(*P[t].begin(), Pstring[t])]++;
    }
  } else {
    for (int v : G[u]) dfs(t, v, k - 1);
  }
  Pset[t].erase(u);
  Pstring[t].pop_back();
  P[t].pop_back();
}

// Dynamic Programming
map<COLORSET, ll> *M[MAXQ + 1];

void processDP() {
  if (verbose_flag) printf("Q = %u\n", 1);

  #pragma omp parallel for schedule(static, 1)
  for (unsigned int u = 0; u < N; u++)
    M[1][u][setBit(0, color[u])] = 1ll;

  for (unsigned int i = 2; i <= q; i++) {
    #pragma omp parallel for schedule(static, 1)
    for (unsigned int u = 0; u < N; u++) {
      for (int v : G[u]) {
        for (auto d : M[i - 1][v]) {
          COLORSET s = d.first;
          ll f = d.second;
          if (getBit(s, color[u])) continue;
          ll fp = M[i][u][setBit(s, color[u])];
          M[i][u][setBit(s, color[u])] = f + fp;
        }
      }
    }
  }
}

bool isPrefix(set<string> W, string x) {
  auto it = W.lower_bound(x);
  if (it == W.end()) return false;
  return mismatch(x.begin(), x.end(), (*it).begin()).first == x.end();
}

map<string, ll> processFrequency(set<string> W, multiset<int> X) {
  set<string> WR;
  for (string w : W) {
    reverse(w.begin(), w.end());
    WR.insert(w);
  }
  vector<tuple<int, string, COLORSET>> old;

  for (int x : X)
    if (isPrefix(WR, string(&label[x], 1)))
      old.push_back(make_tuple(x, string(&label[x], 1), setBit(0ll, color[x])));

  for (int i = q - 1; i > 0; i--) {
    vector<tuple<int, string, COLORSET>> current;
    #pragma omp parallel for schedule(static, 1)
    for (int j = 0; j < (int)old.size(); j++) {
      auto o = old[j];
      int u = get<0>(o);
      string LP = get<1>(o);
      COLORSET CP = get<2>(o);
      for (int v : G[u]) {
        if (getBit(CP, color[v])) continue;
        COLORSET CPv = setBit(CP, color[v]);
        string LPv = LP + label[v];
        if (!isPrefix(WR, LPv)) continue;
        #pragma omp critical
        { current.push_back(make_tuple(v, LPv, CPv)); }
      }
    }
    old = current;
  }
  map<string, ll> frequency;
  for (auto c : old) {
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
    for (int v : G[u]) freq.push_back(M[i][v][D]);
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
  for (int x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while (R.size() < (size_t)r) {
    int u = X[distribution(eng)];
    vector<int> P = randomPathTo(u);
    if (R.find(P) == R.end()) R.insert(P);
  }
  for (auto r : R) {
    reverse(r.begin(), r.end());
    W.insert(L(r));
  }
  return W;
}

map<pair<int, string>, ll> randomColorfulSamplePlus(vector<int> X, int r) {
  map<pair<int, string>, ll> W;
  set<vector<int>> R;
  vector<ll> freqX;
  for (int x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while (R.size() < (size_t)r) {
    int u = X[distribution(eng)];
    vector<int> P = randomPathTo(u);
    if (R.find(P) == R.end()) R.insert(P);
  }
  for (auto r : R) {
    reverse(r.begin(), r.end());
    W[make_pair(*r.begin(), L(r))]++;
  }
  return W;
}

set<string> BCSampler(set<int> A, set<int> B, int r) {
  vector<int> X;
  for (int a : A) X.push_back(a);
  for (int b : B) X.push_back(b);
  return randomColorfulSample(X, r);
}

vector<int> naiveRandomPathTo(int u) {
  vector<int> P;
  set<int> Ps;
  P.push_back(u);
  Ps.insert(u);
  for (int i = q - 1; i > 0; i--) {
    vector<int> Nu;
    for (int j : G[u])
      if (Ps.find(j) == Ps.end()) Nu.push_back(j);
    if (Nu.size() == 0) return P;
    int u = Nu[rand() % Nu.size()];
    Ps.insert(u);
    P.push_back(u);
  }
  return P;
}

map<pair<int, string>, ll> baselineSampler(vector<int> X, int r) {
  set<vector<int>> R;
  while (R.size() < (size_t)r) {
    int u = X[rand() % X.size()];
    vector<int> P = naiveRandomPathTo(u);
    if (P.size() == q && R.find(P) == R.end()) R.insert(P);
  }
  map<pair<int, string>, ll> fx;
  for (auto P : R) fx[make_pair(*P.begin(), L(P))]++;
  return fx;
}

double FJW(set<string> W, set<int> A, set<int> B) {
  multiset<int> AiB, AB;
  for (int a : A) AB.insert(a);
  for (int b : B)
    if (AB.find(b) == AB.end()) AB.insert(b);
  for (int a : A)
    if (B.find(a) != B.end()) AiB.insert(a);
  long long num = 0ll;
  long long den = 0ll;
  map<string, ll> freqAiB = processFrequency(W, AiB);
  map<string, ll> freqAB = processFrequency(W, AB);
  for (string w : W) {
    num += freqAiB[w];
    den += freqAB[w];
  }
  return (double)num / (double)den;
}

double BCW(set<string> W, map<string, ll> freqA, map<string, ll> freqB,
           long long R) {
  ll num = 0ll;
  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    num += 2 * min(fax, fbx);
  }
  return (double)num / (double)R;
}

double BCW(set<string> W, map<string, ll> freqA, map<string, ll> freqB) {
  ll num = 0ll;
  ll den = 0ll;
  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    num += 2 * min(fax, fbx);
    den += fax + fbx;
  }
  return (double)num / (double)den;
}

double FJW(set<string> W, map<string, ll> freqA, map<string, ll> freqB, long long R) {
  ll num = 0ll;
  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    num += min(fax, fbx);
  }
  return (double)num / (double)R;
}

double BCW(set<string> W, set<int> A, set<int> B) {
  ll num = 0ll;
  ll den = 0ll;
  multiset<int> mA, mB;
  for (int a : A) mA.insert(a);
  for (int b : B) mB.insert(b);
  map<string, ll> freqA = processFrequency(W, mA);
  map<string, ll> freqB = processFrequency(W, mB);
  vector<string> vW = vector<string>(W.begin(), W.end());
  for (int i = 0; i < (int)vW.size(); i++) {
    string w = vW[i];
    long long fax = freqA[w];
    long long fbx = freqB[w];
    num += 2 * min(fax, fbx);
    den += fax + fbx;
  }
  return (double)num / (double)den;
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
  gettimeofday(&te, NULL);  // get current time
  long long milliseconds = te.tv_sec * 1000LL + te.tv_usec / 1000;
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

  if (q > MAXQ) {
    printf("q to high! (max value: %d\n", MAXQ);
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
    read(input_fd, &E, sizeof(int));
    if (verbose_flag) printf("N = %d | E = %d\n", N, E);

    label = new char[N + 1];
    color = new int[N + 1];
    int *intLabel = new int[N + 1];

    if (verbose_flag) printf("Reading labels...\n");
    read(input_fd, intLabel, N * sizeof(int));
    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + intLabel[i];

    if (verbose_flag) printf("Reading edges...\n");
    G = new vector<int>[N + 1];
    int *ab = new int[2 * E];
    read(input_fd, ab, 2 * E * sizeof(int));
    for (unsigned int i = 0; i < E; i++) {
      G[ab[2 * i]].push_back(ab[2 * i + 1]);
      G[ab[2 * i + 1]].push_back(ab[2 * i]);
    }
    free(ab);

  } else {
    // Read from stdin, nme format
    N = nextInt();
    E = nextInt();
    if (verbose_flag) printf("N = %d | E = %d\n", N, E);

    label = new char[N + 1];
    color = new int[N + 1];
    G = new vector<int>[N + 1];

    if (verbose_flag) printf("Reading labels...\n");
    for (unsigned int i = 0; i < N; i++) label[i] = 'A' + nextInt();

    if (verbose_flag) printf("Reading edges...\n");
    for (unsigned int i = 0; i < E; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  if (verbose_flag) printf("N = %d | E = %d\n", N, E);

  if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

  // Create DP Table
  for (unsigned int i = 0; i <= q + 1; i++)
    M[i] = new map<COLORSET, ll>[N + 1];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Add fake node N connected to all the nodes
  //  for (unsigned int i = 0; i < N; i++) {
  //    G[N].push_back(i);
  //    G[i].push_back(N);
  //  }
  //  color[N] = q;
  //  N++;
  //  q++;

  // Fill dynamic programming table
  if (verbose_flag) printf("Processing DP table...\n");
  ll time_a = current_timestamp();
  processDP();
  ll time_b = current_timestamp() - time_a;
  if (verbose_flag) printf("End processing DP table [%llu]ms\n", time_b);

  ll time_dp = time_b;
  // N--;
  // q--;
  // for (unsigned int i = 0; i < N; i++) G[i].pop_back();

  vector<int> sampleV;
  for (unsigned int i = 0; i < N; i++) sampleV.push_back(i);

  vector<pair<int, int>> ABsize;
  // ABsize.push_back(make_pair(10,10));
  // ABsize.push_back(make_pair(10, 100));
  ABsize.push_back(make_pair(100, 100));

  vector<double> epsilonV;
  epsilonV.push_back(0.05);
  // epsilonV.push_back(0.4);
  // epsilonV.push_back(0.6);

  dict.clear();
  freqBrute.clear();

  double bc_brute;  // BC-BRUTE
  double bc_alg3;   // BC-2PLUS
  double bc_2plus;  // BC-2PLUS
  double bc_base;   // BC-2PLUS

  double bc_alg3_rel;   // BC-2PLUS_REL
  double bc_2plus_rel;  // BC-2PLUS_REL
  double bc_base_rel;   // BC-2PLUS_REL

  double fj_brute;  // FJ-BRUTE_REL
  double fj_alg3;   // FJ-2PLUS
  double fj_2plus;  // FJ-2PLUS
  double fj_base;   // FJ-2PLUS

  double fj_2plus_rel;  // FJ-2PLUS_REL
  double fj_alg3_rel;   // FJ-2PLUS_REL
  double fj_base_rel;   // FJ-2PLUS_REL

  int tau_brute;  // TAU
  int tau_alg3;   // TAU
  int tau_base;   // TAU
  int tau_2plus;  // TAU

  ll time_brute = 0ll;  // TIME
  ll time_alg3 = 0ll;   // TIME
  ll time_2plus = 0ll;  // TIME
  ll time_base = 0ll;   // TIME

  srand(42);

  for(auto ABs : ABsize)
  {
    random_shuffle(sampleV.begin(), sampleV.end());
    set<int> A = set<int>(sampleV.begin(), sampleV.begin() + ABs.first);
    random_shuffle(sampleV.begin(), sampleV.end());
    set<int> B = set<int>(sampleV.begin(), sampleV.begin() + ABs.second);
    vector<int> X;
    for (int a : A) X.push_back(a);
    for (int b : B) X.push_back(b);
    set<int> AB;
    for (int a : A) AB.insert(a);
    for (int b : B) AB.insert(b);
    vector<int> ABv = vector<int>(AB.begin(), AB.end());

    printf("|A| = %4d |B| = %4d\n",ABs.first, ABs.second);

    for(double epsilon : epsilonV)
    {
      printf("epsilon=%0.2f\n", epsilon);

      int R = 0;
      long long PAB = 0;
      for(int x : AB)
        for(auto y : M[q][x])
          PAB += y.second;

      R = log((double)PAB)/(epsilon*epsilon);

      map<string, ll> freqA, freqB;
      set<string> W;
      double bcw, fjw;
      long long Rp = 0ll;

      // BRUTE-FORCE
      printf("BRUTEFORCE\n");
      dict.clear();
      freqA.clear();
      freqB.clear();
      freqBrute.clear();
      time_brute = current_timestamp();
      #pragma omp parallel for schedule(static, 1)
      for (size_t i = 0; i < ABv.size(); i++) {
        int tid = omp_get_thread_num();
        dfs(tid, ABv[i], q - 1);
      }

      double realBC, realFJ;
      for (auto w : freqBrute) {
        int u = w.first.first;
        string s = w.first.second;
        ll freq = w.second;
        if (A.find(u) != A.end()) {
          Rp += freq;
          freqA[s] += freq;
        }
        if (B.find(u) != B.end()) {
          Rp += freq;
          freqB[s] += freq;
        }
      }

      time_brute = current_timestamp() - time_brute;
      tau_brute = dict.size();
      bc_brute = realBC = bcw = BCW(dict, freqA, freqB);
      fj_brute = realFJ = fjw = FJW(dict, freqA, freqB, Rp);

      for(int exp = 0 ; exp < 100 ; exp++)
      {

        // BASELINE
        W.clear();
        freqA.clear();
        freqB.clear();

        time_base = current_timestamp();
        map<pair<int, string>, ll> BLsampling = baselineSampler(X, R);
        for (auto w : BLsampling) {
          int u = w.first.first;
          W.insert(w.first.second);
          if (A.find(u) != A.end()) freqA[w.first.second] += w.second;
          if (B.find(u) != B.end()) freqB[w.first.second] += w.second;
        }

        tau_base = W.size();
        time_base = current_timestamp() - time_base;
        bc_base = BCW(W, freqA, freqB);
        bc_base_rel = abs(bcw - realBC) / realBC;

        W.clear();
        freqA.clear();
        freqB.clear();
        Rp = 0;

        BLsampling = baselineSampler(ABv, R);
        for (auto w : BLsampling) {
          int u = w.first.first;
          W.insert(w.first.second);
          if (A.find(u) != A.end()) freqA[w.first.second] += w.second;
          if (B.find(u) != B.end()) freqB[w.first.second] += w.second;
        }

        for (auto a : freqA) Rp += a.second;
        for (auto b : freqB) Rp += b.second;

        fj_base = FJW(W, freqA, freqB, (long long)Rp);
        fj_base_rel = abs(fjw - realFJ) / realFJ;

        // COLORFUL SAMPLER OK
        W.clear();
        freqA.clear();
        freqB.clear();

        time_alg3 = current_timestamp();
        set<string> Sample = randomColorfulSample(X, R);
        freqA = processFrequency(Sample, multiset<int>(A.begin(), A.end()));
        freqB = processFrequency(Sample, multiset<int>(B.begin(), B.end()));

        time_alg3 = current_timestamp() - time_alg3;
        tau_alg3 = Sample.size();
        bc_alg3 = BCW(Sample, freqA, freqB);
        bc_alg3_rel = abs(bcw - realBC) / realBC;

        W.clear();
        freqA.clear();
        freqB.clear();

        Sample = randomColorfulSample(vector<int>(AB.begin(), AB.end()), R);
        freqA = processFrequency(Sample, multiset<int>(A.begin(), A.end()));
        freqB = processFrequency(Sample, multiset<int>(B.begin(), B.end()));
        Rp = 0;
        for (auto a : freqA) Rp += a.second;
        for (auto b : freqB) Rp += b.second;
        fj_alg3 = FJW(Sample, freqA, freqB, Rp);
        fj_alg3_rel = abs(fjw - realFJ) / realFJ;

        // COLORFULSAMPLER PLUS
        W.clear();
        freqA.clear();
        freqB.clear();
        Rp = 0;
        time_2plus = current_timestamp();

        map<pair<int, string>, ll> SamplePlus = randomColorfulSamplePlus(X, R);

        for (auto w : SamplePlus) {
          int u = w.first.first;
          W.insert(w.first.second);
          if (A.find(u) != A.end()) {
            freqA[w.first.second] += w.second;
            Rp += w.second;
          }
          if (B.find(u) != B.end()) {
            freqB[w.first.second] += w.second;
            Rp += w.second;
          }
        }

        time_2plus = current_timestamp() - time_2plus;
        tau_2plus = W.size();
        bc_2plus = BCW(W, freqA, freqB, Rp);
        bc_2plus_rel = abs(bcw - realBC) / realBC;

        W.clear();
        freqA.clear();
        freqB.clear();
        Rp = 0ll;

        SamplePlus = randomColorfulSamplePlus(vector<int>(AB.begin(), AB.end()), R);

        for (auto w : SamplePlus) {
          int u = w.first.first;
          W.insert(w.first.second);
          if (A.find(u) != A.end()) {
            freqA[w.first.second] += w.second;
            Rp += w.second;
          }
          if (B.find(u) != B.end()) {
            freqB[w.first.second] += w.second;
            Rp += w.second;
          }
        }
        fj_2plus = fjw = FJW(W, freqA, freqB, (long long)Rp);
        fj_2plus_rel = abs(fjw - realFJ) / realFJ;

        printf("%.1f,", epsilon);               // Q
        printf( "%2d,", q);                     // Q
        printf( "%4d,", R);                     // R
        printf("%4zu,", A.size());              // HA
        printf("%4zu,", B.size());              // HB

        printf( "%.6f,", bc_brute);             // BC-BRUTE
        printf( "%.6f,", fj_brute);             // FJ-BRUTE
        printf(  "%4d,", tau_brute);            // TAU
        printf("%4llu,", time_brute);           // TIME

        printf( "%.6f,", bc_2plus);             // BC-2PLUS
        printf( "%.6f,", bc_2plus_rel);         // BC-2PLUS
        printf( "%.6f,", fj_2plus);             // FJ-2PLUS
        printf( "%.6f,", fj_2plus_rel);         // FJ-2PLUS
        printf(  "%4d,", tau_2plus);            // TAU
        printf("%4llu,", time_dp + time_2plus); // TIME

        printf( "%.6f,", bc_base);              // BC-BASE
        printf( "%.6f,", bc_base_rel);          // BC-BASE
        printf( "%.6f,", fj_base);              // FJ-BASE
        printf( "%.6f,", fj_base_rel);          // FJ-BASE
        printf(  "%4d,", tau_base);             // TAU
        printf("%4llu,", time_base);            // TIME

        printf( "%.6f,", bc_alg3);              // BC-ALG3
        printf( "%.6f,", bc_alg3_rel);          // BC-ALG3
        printf( "%.6f,", fj_alg3);              // FJ-ALG3
        printf( "%.6f,", fj_alg3_rel);          // FJ-ALG3
        printf(  "%4d,", tau_alg3);             // TAU
        printf( "%4llu", time_dp + time_alg3);  // TIME

        printf("\n");

      }
    }
  }
  return 0;
}
