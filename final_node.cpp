/*
Author: Gaspare Ferraro
Count and find simple k colorful-path in a graph
using the color-coding technique (parallel version)
*/
#include <vector>
#include <set>
#include <map>
#include <string>
#include <list>
#include <algorithm>
#include <iterator>
#include <random>
#include <queue>
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
static int verbose_flag, help_flag, bruteforce_flag, fcount_flag, fsample_flag, baseline_flag;

ll cont = 0;
int *color;
char *label;
vector<int> *G;
int *A, *B;

// parameter
unsigned int q = 0;
unsigned thread_count = 0;
unsigned int seed = 42;
unsigned int R = 10;
unsigned int Sa = 10, Sb = 10;
unsigned int mod = 0;
unsigned int experiment = 10;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Random generator
mt19937_64 eng;
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
  #pragma omp parallel for schedule(guided)
  for (unsigned int u = 0; u < N; u++) M[1][u][setBit(0, color[u])] = 1ll;

  for (unsigned int i = 2; i <= q; i++) {
    #pragma omp parallel for schedule(guided)
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
    current.clear();
    #pragma omp parallel for schedule(guided)
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
  list<int> P;
  P.push_front(u);
  COLORSET D = getCompl(setBit(0l, color[u]));
  for (int i = q - 1; i > 0; i--) {
    vector<ll> freq;
    for (int v : G[u]) freq.push_back(M[i][v][D]);
    discrete_distribution<int> distribution(freq.begin(), freq.end());
    #pragma omp critical
    {
      u = G[u][distribution(eng)];
    }
    P.push_front(u);
    D = clearBit(D, color[u]);
  }
  vector<int> ret;
  ret.clear();
  ret = vector<int>(begin(P), end(P));
  return ret;
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
  freqX.clear();
  for (int x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while( R.size() < (size_t)r)
  {
    int rem = r - R.size();
    // #pragma omp parallel for schedule(guided)
    for(int i=0; i<rem; i++)
    {
      int u;
//      #pragma omp critical
//      {
        u = X[distribution(eng)];
//      }
      vector<int> P = randomPathTo(u);
      // #pragma omp critical
      // {
        R.insert(P);
      // }
    }
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
    int u;
    #pragma omp critical
    {
      u = Nu[eng() % Nu.size()];
    }
    Ps.insert(u);
    P.push_back(u);
  }
  return P;
}

map<pair<int, string>, ll> baselineSampler(vector<int> X, int r) {
  set<vector<int>> R;

  while( R.size() < (size_t)r)
  {
    int rem = r - R.size();
    #pragma omp parallel for schedule(guided)
    for(int i=0; i<rem; i++)
    {
      int u;
      #pragma omp critical
      {
        u = X[eng() % X.size()];
      }
      vector<int> P = naiveRandomPathTo(u);
      #pragma omp critical
      {
        R.insert(P);
      }
    }
  }
  map<pair<int, string>, ll> fx;
  for (auto P : R) fx[make_pair(*P.begin(), L(P))]++;
  return fx;
}

double BCW(set<string> W, map<string, ll> freqA, map<string, ll> freqB) {
  double ret = 0.;

  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    if( fax + fbx == 0 ) continue;
    ret += (double) min(fax, fbx) / ( fax + fbx );
  }
  return ((double) 2 / W.size()) * ret;
}

double BCW_old(set<string> W, map<string, ll> freqA, map<string, ll> freqB) {
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
  return (double)num / (double) R;
}

  vector<int> sampleV;
  set<int> randomChoose(int s, int mod)
  {
    set<int> out ;
    // Random set
    if( mod == 0 )
    {
      shuffle(sampleV.begin(), sampleV.end(), eng);
      out = set<int>(sampleV.begin(), sampleV.begin()+s);
    }
    // BFS from random node
    else if( mod == 1 )
    {
      int from = eng() % N;
      queue<int> bfsQ;
      bfsQ.push(from);
      while( out.size() < (size_t) s && !bfsQ.empty() )
      {
        int att = bfsQ.front();
        bfsQ.pop();
        if( out.find(att) == out.end() )
        {
          out.insert(att);
          for(int v : G[att])
          {
            if( out.find(v) != out.end() ) continue;
            bfsQ.push(v);
          }
        }
      }
    }
    return out;
  }

  void print_usage(char *filename) {
    printf(
      "Usage: ./%s -q length -g filename -p threadcount -s filename"
      "--help --verbose\n",
      filename);
      printf("Valid arguments:\n");

      printf("-g, --input filename\n");
      printf("\tInput file of labeled graph in nmle format (default stdin)\n");

      printf("-p, --parallel threadcount\n");
      printf("\tNumber of threads to use (default maximum thread avaiable)\n");

      // Experiment option
      printf("-Q, --path length\n");
      printf("\tLength of the path.\n");

      printf("-S, --seed number\n");
      printf("\tSeed for random generator\n");

      printf("-R, --samplesize number\n");
      printf("\tSample size\n");

      printf("-A, --asize number\n");
      printf("\tSize of A\n");

      printf("-B, --bsize number\n");
      printf("\tSize of B\n");

      printf("-M, --modality number\n");
      printf("\tA & B choose (0=random, 1=component)\n");

      printf("-E, --experiment number\n");
      printf("\tNumber of experiments\n");

      printf("--bruteforce\n");
      printf("\tExecute bruteforce algorithm\n");

      printf("--baseline\n");
      printf("\tExecute baseline algorithm\n");

      printf("--fsample\n");
      printf("\tExecute f-sample algorithm\n");

      printf("--fcount\n");
      printf("\tExecute f-count algorithm\n");

      printf("--help\n");
      printf("\tDisplay help text and exit.\n");

      printf("--verbose\n");
      printf("\tPrint status messages.\n");
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

        // Execution option
        {   "input", required_argument, 0, 'g'},
        {"parallel", required_argument, 0, 'p'},

        // Test option
        {      "path", required_argument, 0, 'Q'},
        {      "seed", required_argument, 0, 'S'},
        {"samplesize", required_argument, 0, 'R'},
        {     "asize", required_argument, 0, 'A'},
        {     "bsize", required_argument, 0, 'B'},
        {  "modality", required_argument, 0, 'M'},
        {"experiment", required_argument, 0, 'E'},

        // Info flag
        {   "help", no_argument, &help_flag   , 1},
        {"verbose", no_argument, &verbose_flag, 1},

        // Algorithms flag
        {"bruteforce", no_argument, &bruteforce_flag, 1},
        {"fcount"    , no_argument, &fcount_flag, 1},
        {"fsample"   , no_argument, &fsample_flag, 1},
        {"baseline"  , no_argument, &baseline_flag, 1},

        {0, 0, 0, 0}
      };

      int option_index = 0;
      int c;
      while (1) {
        c = getopt_long(argc, argv, "g:q:p:Q:S:R:A:B:M:E:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
          case 'g':
          input_graph_flag = true;
          if (optarg != NULL) input_graph = optarg;
          break;
          case 'p':
          if (optarg != NULL) thread_count = atoi(optarg);
          break;
          case 'Q':
          if (optarg != NULL) q = atoi(optarg);
          break;
          case 'S':
          if (optarg != NULL) seed = atoi(optarg);
          break;
          case 'R':
          if (optarg != NULL) R = atoi(optarg);
          break;
          case 'A':
          if (optarg != NULL) Sa = atoi(optarg);
          break;
          case 'B':
          if (optarg != NULL) Sb = atoi(optarg);
          break;
          case 'M':
          if (optarg != NULL) mod = atoi(optarg);
          break;
          case 'E':
          if (optarg != NULL) experiment = atoi(optarg);
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
        printf("thread = %d\n", thread_count);
        printf("input_graph = %s\n", input_graph != NULL ? input_graph : "stdin");
        printf("Q = %d\n", q);
        printf("S = %d\n", q);
        printf("R = %d\n", R);
        printf("A = %d\n", Sa);
        printf("B = %d\n", Sb);
        printf("M = %d\n", mod);
        printf("E = %d\n", experiment);
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

      // if (verbose_flag) printf("|A| = %d | |B| = %d\n", Sa, Sb);

      for(unsigned int i=0; i<N; i++) sampleV.push_back(i);

      // Create DP Table
      for (unsigned int i = 0; i <= q + 1; i++) M[i] = new map<COLORSET, ll>[N + 1];

      // Random color graph
      if (verbose_flag) printf("Random coloring graph...\n");
      randomColor();

      // Fill dynamic programming table
      if (verbose_flag) printf("Processing DP table...\n");
      ll time_a = current_timestamp();
      processDP();
      ll time_b = current_timestamp() - time_a;
      if (verbose_flag) printf("End processing DP table [%llu]ms\n", time_b);

      ll time_dp = time_b;

      long long entry = 0;
      for(int i=1; i<=q; i++)
      for(int j=0; j<N; j++)
      {
        entry += M[i][j].size();
      }
      printf("DP ENTRY: [%lld]\n", entry);
      double bc_brute;
      double bc_fcount;
      double bc_fsample;
      double bc_base;

      double bc_fcount_rel;
      double bc_fsample_rel;
      double bc_base_rel;

      double fj_brute;
      double fj_fcount;
      double fj_fsample;
      double fj_base;

      double fj_fsample_rel;
      double fj_fcount_rel;
      double fj_base_rel;

      int tau_brute;
      int tau_fcount;
      int tau_base;
      int tau_fsample;

      ll time_brute = 0ll;
      ll time_fcount = 0ll;
      ll time_fsample = 0ll;
      ll time_base = 0ll;

      eng = mt19937_64(seed);
      srand(seed);

      set<int> A = randomChoose(Sa, mod);
      set<int> B = randomChoose(Sb, mod);
/*      set<int> A,B;
      int chooseA = rand() % N;
      int chooseB = -1;

      int choose = -1;
      do {
	chooseB = G[chooseA][rand() % G[chooseA].size()];
      }
      while( chooseA == chooseB || label[chooseA] != label[chooseB] );
      A.insert( chooseA );
      B.insert( chooseB );
      printf("CHOOSE[%d][%d] LABEL[%d][%d]\n", chooseA, chooseB, label[chooseA], label[chooseB]);
*/
//      A.clear();
//      A.insert(620);

//      B.clear();
//      B.insert(82);

      vector<int> X;
      for (int a : A) X.push_back(a);
      for (int b : B) X.push_back(b);
      set<int> AB;
      for (int a : A) AB.insert(a);
      for (int b : B) AB.insert(b);
      vector<int> ABv = vector<int>(AB.begin(), AB.end());

      // HEADER
      printf("Q,R,HA,HB,");
      if( bruteforce_flag ) printf("BC_BRUTE,FJ_BRUTE,TAU,TIME,");

      if( fcount_flag )
      printf("BC_FCOUNT,BC_REL_FCOUNT,FJ_FCOUNT,FJ_REL_FCOUNT,TAU_FCOUNT,TIME_FCOUNT,");

      if( fsample_flag )
      printf("BC_FSAMPLE,BC_REL_FSAMPLE,FJ_FSAMPLE,FJ_REL_FSAMPLE,TAU_FSAMPLE,TIME_FSAMPLE,");

      if( baseline_flag )
      printf("BC_BASE,BC_REL_BASE,FJ_BASE,FJ_REL_BASE,TAU_BASE,TIME_BASE,");


      printf("\n");

      map<string, ll> freqA, freqB, freqAB;
      set<string> W;
      double bcw, fjw;
      long long Rp = 0ll;
      long long Rpp = 0ll;

      double realBC = 1 , realFJ = 1;
      if( bruteforce_flag )
      {
        // BRUTE-FORCE
        Rp = 0;
        Rpp = 0;
        W.clear();
        freqA.clear();
        freqB.clear();
        dict.clear();
        freqBrute.clear();
        time_brute = current_timestamp();
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < (int)ABv.size(); i++) {
          int tid = omp_get_thread_num();
          dfs(tid, ABv[i], q - 1);
        }
        for (auto w : freqBrute) {
          int u = w.first.first;
          string s = w.first.second;
          ll freq = w.second;
          Rp += freq;
          if (A.find(u) != A.end()) {
            freqA[s] += freq;
          }
          if (B.find(u) != B.end()) {
            freqB[s] += freq;
          }
        }
        time_brute = current_timestamp() - time_brute;
        tau_brute = dict.size();
        bc_brute = realBC = bcw = BCW(dict, freqA, freqB);
        fj_brute = realFJ = fjw = FJW(dict, freqA, freqB, (long long)Rp);
      }

      // Experiments
      for(unsigned int i = 0; i < experiment; i++)
      {
        /**************************************************************************/
        /**************************************************************************/
        // FCOUNT OK
        if( fcount_flag )
        {
          // BRAY-CURTIS OK
          Rp = 0;
          W.clear();
          freqA.clear();
          freqB.clear();
          time_fcount = current_timestamp();
          set<string> Sample = randomColorfulSample(X, R);
          freqA = processFrequency(Sample, multiset<int>(A.begin(), A.end()));
          freqB = processFrequency(Sample, multiset<int>(B.begin(), B.end()));
          time_fcount = current_timestamp() - time_fcount;
          tau_fcount = Sample.size();
          bc_fcount = bcw = BCW(Sample, freqA, freqB);
          if( bruteforce_flag )bc_fcount_rel = abs(bcw - realBC) / realBC;

          // JACCARD OK
          Rp = 0;
          W.clear();
          freqA.clear();
          freqB.clear();
          Sample = randomColorfulSample(vector<int>(AB.begin(), AB.end()), R);

          for(int a : AB)
          {
            multiset<int> aa;
            aa.insert(a);
            freqAB = processFrequency(Sample, aa);
            for (auto w : freqAB) {
              Rp += w.second;
              if (A.find(a) != A.end()) freqA[w.first] += w.second;
              if (B.find(a) != B.end()) freqB[w.first] += w.second;
            }
          }
          fj_fcount = fjw = FJW(Sample, freqA, freqB, Rp);
          if( bruteforce_flag )fj_fcount_rel = abs(fjw - realFJ) / realFJ;
        }
        /**************************************************************************/
        /**************************************************************************/
        // BASELINE
        if( baseline_flag )
        {
          // BRAY-CURTIS
          Rp = 0;
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
          bc_base = bcw = BCW(W, freqA, freqB);
          if( bruteforce_flag ) bc_base_rel = abs(bcw - realBC) / realBC;

          // JACCARD
          Rp = 0;
          W.clear();
          freqA.clear();
          freqB.clear();
          BLsampling = baselineSampler(ABv, R);
          for (auto w : BLsampling) {
            int u = w.first.first;
            W.insert(w.first.second);
            Rp += w.second;
            if (A.find(u) != A.end()) freqA[w.first.second] += w.second;
            if (B.find(u) != B.end()) freqB[w.first.second] += w.second;
          }
          fj_base = fjw = FJW(W, freqA, freqB, (long long)Rp); // ); OK
          if( bruteforce_flag ) fj_base_rel = abs(fjw - realFJ) / realFJ;
        }
        /**************************************************************************/
        /**************************************************************************/
        // FSAMPLE
        if( fsample_flag )
        {
          // BRAY-CURTIS
          W.clear();
          freqA.clear();
          freqB.clear();
          Rp = 0ll;

          time_fsample = current_timestamp();
          map<pair<int, string>, ll> SamplePlus = randomColorfulSamplePlus(X, R);

          for (auto w : SamplePlus) {
            int u = w.first.first;
            W.insert(w.first.second);
            if (A.find(u) != A.end()) freqA[w.first.second] += w.second;
            if (B.find(u) != B.end()) freqB[w.first.second] += w.second;
          }

          time_fsample = current_timestamp() - time_fsample;
          tau_fsample = W.size();
          bc_fsample = bcw = BCW(W, freqA, freqB);
          if( bruteforce_flag ) bc_fsample_rel = abs(bcw - realBC) / realBC;

          // JACCARD
          W.clear();
          freqA.clear();
          freqB.clear();
          Rp = 0ll;

          SamplePlus = randomColorfulSamplePlus(ABv, R);

          for (auto w : SamplePlus) {
            int u = w.first.first;
            W.insert(w.first.second);
            Rp += w.second;
            if (A.find(u) != A.end()) {
              freqA[w.first.second] += w.second;
//              Rp += w.second;
            }
            if (B.find(u) != B.end()) {
              freqB[w.first.second] += w.second;
              // Rp += w.second;
            }
          }
          fj_fsample = fjw = FJW(W, freqA, freqB, (long long)Rp); // p);
          if( bruteforce_flag ) fj_fsample_rel = abs(fjw - realFJ) / realFJ;

        }
        /**************************************************************************/
        /**************************************************************************/

        // OUTPUT
        printf("%2d,", q);          // Q
        printf("%4d,", R);          // R
        printf("%4zu,", A.size());  // HA
        printf("%4zu,", B.size());  // HB
        if( bruteforce_flag )
        {
          printf("%.6f,", bc_brute);     // BC-BRUTE
          printf("%.6f,", fj_brute);     // FJ-BRUTE
          printf("%4d,", tau_brute);     // TAU
          printf("%4llu,", time_brute);  // TIME
        }
        if( fcount_flag )
        {
          printf("%.6f,", bc_fcount);      // BC-fcount
          printf("%.6f,", bc_fcount_rel);  // BC-fcount
          printf("%.6f,", fj_fcount);      // FJ-fcount
          printf("%.6f,", fj_fcount_rel);  // FJ-fcount
          printf("%4d,", tau_fcount);      // TAU
          printf("%4llu,", time_fcount);   // TIME
        }
        if( fsample_flag )
        {
          printf("%.6f,", bc_fsample);      // BC-fsample
          printf("%.6f,", bc_fsample_rel);  // BC-fsample
          printf("%.6f,", fj_fsample);      // FJ-fsample
          printf("%.6f,", fj_fsample_rel);  // FJ-fsample
          printf("%4d,", tau_fsample);      // TAU
          printf("%4llu,", time_fsample);   // TIME
        }
        if( baseline_flag )
        {
          printf("%.6f,", bc_base);      // BC-BASE
          printf("%.6f,", bc_base_rel);  // BC-BASE
          printf("%.6f,", fj_base);      // FJ-BASE
          printf("%.6f,", fj_base_rel);  // FJ-BASE
          printf("%4d,", tau_base);      // TAU
          printf("%4llu,", time_base);   // TIME
        }
        printf("\n");

      }

      return 0;
    }

