#include <bits/stdc++.h>
#include <omp.h>
#include "cxxopts.hpp"
#include "bloomfilter.hpp"

#define ERROR(c,s) if( c ){ perror(s); return -1; }

typedef uint64_t ll;
typedef std::vector<int> path;
typedef std::string qpath;
typedef std::set<qpath> dict_t;     // dictionary of string (q-paths)
typedef std::map<qpath, ll> fdict_t;// frequency dictionary of string (q-paths)

// Graph data struct
typedef struct graph {
  int N;

  std::vector<char> label;
  std::vector<bf_t> filter;
  std::vector<std::vector<int>> edges;

  graph(int n) {
    N = n;
    label.resize(N, 0);
    filter.resize(N, 0);
    edges.resize(N, std::vector<int>());
  };

  void addEdge(int i, int j) {
    edges[i].push_back(j);
    edges[j].push_back(i);
  };

} graph;

// Graph
graph G(1);
int N, M;

std::mt19937 rng;

// Execution flag -> refactor made struct
bool verbose = false;
bool help = false;

bool bruteforce_f = false;
bool fcount_f = false;
bool fsample_f = false;
bool baseline_f = false;
bool all = false;

std::string input = "";
std::string output = "";

int experiments = 1;
size_t Rsize = 1000;
int Nthreads = -1;

bool first = false;

size_t Q = 4;
int A = -1;
int B = -1;
int H = 8;
int Z = 64;
int seed = 42;

// Time functions
auto timer_start()
{
  return std::chrono::high_resolution_clock::now();
}

auto timer_step(auto timer_now)
{
  auto elapsed = timer_start() - timer_now;
  auto msec    = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  return msec;
}

// Bray Curtis weighted
double bcw(const dict_t& W, fdict_t& freqA, fdict_t& freqB) {
  double ret = 0.;
  for (const qpath& x : W)
  {
    double fax = static_cast<double>(freqA[x]);
    double fbx = static_cast<double>(freqB[x]);
    if (fax == 0. && fbx == 0.) continue; // should not happens
    ret += std::min(fax, fbx) / (fax + fbx);
  }
  ret *= 2.;
  ret /= static_cast<double>(W.size());
  return ret;
}

// Frequency Jaccard Weighted
double fjw(const dict_t& W, fdict_t& freqA, fdict_t& freqB, ll R) {
  ll num = 0ll;
  for (qpath x : W) num += std::min(freqA[x], freqB[x]);
  return static_cast<double>(num) / static_cast<double>(R);
}

// path -> qpath
qpath L(const path& p)
{
  qpath out;
  out.reserve(p.size());
  for(int v : p) out += G.label[v];
  if(first) out = out.substr(1);
  return out;
}

// Bruteforce
std::vector<path> bfs(int source)
{
  std::vector<path> out;
  std::deque<path> pqueue;

  path start = {source};
  pqueue.push_back(start);

  while(!pqueue.empty())
  {
    path att = pqueue.front();
    pqueue.pop_front();

    if(att.size() == Q)
    {
      out.push_back(att);
      continue;
    }

    for(int next : G.edges[att[att.size()-1]])
    {
      if(std::find(std::begin(att), std::end(att), next) != std::end(att)) continue;
      att.push_back(next);
      pqueue.push_back(att);
      att.pop_back();
    }
  }
  return out;
}

std::tuple<double, double> bruteforce()
{
  double fj = 0.;
  double bc = 0.;

  fdict_t freqA, freqB;
  dict_t W;

  auto start = timer_start();

  std::vector<path> pathA = bfs(A); // generate all paths from A
  std::vector<path> pathB = bfs(B); // generate all paths from B

  std::cerr << "Found " << pathA.size() << " " << Q << "-paths from " << A << std::endl;
  std::cerr << "Found " << pathB.size() << " " << Q << "-paths from " << B << std::endl;

  for(const path& w : pathA) freqA[L(w)]++; // compute qpaths frequencies in A
  for(const path& w : pathB) freqB[L(w)]++; // compute qpaths frequencies in B

  for(const auto& [w, v] : freqA) W.insert(w); // build dictionary W from qpaths in A
  for(const auto& [w, v] : freqB) W.insert(w); // build dictionary W from qpaths in B

  ll Rbruteforce = 0;
  for(const auto& [w, v] : freqA) Rbruteforce += v;
  for(const auto& [w, v] : freqB) Rbruteforce += v;

  fj = fjw(W, freqA, freqB, Rbruteforce);
  bc = bcw(W, freqA, freqB);

  std::cerr << "Bruteforce similarity in " << timer_step(start) << "msec" << std::endl;

  return std::make_tuple(fj, bc);
}

// DP Preprocessing
std::vector<std::vector<std::map<bf_t, ll>>> dp;

void processDP()
{

  // create matrix (Q+1)*N
  dp.resize(Q+1, std::vector<std::map<bf_t, ll>>(N));

  auto timer_now = timer_start();

  // Base case
  std::cerr << "DP preprocessing 1/" << Q << std::endl;
  #pragma omp parallel for schedule(guided)
  for(int u = 0; u < N; u++) dp[1][u][G.filter[u]] = 1ll;

  // For each level
  for(size_t i = 2; i <= Q; i++)
  {
    std::cerr << "DP preprocessing " << i << "/" << Q << " (" << timer_step(timer_now) << "ms)" << std::endl;
    // For each node
    #pragma omp parallel for schedule(guided)
    for (int u = 0; u < N; u++)
    {
      // For each neighbor
      for (int v : G.edges[u])
      {
        // For each entry in dp table
        for (auto d : dp[i-1][v])
        {
          bf_t s = d.first;
          bf_t s2 = G.filter[u];
          ll f = d.second;
          if((s | s2) == s) continue;
          dp[i][u][s|s2] += f;
        }
      }
    }
  }
  std::cerr << "DP table processed in " << timer_step(timer_now) << "ms" << std::endl;

}


// Fcount
std::tuple<double, double> fcount()
{
  double fj = 0.;
  double bc = 0.;

  // TODO

  return std::make_tuple(fj, bc);
}

// Fsample
std::tuple<double, double> fsample()
{
  double fj = 0.;
  double bc = 0.;

  // TODO

  return std::make_tuple(fj, bc);
}

// Baseline
path naiveRandomPath()
{
  path rPath = { (rng()%2==0)?A:B };
  std::set<int> seen;

  while(rPath.size() < Q)
  {
    int last = rPath[rPath.size()-1];
    seen.insert(last);

    std::vector<int> valid;
    for(int n : G.edges[last])
      if(seen.find(n) == seen.end())
        valid.push_back(n);

    if(valid.size() == 0) break;

    rPath.push_back(valid[rng()%valid.size()]);
  }

  return rPath;
}

std::tuple<double, double> baseline()
{
  double fj = 0.;
  double bc = 0.;

  fdict_t freqA;
  fdict_t freqB;
  dict_t W;
  std::set<path> r;

  while(r.size() < Rsize)
  {
    path randomPath = naiveRandomPath();
    bool fromA = *std::begin(randomPath) == A;

    if(randomPath.size() < Q) continue;
    if(r.find(randomPath) != r.end()) continue;
    r.insert(randomPath);

    if(fromA) freqA[L(randomPath)]++;
    else freqB[L(randomPath)]++;
  }

  for(const auto& [w, v] : freqA) W.insert(w);
  for(const auto& [w, v] : freqB) W.insert(w);

  ll Rbaseline= 0;
  for(const auto& [w, v] : freqA) Rbaseline += v;
  for(const auto& [w, v] : freqB) Rbaseline += v;

  fj = fjw(W, freqA, freqB, Rsize);
  bc = bcw(W, freqA, freqB);

  return std::make_tuple(fj, bc);
}

int main(int argc, char** argv) {


  // Parse arguments
  cxxopts::Options options("nSimGram-BloomFilter","Compute node similarity using bloom filter");
  options.add_options()
    // Help flag
    (       "h,help", "Print help",                                           cxxopts::value(help))
    // Stream options
    (    "v,verbose", "Verbose log",                                          cxxopts::value(verbose))
    (      "i,input", "Input file name (default: stdin)",                     cxxopts::value(input))
    (     "o,output", "Output file name (default: stdout)",                   cxxopts::value(output))
    // Algorithms to run
    (   "bruteforce", "Compute similarity with bruteforce",                   cxxopts::value(bruteforce_f))
    (       "fcount", "Compute similarity with fCount",                       cxxopts::value(fcount_f))
    (      "fsample", "Compute similarity with fSample",                      cxxopts::value(fsample_f))
    (     "baseline", "Compute similarity with baseline",                     cxxopts::value(baseline_f))
    (          "all", "Compute similarity with all algorithms",               cxxopts::value(all))
    // Experiments parameters
    ("e,experiments", "Number of experiments to run (default: 1)",            cxxopts::value(experiments))
    (      "r,rsize", "Size of the sample",                                   cxxopts::value(Rsize))
    (    "t,threads", "Number of parallel threads to run (default: all)",     cxxopts::value(Nthreads))
    (      "f,first", "Drop first label in path (default: false)",            cxxopts::value(first))
    // Algorithms parameters
    (            "Q", "Length of the paths (default: 4)",                     cxxopts::value(Q))
    (            "A", "First node of similarity (default: random node)",      cxxopts::value(A))
    (            "B", "Second node of similarity (default: random node)",     cxxopts::value(B))
    (            "H", "Number of hash function in bloom filter (default: 8)", cxxopts::value(H))
    (            "Z", "Number of bits in bloom filter (default: 64)",         cxxopts::value(Z))
    (            "S", "Seed for random number generation (default: 42)",      cxxopts::value(seed));

  auto result = options.parse(argc, argv);

  // Print help
  if (help || !(bruteforce_f || fcount_f || fsample_f || baseline_f || all)) {
    std::cout << options.help();
    return 0;
  }

  // Check parameters validity
  ERROR(experiments < 1,"Number of experiments too low");
  ERROR(Rsize < 1,"Size of sample too low");
  ERROR(Nthreads < 1,"Number of threads too low");
  ERROR(Nthreads > omp_get_max_threads(),"Number of threads too high");
  ERROR(Q < 2,"Length of path too low");
  ERROR(Q > 16,"Length of path too high");
  ERROR(A < -1,"Invalid node A");
  ERROR(B < -1,"Invalid node B");
  ERROR(H < 1,"Number of hash functions too low");
  ERROR(Z < 1, "Number of bits in bloom filter too low");
  ERROR(H >= Z, "Too many hash functions (H >= Z)");

  // Set number of threads
  omp_set_num_threads(Nthreads);

  // Set all the experiments
  if (all) {
    bruteforce_f = true;
    fcount_f = true;
    fsample_f = true;
    baseline_f = true;
  }

  // Random number generator
  rng = std::mt19937(seed);

  // Redirect cerr buffer
  if (!verbose) {
    std::ofstream fnull("/dev/null");
    std::cerr.rdbuf(fnull.rdbuf());
  }

  // Redirect cin buffer
  if (input.size() != 0) {
    std::ifstream fin(input);
    std::cin.rdbuf(fin.rdbuf());
  }

  // Redirect cout buffer
  if (output.size() != 0) {
    std::ofstream fout(output);
    std::cout.rdbuf(fout.rdbuf());
  }

  // Output precision
  std::cout.setf(std::ios::fixed, std:: ios::floatfield);
  std::cerr.setf(std::ios::fixed, std:: ios::floatfield);

  std::cout.precision(6);
  std::cerr.precision(6);

  // Read input
  std::cin >> N >> M;
  std::cerr << "Read graph N = " << N << " M = " << M << std::endl;

  G = graph(N);

  // Check starting node
  ERROR(A >= N, "Invalid node A");
  ERROR(B >= N, "Invalid node B");

  if(A == -1) while(A == B) A = rng()%N;
  if(B == -1) while(A == B) B = rng()%N;

  // Reading labels
  std::cerr << "Reading labels..." << std::endl;
  for (int i = 0; i < N; i++) std::cin >> G.label[i];
  for (int i = 0; i < N; i++) G.label[i] += 33; // make printable, maybe change?
  std::cerr << "end" << std::endl;

  // Reading nodes
  std::cerr << "Reading edges..." << std::endl;
  for (int i = 0; i < M; i++) {
    int a, b;
    std::cin >> a >> b;
    G.addEdge(a, b);
  }
  std::cerr << "end" << std::endl;

  // Create filter
  std::cerr << "Create bloomfilter..." << std::endl;
  for(int i=0; i<N; i++) G.filter[i] = 0; // TODO createBloomFilter(G.label[i], Z, H, rng);
  std::cerr << "end" << std::endl;

  // Bruteforce similarity
  double realFJ = 0.;
  double realBC = 0.;

  if(bruteforce_f)
  {
    std::cerr << "Bruteforce..." << std::endl;
    auto bf = bruteforce();
    realFJ = std::get<0>(bf);
    realBC = std::get<1>(bf);
    std::cerr << "End bruteforce" << std::endl;
    std::cerr << "Real FJ(A,B) = " << realFJ << std::endl;
    std::cerr << "Real BC(A,B) = " << realBC << std::endl;
  }

  // Process DP only if fcount or fsample are enabled
  if(fcount_f || fsample_f)
  {
    std::cerr << "Start processing DP Table..." << std::endl;
    processDP();
    std::cerr << "end" << std::endl;
  }

  // Print CSV headers
  if(bruteforce_f)
  {
    std::cout << "bruteforce_FJ" << ",";
    std::cout << "bruteforce_BC" << ",";
  }

  if(fcount_f)
  {
    std::cout << "fcount_FJ" << ",";
    std::cout << "fcount_BC" << ",";
  }

  if(fsample_f)
  {
    std::cout << "fsample_FJ" << ",";
    std::cout << "fsample_BC" << ",";
  }

  if(baseline_f)
  {
    std::cout << "baseline_FJ" << ",";
    std::cout << "baseline_BC" << ",";
  }

  std::cout << std::endl;

  // Average values for statistics
  double fcount_fj_avg = 0.;
  double fcount_bc_avg = 0.;
  double fsample_fj_avg = 0.;
  double fsample_bc_avg = 0.;
  double baseline_fj_avg = 0.;
  double baseline_bc_avg = 0.;

  // Start experiments
  for(int i=0; i<experiments; i++)
  {
    if(bruteforce_f)
    {
      std::cout << realFJ << ",";
      std::cout << realBC << ",";
    }

    if(fcount_f)
    {
      auto fc = fcount();
      std::cout << std::get<0>(fc) << ",";
      std::cout << std::get<1>(fc) << ",";
      fcount_fj_avg += static_cast<double>(std::get<0>(fc));
      fcount_bc_avg += static_cast<double>(std::get<1>(fc));
    }

    if(fsample_f)
    {
      auto fs = fsample();
      std::cout << std::get<0>(fs) << ",";
      std::cout << std::get<1>(fs) << ",";
      fsample_fj_avg += static_cast<double>(std::get<0>(fs));
      fsample_bc_avg += static_cast<double>(std::get<1>(fs));
    }

    if(baseline_f)
    {
      auto bl = baseline();
      std::cout << std::get<0>(bl) << ",";
      std::cout << std::get<1>(bl) << ",";
      baseline_fj_avg += static_cast<double>(std::get<0>(bl));
      baseline_bc_avg += static_cast<double>(std::get<1>(bl));
    }
    std::cout << std::endl;
  }

  fcount_fj_avg /= static_cast<double>(experiments);
  fcount_bc_avg /= static_cast<double>(experiments);
  fsample_fj_avg /= static_cast<double>(experiments);
  fsample_bc_avg /= static_cast<double>(experiments);
  baseline_fj_avg /= static_cast<double>(experiments);
  baseline_bc_avg /= static_cast<double>(experiments);

  if(fcount_f || fsample_f || baseline_f)
  {
    std::cerr << std::endl;
    std::cerr << "Average of " << experiments << " experiments:" << std::endl;
    if(bruteforce_f) std::cerr << "Bruteforce FJ(A,B) = " <<          realFJ << " BC(A,B) = " <<          realBC << std::endl;
    if(    fcount_f) std::cerr << "FCount:    FJ(A,B) = " <<   fcount_fj_avg << " BC(A,B) = " <<   fcount_bc_avg << std::endl;
    if(   fsample_f) std::cerr << "FSample:   FJ(A,B) = " <<  fsample_fj_avg << " BC(A,B) = " <<  fsample_bc_avg << std::endl;
    if(  baseline_f) std::cerr << "Baseline:  FJ(A,B) = " << baseline_fj_avg << " BC(A,B) = " << baseline_bc_avg << std::endl;
  }

  // Flush output buffers
  std::cout.flush();
  std::cerr.flush();

  return 0;
}
