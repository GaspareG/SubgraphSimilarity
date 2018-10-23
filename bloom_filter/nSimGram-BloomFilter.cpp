#include <bits/stdc++.h>
#include <omp.h>
#include "cxxopts.hpp"
#include "bloomfilter.hpp"

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

// Bray Curtis weighted
double BCW(const dict_t& W, fdict_t& freqA, fdict_t& freqB) {
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
double FJW(const dict_t& W, fdict_t& freqA, fdict_t& freqB, ll R) {
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
std::vector<path> dfs(int source)
{
  std::vector<path> out;
  // do things
  return out;
}

std::tuple<double, double> bruteforce()
{
  double fj = 0.;
  double bc = 0.;

  fdict_t freqA, freqB;
  dict_t W;

  std::vector<path> pathA = dfs(A); // generate all paths from A
  std::vector<path> pathB = dfs(B); // generate all paths from B

  for(const path& w : pathA) freqA[L(w)]++; // compute qpaths frequencies in A
  for(const path& w : pathB) freqB[L(w)]++; // compute qpaths frequencies in B

  for(const auto& [w, v] : freqA) W.insert(w); // build dictionary W from qpaths in A
  for(const auto& [w, v] : freqB) W.insert(w); // build dictionary W from qpaths in B

  fj = FJW(W, freqA, freqB, Rsize);
  bc = BCW(W, freqA, freqB);

  return std::make_tuple(fj, bc);
}

// DP Preprocessing
void processDP()
{

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
    int next = -1;

    seen.insert(last);
    while(next < 0)
    {
        // G.edges[last];
        // TODO
    }

    rPath.push_back(last);
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

  fj = FJW(W, freqA, freqB, Rsize);
  bc = BCW(W, freqA, freqB);

  return std::make_tuple(fj, bc);
}

int main(int argc, char** argv) {

  // Parse arguments
  cxxopts::Options options("nSimGram-BloomFilter","Compute node similarity using bloom filter");
  options.add_options()
    (       "h,help", "Print help",                                           cxxopts::value(help))
    (    "v,verbose", "Verbose log",                                          cxxopts::value(verbose))
    (      "i,input", "Input file name (default: stdin)",                     cxxopts::value(input))
    (     "o,output", "Output file name (default: stdout)",                   cxxopts::value(output))
    (   "bruteforce", "Compute similarity with bruteforce",                   cxxopts::value(bruteforce_f))
    (       "fcount", "Compute similarity with fCount",                       cxxopts::value(fcount_f))
    (      "fsample", "Compute similarity with fSample",                      cxxopts::value(fsample_f))
    (     "baseline", "Compute similarity with baseline",                     cxxopts::value(baseline_f))
    (          "all", "Compute similarity with all algorithms",               cxxopts::value(all))
    ("e,experiments", "Number of experiments to run (default: 1)",            cxxopts::value(experiments))
    (      "r,rsize", "Size of the sample",                                   cxxopts::value(Rsize))
    (    "t,threads", "Number of parallel threads to run (default: all)",     cxxopts::value(Nthreads))
    (      "f,first", "Drop first label in path (default: false)",            cxxopts::value(first))
    (            "Q", "Length of the paths (default: 4)",                     cxxopts::value(Q))
    (            "A", "First node of similarity (default: random node)",      cxxopts::value(A))
    (            "B", "Second node of similarity (default: random node)",     cxxopts::value(B))
    (            "H", "Number of hash function in bloom filter (default: 8)", cxxopts::value(H))
    (            "Z", "Number of bits in bloom filter (default: 64)",         cxxopts::value(Z))
    (            "S", "Seed for random number generation (default: 42)",      cxxopts::value(seed));

  auto result = options.parse(argc, argv);

  if (help || !(bruteforce_f || fcount_f || fsample_f || baseline_f || all)) {
    std::cout << options.help();
    return 0;
  }

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

  // Read input
  std::cin >> N >> M;
  std::cerr << "Read graph N = " << N << " M = " << M << std::endl;

  G = graph(N);

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
  for(int i=0; i<N; i++)
    G.filter[i] = 0; // TODO createBloomFilter(G.label[i], Z, H, rng);
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

  // Print csv header
  if(bruteforce_f)
  {
    std::cout << "bruteforce_FJ" << "," << std::endl;
    std::cout << "bruteforce_BC" << "," << std::endl;
  }

  if(fcount_f)
  {
    std::cout << "fcount_FJ" << "," << std::endl;
    std::cout << "fcount_BC" << "," << std::endl;
  }

  if(fsample_f)
  {
    std::cout << "fsample_FJ" << "," << std::endl;
    std::cout << "fsample_BC" << "," << std::endl;
  }

  if(baseline_f)
  {
    std::cout << "baseline_FJ" << "," << std::endl;
    std::cout << "baseline_BC" << "," << std::endl;
  }

  // Start experiments
  for(int i=0; i<experiments; i++)
  {
    if(bruteforce_f)
    {
      std::cout << realFJ << "," << std::endl;
      std::cout << realBC << "," << std::endl;
    }

    if(fcount_f)
    {
      auto fc = fcount();
      std::cout << std::get<0>(fc) << "," << std::endl;
      std::cout << std::get<1>(fc) << "," << std::endl;
    }

    if(fsample_f)
    {
      auto fs = fsample();
      std::cout << std::get<0>(fs) << "," << std::endl;
      std::cout << std::get<1>(fs) << "," << std::endl;
    }

    if(baseline_f)
    {
      auto bl = baseline();
      std::cout << std::get<0>(bl) << "," << std::endl;
      std::cout << std::get<1>(bl) << "," << std::endl;
    }
  }

  // Flush output buffers
  std::cout.flush();
  std::cerr.flush();

  return 0;
}
