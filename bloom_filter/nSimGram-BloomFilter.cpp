#include <bits/stdc++.h>
#include <omp.h>
#include "cxxopts.hpp"

#define ERROR(c,s) if(c){perror(s); return -1;}

typedef uint64_t bf_t;
typedef uint64_t ll;
typedef std::vector<int> path;
typedef std::string qpath;
typedef std::set<qpath> dict_t;     // dictionary of string (q-paths)
typedef std::map<qpath, ll> fdict_t;// frequency dictionary of string (q-paths)

// Bloom filter
namespace bf
{
  bf_t create(int z, int h, std::mt19937& rng)
  {
    bf_t out = bf_t(0);
    for(int i=0; i<h; i++) out |= 1<<(rng()%z);
    /*
    	maybe:
    		bitset<> {1, ..., 1, 0, ..., 0}
    				  1 x h + 0 x (z-h)
    		random shuffle
    		bit -> int
    */
    return out;
  }

  bf_t join(const bf_t a, const bf_t b)
  {
    return (a|b);
  }
}

/*
struct bloom_filter 
{
  bf_t data;
	
  bloom_filter(int z, int h, std::mt19937& rng)
  {
    std::vector<bool> vb(z, false);
    std::fill(vb.begin(), vb.begin()+h, true);
    std::shuffle(vb.begin(), vb.end(), rng);
    data = bf_t(0);
    for(bool b : vb) data = (data<<1) | b;
  }
  
  bloom_filter(bf_t bdata) : data(bdata) {}
  
  struct bloom_filter operator+(const struct bloom_filter& a) const
  {
    return struct bloom_filter(data|a.data);
  }

  // equality comparison. doesn't modify object. therefore const.
  bool operator==(const struct bloom_filter& a) const
  {
    return (data == a.data);
  }  
};
*/

// Graph data struct
typedef struct graph 
{
  int N;

  std::vector<char> label;
  std::vector<bf_t> filter;
  std::vector<std::vector<int>> edges;

  graph(int n) 
  {
    N = n;
    label.resize(N, 0);
    filter.resize(N, 0);
    edges.resize(N, std::vector<int>());
  };

  void addEdge(int i, int j) 
  {
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

bool first = true;

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
  std::sort(out.begin(), out.end());
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
std::vector< std::vector<std::map<bf_t, long long int>> > dp;
void processDP()
{

  // create matrix (Q+1)*N
  dp.resize(Q+1, std::vector<std::map<bf_t, long long int>>(N));

  auto timer_now = timer_start();

  // Base case
  std::cerr << "DP preprocessing 1/" << Q << std::endl;
  #pragma omp parallel for schedule(guided)
  for(int u = 0; u < N; u++) dp[1][u][G.filter[u]] = 1;

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
        for (auto [s, f] : dp[i-1][v])
        {
          bf_t s2 = G.filter[u];
          if(bf::join(s, s2) == s) continue;
          dp[i][u][bf::join(s, s2)] = f+1;
        }
      }
    }
    ll count = 0;
    for(int u=0; u<N; u++) count += dp[i][u].size();
    std::cerr << "\t" << count << " bf at level " << i << std::endl;
  }
  std::cerr << "DP table processed in " << timer_step(timer_now) << "ms" << std::endl;

}


// Fcount
bool isPrefix(dict_t& W, std::string& x) {
  auto it = W.lower_bound(x);
  if (it == W.end()) return false;
  return std::mismatch(x.begin(), x.end(), (*it).begin()).first == x.end();
}

path randomPathTo(int u)
{
  path p = {u};
    
  for(size_t i=1; i<Q; i++)
  {
    std::vector<ll> freq;
    bf_t pref(0);
    for(size_t j=1; j<p.size(); j++) pref = bf::join(pref, G.filter[p[j]]);

    for (int v : G.edges[u])
    {
      bf_t temp  = bf::join(pref, G.filter[v]);
      bf_t temp2 = bf::join(temp, G.filter[p[0]]);
      if(temp == temp2) continue;
      freq.push_back(dp[i][v][temp]);
    }
    std::discrete_distribution<ll> distribution(freq.begin(), freq.end());

    #pragma omp critical
    {
      u = G.edges[u][distribution(rng)];
    }
    p.push_back(u);

  }
  return p;
}

dict_t randomSample()
{
  dict_t W;

  ll fA = static_cast<ll>(dp[Q][A].size());
  ll fB = static_cast<ll>(dp[Q][B].size());
  std::vector<ll> fV = {fA, fB};
  
  std::discrete_distribution<ll> distribution(fV.begin(), fV.end());
  std::vector<int> X = {A,B};

  std::set<path> R;
  while(R.size() < Rsize)
  {
    int start = X[distribution(rng)];
    path toAdd = randomPathTo(start);
    if(toAdd.size() < Q) continue;
    R.insert(toAdd);
  }

  for(const path& r : R) W.insert(L(r));

  return W;
}

fdict_t processFrequency(dict_t& W, int X)
{
  fdict_t ret;
   
  return ret;
}

std::tuple<double, double> fcount()
{
  double fj = 0.;
  double bc = 0.;

  dict_t W = randomSample();

  fdict_t freqA = processFrequency(W, A);
  fdict_t freqB = processFrequency(W, B);

  ll Rcount = 0;
  for(const auto& [w, v] : freqA) Rcount += v;
  for(const auto& [w, v] : freqB) Rcount += v;

  fj = fjw(W, freqA, freqB, Rcount);
  bc = bcw(W, freqA, freqB);

  return std::make_tuple(fj, bc);
}

// Fsample
std::tuple<dict_t, fdict_t, fdict_t> randomSamplePlus()
{
  dict_t W;
  fdict_t freqA;
  fdict_t freqB;

  ll fA = dp[Q][A].size();
  ll fB = dp[Q][B].size();
  std::vector<ll> fV = {fA, fB};
  std::discrete_distribution<ll> distribution(fV.begin(), fV.end());
  std::vector<int> X = {A,B};

  std::set<path> R;
  while(R.size() < Rsize)
  {
    int start = X[distribution(rng)];
    path toAdd = randomPathTo(start);
    if(toAdd.size() < Q) continue;
    R.insert(toAdd);
  }

  for(const path& p : R)
  {
    qpath q = L(p);
    W.insert(q);
    if(p[0] == A) freqA[q]++;
    else freqB[q]++;
  }

  return std::make_tuple(W, freqA, freqB);
}

std::tuple<double, double> fsample()
{
  double fj = 0.;
  double bc = 0.;

  std::tuple<dict_t, fdict_t, fdict_t> sample = randomSamplePlus();

  dict_t W = std::get<0>(sample);

  fdict_t freqA = std::get<1>(sample);
  fdict_t freqB = std::get<2>(sample);

  ll Rsample = 0;
  for(const auto& [w, v] : freqA) Rsample += v;
  for(const auto& [w, v] : freqB) Rsample += v;

  fj = fjw(W, freqA, freqB, Rsample);
  bc = bcw(W, freqA, freqB);

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
  for(int i=0; i<N; i++) G.filter[i] = bf::create(Z, H, rng);
  std::cerr << "end" << std::endl;

  // Bruteforce similarity
  double real_fj = 0.;
  double real_bc = 0.;

  if(bruteforce_f)
  {
    std::cerr << "Bruteforce..." << std::endl;
    auto bf = bruteforce();
    real_fj = std::get<0>(bf);
    real_bc = std::get<1>(bf);
    std::cerr << "End bruteforce" << std::endl;
    std::cerr << "Real FJ(A,B) = " << real_fj << std::endl;
    std::cerr << "Real BC(A,B) = " << real_bc << std::endl;
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
  std::vector<double> fcount_fj;
  std::vector<double> fcount_bc;
  std::vector<double> fsample_fj;
  std::vector<double> fsample_bc;
  std::vector<double> baseline_fj;
  std::vector<double> baseline_bc;

  // Start experiments
  for(int i=0; i<experiments; i++)
  {
    if(bruteforce_f)
    {
      std::cout << real_fj << ",";
      std::cout << real_bc << ",";
    }

    if(fcount_f)
    {
      auto fc = fcount();
      std::cout << std::get<0>(fc) << ",";
      std::cout << std::get<1>(fc) << ",";
      fcount_fj.push_back(static_cast<double>(std::get<0>(fc)));
      fcount_bc.push_back(static_cast<double>(std::get<1>(fc)));
    }

    if(fsample_f)
    {
      auto fs = fsample();
      std::cout << std::get<0>(fs) << ",";
      std::cout << std::get<1>(fs) << ",";
      fsample_fj.push_back(static_cast<double>(std::get<0>(fs)));
      fsample_bc.push_back(static_cast<double>(std::get<1>(fs)));
    }

    if(baseline_f)
    {
      auto bl = baseline();
      std::cout << std::get<0>(bl) << ",";
      std::cout << std::get<1>(bl) << ",";
      baseline_fj.push_back(static_cast<double>(std::get<0>(bl)));
      baseline_bc.push_back(static_cast<double>(std::get<1>(bl)));
    }
    std::cout << std::endl;
  }

  auto avg = [](const std::vector<double> &V)
  {
    if(V.size() == 0) return 0.;
    return std::accumulate(V.begin(), V.end(), 0.) / V.size();
  };

  auto var = [](const std::vector<double> &V, double avg)
  {
    if(V.size() == 0) return 0.;
    return std::accumulate(V.begin(), V.end(), 0., [avg](const double sum, const double x)
    {
      return sum + (x-avg)*(x-avg);
    }) / V.size();
  };

  double   fcount_fj_avg = avg(fcount_fj);
  double   fcount_bc_avg = avg(fcount_bc);
  double  fsample_fj_avg = avg(fsample_fj);
  double  fsample_bc_avg = avg(fsample_bc);
  double baseline_fj_avg = avg(baseline_fj);
  double baseline_bc_avg = avg(baseline_bc);

  double   fcount_fj_var = var(fcount_fj, fcount_fj_avg);
  double   fcount_bc_var = var(fcount_bc, fcount_bc_avg);
  double  fsample_fj_var = var(fsample_fj, fsample_fj_avg);
  double  fsample_bc_var = var(fsample_bc, fsample_bc_avg);
  double baseline_fj_var = var(baseline_fj, baseline_fj_avg);
  double baseline_bc_var = var(baseline_bc, baseline_bc_avg);

  std::cout << std::endl;
  std::cout << "algorithm,fj_average,fj_variance,bc_average,bc_variance" << std::endl;
  
  if(bruteforce_f) std::cout << "bruteforce," <<     real_fj     << "," <<      "0.000000" << "," <<     real_bc     << "," <<      "0.000000" << std::endl;
  if(    fcount_f) std::cout << "    fcount," <<   fcount_fj_avg << "," <<   fcount_fj_var << "," <<   fcount_bc_avg << "," <<   fcount_bc_var << std::endl;
  if(   fsample_f) std::cout << "   fsample," <<  fsample_fj_avg << "," <<  fsample_fj_var << "," <<  fsample_bc_avg << "," <<  fsample_bc_var << std::endl;
  if(  baseline_f) std::cout << "  baseline," << baseline_fj_avg << "," << baseline_fj_var << "," << baseline_bc_avg << "," << baseline_bc_var << std::endl;

  // Flush output buffers
  std::cout.flush();
  std::cerr.flush();

  return 0;
}
