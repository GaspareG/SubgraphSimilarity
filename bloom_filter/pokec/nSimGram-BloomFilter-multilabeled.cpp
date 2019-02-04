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

typedef struct bloom_filter
{
  bf_t data;

  bloom_filter(){ data = 0; }
  bloom_filter(const bf_t d){ data = d; }
  bloom_filter(const struct bloom_filter &d){ data = d.data; }
  bloom_filter(const int z, const int h, std::mt19937& rng)
  {
    std::vector<bool> vb(z, false);
    std::fill(vb.begin(), vb.begin()+h, true);
    std::shuffle(vb.begin(), vb.end(), rng);
    data = bf_t(0);
    for(bool b : vb) data = (data<<1) | b;
  }

  bloom_filter operator+(const struct bloom_filter& a) const
  {
    return bloom_filter(data|a.data);
  }

  // equality comparison. doesn't modify object. therefore const.
  bool operator==(const struct bloom_filter& a) const
  {
    return (data == a.data);
  }

  // equality comparison. doesn't modify object. therefore const.
  bool operator<(const struct bloom_filter& a) const
  {
    return (data < a.data);
  }

} bloom_filter;

// Graph data struct
typedef struct graph
{
  int N;

  std::vector<char> label;
  std::vector<bloom_filter> filter;
  std::vector<std::vector<int>> edges;
  std::vector<std::multiset<int>> attributes;
  std::vector<std::string> attnames;

  graph(int n)
  {
    N = n;
    label.resize(N, 0);
    filter.resize(N, bloom_filter());
    edges.resize(N, std::vector<int>());
    attributes.resize(N, std::multiset<int>());
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
bool sort = true;
bool multi = false;

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
  if(sort) std::sort(out.begin(), out.end());
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

std::multiset<int> multiUnion(const std::multiset<int> &m1, const std::multiset<int> &m2)
{
  std::multiset<int> ret;
  for(auto v = m1.begin(); v != m1.end(); v=m1.upper_bound(*v))
  {
    int c = multi ? m1.count(*v) : 1;
    while(c--) ret.insert(*v);
  }
  for(auto v = m2.begin(); v != m2.end(); v=m2.upper_bound(*v))
  {
    if(ret.count(*v) > 0) continue;
    int c = multi ? m2.count(*v) : 1;
    if(multi && ret.count(*v) > 0) c=0;
    while(c--) ret.insert(*v);
  }
  return ret;
}

std::multiset<int> multiIntersect(const std::multiset<int> &m1, const std::multiset<int> &m2)
{
  std::multiset<int> ret;
  for(auto v = m1.begin(); v != m1.end(); v=m1.upper_bound(*v))
  {
    int c = std::min(m1.count(*v), m2.count(*v));
    if(multi) c = c > 0 ? 1 : 0;
    while(c--) ret.insert(*v);
  }
  return ret;
}

std::multiset<int> fingerprint(const std::vector<path> &paths)
{
  auto pathToAttr = [](const path& p){
    if(p.size() == 1) return std::multiset<int>();
    std::multiset<int> attr = G.attributes[p[1]];
    for(size_t i=2; i<p.size(); i++) attr = multiIntersect(attr, G.attributes[p[i]]);
    return attr;
  };

  auto AttrToFinger = [](const std::vector<std::multiset<int>> &fingers){
    if(fingers.size() == 0) return std::multiset<int>();
    std::multiset<int> ret = fingers[0];
    for(size_t i=1; i<fingers.size(); i++)
    {
      ret = multiUnion(ret, fingers[i]);
    }
    return ret;
  };

  std::vector<std::multiset<int>> fingers;
  for(auto &p : paths) fingers.push_back(pathToAttr(p));
  return AttrToFinger(fingers);
}

// DP Preprocessing
std::vector< std::vector<std::map<bloom_filter, long long int>> > dp;
void processDP()
{
  // create matrix (Q+1)*N
  dp.resize(Q+1, std::vector<std::map<bloom_filter, long long int>>(N));

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
        bloom_filter s2 = G.filter[u];
        // For each entry in dp table
        for (auto [s, f] : dp[i-1][v])
        {
          if((s+s2) == s) continue;
          dp[i][u][s+s2] = f+1;
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
bool isPrefix(dict_t& W, qpath& x) {

  if(!sort)
  {
    auto it = W.lower_bound(x);
    if (it == W.end()) return false;
    else return std::mismatch(x.begin(), x.end(), (*it).begin()).first == x.end();
  }

  bool found = false;
  for(const qpath &q : W)
  {
    size_t l=0;
    size_t r=0;
    while(l < q.size() && r < x.size())
    {
      if(q[l] == x[r]) l++, r++;
      else if(q[l] < x[r]) l++;
      else break;
    }

    found = (r == x.size());
     if(found) break;
  }

  return found;
}

path randomPathTo(int u)
{
  path p = {u};
  bloom_filter cur = G.filter[u];

  for(size_t i=2; i<=Q; i++)
  {
    std::vector<ll> freq;
    bool one = false;
    for (int v : G.edges[u])
    {
      bool valid = !((cur+G.filter[v]) == cur);
      freq.push_back(valid ? dp[i][v][cur+G.filter[v]] : 0);
      one = one || valid;
    }

    if(!one) return p;

    std::discrete_distribution<ll> distribution(freq.begin(), freq.end());

    int next = 0;
    #pragma omp critical
    {
      next = distribution(rng);
    }
    u = G.edges[u][next];
    cur = cur+G.filter[u];
    p.push_back(u);
  }
  return p;
}

// Baseline
path naiveRandomPath()
{
  path rPath = {(rng()%2==0)? A : B};
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

  dict_t W;
  std::set<path> r;

  std::vector<path> pathA, pathB;

  while(r.size() < Rsize)
  {
    path randomPath = naiveRandomPath();
    bool fromA = *std::begin(randomPath) == A;

    if(randomPath.size() < Q) continue;
    if(r.find(randomPath) != r.end()) continue;
    r.insert(randomPath);

    if(fromA) pathA.push_back(randomPath);
    else pathA.push_back(randomPath);
  }

  std::multiset<int> fingerprintA = fingerprint(pathA);
  std::multiset<int> fingerprintB = fingerprint(pathB);

  std::multiset<int> intersect = multiIntersect(fingerprintA, fingerprintB);
  std::multiset<int> unions = multiUnion(fingerprintA, fingerprintB);

  fj = static_cast<double>(intersect.size()) / static_cast<double>(unions.size());
  bc = 2.*static_cast<double>(intersect.size()) / static_cast<double>(fingerprintA.size() + fingerprintB.size());

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
    (      "f,first", "Drop first label in path (default: true)",             cxxopts::value(first))
    (       "s,sort", "Consider QPath as set instead of list (default: true)",cxxopts::value(sort))
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
  input = "dataset/pokec.graph";
  std::ifstream fin(input);
  std::cin.rdbuf(fin.rdbuf());

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

  N++;
  G = graph(N);

  // Check starting node
  ERROR(A >= N, "Invalid node A");
  ERROR(B >= N, "Invalid node B");

  if(A == -1) while(A == B) A = rng()%N;
  if(B == -1) while(A == B) B = rng()%N;

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
  for(int i=0; i<N; i++) G.filter[i] = bloom_filter(Z, H, rng);
  std::cerr << "end" << std::endl;


  std::cerr << "Read attributes..." << std::endl;
  input = "dataset/pokec.att";
  fin = std::ifstream(input);
  std::cin.rdbuf(fin.rdbuf());
  for(int i=1; i<N; i++)
  {
    std::string line;
    std::getline(std::cin, line);
    std::istringstream is(line);
    int id, att;
    is >> id;
    while(is >> att) G.attributes[id].insert(att);
  }
  std::cerr << "end" << std::endl;
  
  std::vector<std::multiset<int>> sampledFingerprint(N);

  // Process DP only if fcount or fsample are enabled
  std::cerr << "Start processing DP Table..." << std::endl;
  processDP();
  std::cerr << "end" << std::endl;

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

  auto sse = [](const std::vector<double> &V, double real)
  {
    if(V.size() == 0) return 0.;
    return std::accumulate(V.begin(), V.end(), 0., [real](const double sum, const double x)
    {
      return sum + (x-real)*(x-real);
    }) / V.size();
  };
  
  
  auto get_sim = [](const std::multiset<int> &fingerprintA, const std::multiset<int> &fingerprintB)
  {
    std::multiset<int> intersect = multiIntersect(fingerprintA, fingerprintB);
    std::multiset<int> unions = multiUnion(fingerprintA, fingerprintB);

    double fj = static_cast<double>(intersect.size()) / static_cast<double>(unions.size());
    double bc = 2.*static_cast<double>(intersect.size()) / static_cast<double>(fingerprintA.size() + fingerprintB.size());

    return std::array<double,2>{fj, bc};
  };
  
  // cross(a,b) = (sim(a,N(b)) + sim(b,N(a))) / 2
  auto cross_sim = [&](const int iA, const int iB)
  {
    std::multiset<int> sA = sampledFingerprint[iA];
    std::multiset<int> sB = sampledFingerprint[iB];

    std::multiset<int> nA;
    std::multiset<int> nB;

    for(int v : G.edges[iA]) nA = multiUnion(nA, sampledFingerprint[v]);
    for(int v : G.edges[iB]) nB = multiUnion(nB, sampledFingerprint[v]);

    auto cross_sA_nB = get_sim(sA, nB);
    auto cross_sB_nA = get_sim(sB, nA);
    
    double fj = (cross_sA_nB[0] + cross_sB_nA[0])/2.0;
    double bc = (cross_sA_nB[1] + cross_sB_nA[1])/2.0;
    return std::array<double,2>{fj, bc};  
  };
  
  auto attr_jaccard = [&](const int iA, const int iB)
  {
    std::multiset<int> fingerprintA = G.attributes[iA];
    std::multiset<int> fingerprintB = G.attributes[iB];
    std::multiset<int> intersect = multiIntersect(fingerprintA, fingerprintB);
    std::multiset<int> unions = multiUnion(fingerprintA, fingerprintB);
    return static_cast<double>(intersect.size()) / static_cast<double>(unions.size());
  };
  
  // 
  std::set<int> X;
  std::set<std::array<int, 2>> Y;
  std::set<int> choosenNodes;  
  size_t dim = 10;
  
  while(X.size() < dim)
  {
    int choosen = 1+(rng()%(N-1));
    if(X.find(choosen) != X.end()) continue;
    choosenNodes.insert(choosen);
    X.insert(choosen);
  }
  
  while(Y.size() < dim)
  {
    int u = 1+(rng()%(N-1));
    if(G.edges[u].size() == 0) continue;
    int v = G.edges[u][rng() % G.edges[u].size()];
    if(Y.find({u, v}) != Y.end()) continue;
    Y.insert({u, v});
    choosenNodes.insert(u);
    choosenNodes.insert(v);
  }
  
  // Compute fingerprints
  std::cerr << "Start BFS" << std::endl;
  for(int i : choosenNodes)
  {

    size_t limit = 0;
    for(auto &[k, v] : dp[Q][i]) limit += v;
    limit /= 100;
    
    std::cout << "BFS " << i << "/" << choosenNodes.size() << " " << limit << "\r";
    
    std::set<path> R;
    for(size_t t=0; t<10*limit; t++)
    {
      path toAdd = randomPathTo(i);
      if(toAdd.size() < Q) continue;
      R.insert(toAdd);
      if(R.size() == limit) break;
    }
    
    sampledFingerprint[i] = fingerprint(std::vector<path>(R.begin(), R.end()));
  }
  
  
  //X = 1000 coppie di nodi a caso
  //Y = 1000 archi a caso (che sono coppie di nodi adiacenti)

  std::vector<int> vX(X.begin(), X.end());
  std::vector<std::array<int, 2>> vY(Y.begin(), Y.end());

  std::vector<double> sim_X, cross_X, jaccard_X;
  std::vector<double> sim_Y, cross_Y, jaccard_Y;
  
  // A : similarita' media (con varianza) delle coppie in X, e delle coppie in Y 
  for(size_t i=0; i<vX.size(); i++)
  {
    for(size_t j=i+1; j<vX.size(); j++)
    {
      int u = vX[i];
      int v = vX[j];
      sim_X.push_back(get_sim(sampledFingerprint[u], sampledFingerprint[v])[0]);
      cross_X.push_back(cross_sim(u, v)[0]);
      jaccard_X.push_back(attr_jaccard(u, v));
    }
  }
  
  // B : cross similarity media (con varianza) delle coppie in X, e delle coppie in Y 
  for(auto &[u, v]: vY)
  {
    sim_Y.push_back(get_sim(sampledFingerprint[u], sampledFingerprint[v])[0]);
    cross_Y.push_back(cross_sim(u, v)[0]);
    jaccard_Y.push_back(attr_jaccard(u, v));  
  }
  
  
  double     sim_X_avg = avg(    sim_X);
  double   cross_X_avg = avg(  cross_X);
  double jaccard_X_avg = avg(jaccard_X);
  double     sim_Y_avg = avg(    sim_Y);
  double   cross_Y_avg = avg(  cross_Y);
  double jaccard_Y_avg = avg(jaccard_Y);
  
  double     sim_X_var = var(    sim_X,     sim_X_avg);
  double   cross_X_var = var(  cross_X,   cross_X_avg);
  double jaccard_X_var = var(jaccard_X, jaccard_X_avg);
  double     sim_Y_var = var(    sim_Y,     sim_Y_avg);
  double   cross_Y_var = var(  cross_Y,   cross_Y_avg);
  double jaccard_Y_var = var(jaccard_Y, jaccard_Y_avg);

  std::cout << "RESULT X" << std::endl;
  std::cout << "sim_avg, sim_var, cross_avg, cross_var, jaccard_avg, jaccard_var" << std::endl;
  std::cout <<     sim_X_avg << ", ";
  std::cout <<     sim_X_var << ", ";
  std::cout <<   cross_X_avg << ", ";
  std::cout <<   cross_X_var << ", ";
  std::cout << jaccard_X_avg << ", ";
  std::cout << jaccard_X_var << std::endl;
  
  std::cout << std::endl;

  std::cout << "RESULT Y" << std::endl;
  std::cout << "sim_avg, sim_var, cross_avg, cross_var, jaccard_avg, jaccard_var" << std::endl;
  std::cout <<     sim_Y_avg << ", ";
  std::cout <<     sim_Y_var << ", ";
  std::cout <<   cross_Y_avg << ", ";
  std::cout <<   cross_Y_var << ", ";
  std::cout << jaccard_Y_avg << ", ";
  std::cout << jaccard_Y_var << std::endl;
  
  // Flush output buffers
  std::cout.flush();
  std::cerr.flush();

  return 0;
}
