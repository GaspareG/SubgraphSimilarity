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
  input = "mep.graph";
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
  input = "mep.att";
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

  // Read infos
  input = "mep.graphnames";
  fin = std::ifstream(input);
  std::cin.rdbuf(fin.rdbuf());

  std::vector<int> voter_id(N), party_id(N); // id -> voter / id -> twitter
  std::vector<ll> twitter_id(N);
  std::vector<std::string> eu_group(N), name(N), party(N); // id -> eu_group, id -> name
  
  std::map<std::string, std::vector<int>> eu_member;// eu_group -> {id}
  std::map<std::string, std::vector<int>> party_member; // party_name -> {id}
  std::map<int, std::string> party_to_name; // party_id -> name
  std::map<std::string, std::string> party_to_eu; // party_name -> eu_name
  
  for(int i=1; i<N; i++)
  {
    // id, 
    // voter id, 
    // twitter id, 
    // EU group, 
    // national party id, 
    // national party name, 
    // name
    
    std::string str;
    std::getline(std::cin, str);
    std::stringstream buffer(str);
    std::string temp;
    std::vector<std::string> values;
    while(getline(buffer, temp, '\t'))
    {
      if(temp == "null") temp = "-1";
      values.push_back(temp);
	  }
	
    int id                 =  std::stoi(values[0]);
    int voter              =  std::stoi(values[1]);
    ll twitter             = std::stoll(values[2]);
    std::string eu_name    =            values[3];
    int party_cod          =  std::stoi(values[4]);
    std::string party_name =            values[5];
    std::string id_name    =            values[6];
    
    voter_id[id] = voter;
    twitter_id[id] = twitter;
    party_id[id] = party_cod;
    eu_group[id] = eu_name;
    name[id] = id_name;
    party[id] = party_name;
    
    eu_member[eu_name].push_back(id);
    party_member[party_name].push_back(id);
    
    party_to_name[party_cod] = party_name;
    
    party_to_eu[party_name] = eu_name;
  }
  

  // Aggiunge partiti large e avg
  int max_party = 0ll;
  int avg_party = 0ll;
  for(auto &[party_name, members] : party_member)
  {
  	max_party = std::max(max_party, static_cast<int>(members.size()));
  	avg_party += members.size();
  }
  avg_party /= party_member.size();

  for(int i=0; i<max_party; i++) party_member["random-large"].push_back(1+(rng()%(N-1)));
  for(int i=0; i<avg_party; i++) party_member["random"].push_back(1+(rng()%(N-1)));


  std::vector<std::multiset<int>> realFingerprint(N);
  std::vector<std::multiset<int>> sampledFingerprint(N);

  // Process DP only if fcount or fsample are enabled
  std::cerr << "Start processing DP Table..." << std::endl;
  processDP();
  std::cerr << "end" << std::endl;

  std::cerr << "Start BFS" << std::endl;
  #pragma omp parallel for
  for(int i=1; i<N; i++)
  {
    //std::vector<path> pathI = bfs(i); // generate all paths from I
    //realFingerprint[i] = fingerprint(pathI);
    
    // size_t limit = realFingerprint[i].size() / 10; // 10% threshold
    size_t limit = 0;
    for(auto &[k, v] : dp[Q][i]) limit += v;
    limit /= 10;
    
    std::set<path> R;
    for(size_t t=0; t<10*limit; t++)
    {
      path toAdd = randomPathTo(i);
      if(toAdd.size() < Q) continue;
      R.insert(toAdd);
      if(R.size() == limit) break;
    }
    
    // std::cerr << "I = " << i << " L = " << limit << " SAMPLED = " << R.size() << std::endl;
    
    sampledFingerprint[i] = fingerprint(std::vector<path>(R.begin(), R.end()));
  }

  auto get_sim = [](const std::multiset<int> &fingerprintA, const std::multiset<int> &fingerprintB)
  {
    std::multiset<int> intersect = multiIntersect(fingerprintA, fingerprintB);
    std::multiset<int> unions = multiUnion(fingerprintA, fingerprintB);

    double fj = static_cast<double>(intersect.size()) / static_cast<double>(unions.size());
    double bc = 2.*static_cast<double>(intersect.size()) / static_cast<double>(fingerprintA.size() + fingerprintB.size());

    return std::array<double,2>{fj, bc};
  };
  
  /*******************************************************************/
  /* EXP 1: per ogni partito, media&varianza delle sim dei membri */
  /*******************************************************************/
  /*
  std::cout << "party\tfj_avg\tfj_var\tbc_avg\tbc_var\treal_fj_avg\treal_fj_var\treal_bc_avg\treal_bc_var" << std::endl;
  for(auto &[party_name, members] : party_member)
  {
    if(members.size() <= 2) continue;
    std::vector<double> fj_vec, bc_vec;
    std::vector<double> real_fj_vec, real_bc_vec;
    for(size_t i=0; i<members.size(); i++)
    {
      for(size_t j=i+1; j<members.size(); j++)
      {
        int A = members[i];
        int B = members[j];
        auto sim_A_B = get_sim(sampledFingerprint[A], sampledFingerprint[B]);
        auto real_sim_A_B = get_sim(realFingerprint[A], realFingerprint[B]);
        fj_vec.push_back(sim_A_B[0]);
        bc_vec.push_back(sim_A_B[1]);
        real_fj_vec.push_back(real_sim_A_B[0]);
        real_bc_vec.push_back(real_sim_A_B[1]);
      }
    }
    
    double fj_avg = avg(fj_vec);
    double bc_avg = avg(bc_vec);
    double real_fj_avg = avg(real_fj_vec);
    double real_bc_avg = avg(real_bc_vec);

    double fj_var = var(fj_vec, fj_avg);
    double bc_var = var(bc_vec, bc_avg);
    double real_fj_var = var(real_fj_vec, real_fj_avg);
    double real_bc_var = var(real_bc_vec, real_bc_avg);

    std::cout << party_name << "\t";
    std::cout << fj_avg << "\t";
    std::cout << fj_var << "\t";
    std::cout << bc_avg << "\t";
    std::cout << bc_var << "\t";
    std::cout << real_fj_avg << "\t";
    std::cout << real_fj_var << "\t";
    std::cout << real_bc_avg << "\t";
    std::cout << real_bc_var << std::endl;
  }
  */
  
  /*******************************************************************/
  /* EXP 2: la similarita' media (con varianza) di un set di coppie prese a caso */
  /*******************************************************************/
  /*
  std::cout << "|A|\t|B|\tfj_avg\tfj_var\tbc_avg\tbc_var\treal_fj_avg\treal_fj_var\treal_bc_avg\treal_bc_var" << std::endl;
  
  for(int i=10; i<=100; i+=10)
  {
    for(int j=10; j<=100; j+=10)
    {
      std::vector<int> idx(N-1);
      std::iota(idx.begin(), idx.end(), 1);
      std::shuffle(idx.begin(), idx.end(), rng);
      
      // A = idx[0:i]
      // B = idx[-j:]
 
      // SAMPLED
      std::multiset<int> party1_finger = sampledFingerprint[idx[0]];
      for(size_t k=1; k<i; k++)
        party1_finger = multiUnion(party1_finger, sampledFingerprint[members1[j]]);
       
      std::multiset<int> party2_finger = sampledFingerprint[members2[0]];
      for(size_t k=1; k<j; k++)
        party2_finger = multiUnion(party2_finger, sampledFingerprint[members2[j]]);

      // REAL
   
      auto sim_A_B = get_sim(setI_finger, setJ_finger);
      auto sim_A_B_real = get_sim(real_setI_finger, real_setJ_finger);

      party_party_sim.emplace_back(-sim_A_B[0], -sim_A_B[1], -sim_A_B_real[0], -sim_A_B_real[1], party_name2);
 
    }  
  }
*/
  
  /*******************************************************************/
  /* EXP 3: C1 : ogni politico -> similarita' con tutti i partiti    */
  /*******************************************************************/
  
  std::cout << "name\tparty\tfj_avg\tbc_avg" << std::endl;
  int validi = 0;
  int validi_eu = 0;
  int tot_validi = 0;
  for(int i=1; i<N; i++)
  {
    std::vector<std::tuple<double, double, std::string>> member_party_sim;
    
    for(auto &[party_name, members] : party_member)
    {
      if(members.size() <= 2) continue;
      if(party_name == "-1") continue;
      if(party_name == "random") continue;
      if(party_name == "random-large") continue;
      double fj_party = 0.;
      double bc_party = 0.;
      
      std::multiset<int> party_fingerprint = sampledFingerprint[members[0]];
      
      for(size_t j=1; j<members.size(); j++)
        party_fingerprint = multiUnion(party_fingerprint, sampledFingerprint[members[j]]);
      
      auto sim_A_B = get_sim(sampledFingerprint[i], party_fingerprint);
      
      fj_party = sim_A_B[0];
      bc_party = sim_A_B[1];
      member_party_sim.emplace_back(-fj_party, -bc_party, party_name);
    }
    
    std::sort(member_party_sim.begin(), member_party_sim.end());
    
    for(size_t j=0; j<member_party_sim.size(); j++)
    {
      if(j == 10) break;
      std::cout << name[i] << "\t";
      std::cout << std::get<2>(member_party_sim[j]) << "\t";
      std::cout << -std::get<0>(member_party_sim[j]) << "\t";
      std::cout << -std::get<1>(member_party_sim[j]) << std::endl;
    }

    if(party_to_name[party_id[i]] == "-1") continue;    
    if(party_to_name[party_id[i]] == "random") continue;    
    if(party_to_name[party_id[i]] == "random-large") continue;    
    if(party_member[party_to_name[party_id[i]]].size() <= 2) continue;
    
    tot_validi++;  
      
    if(std::get<2>(member_party_sim[0]) == party_to_name[party_id[i]]) validi++;
    else if(std::get<2>(member_party_sim[1]) == party_to_name[party_id[i]]) validi++;
    else if(std::get<2>(member_party_sim[2]) == party_to_name[party_id[i]]) validi++;
    else if(std::get<2>(member_party_sim[3]) == party_to_name[party_id[i]]) validi++;
    else if(std::get<2>(member_party_sim[4]) == party_to_name[party_id[i]]) validi++;

    if(party_to_eu[std::get<2>(member_party_sim[0])] == party_to_eu[party_to_name[party_id[i]]]) validi_eu++;
    else if(party_to_eu[std::get<2>(member_party_sim[1])] == party_to_eu[party_to_name[party_id[i]]]) validi_eu++;
    else if(party_to_eu[std::get<2>(member_party_sim[2])] == party_to_eu[party_to_name[party_id[i]]]) validi_eu++;
    else if(party_to_eu[std::get<2>(member_party_sim[3])] == party_to_eu[party_to_name[party_id[i]]]) validi_eu++;
    else if(party_to_eu[std::get<2>(member_party_sim[4])] == party_to_eu[party_to_name[party_id[i]]]) validi_eu++;    
  }
  std::cout << "VALIDI: " << validi << "/" << tot_validi << std::endl;
  std::cout << "EU VALIDI: " << validi_eu << "/" << tot_validi << std::endl;
  
  
  /*******************************************************************/
  /* EXP 3: C2 : ogni politico -> similarita' con tutti i gruppi EU  */
  /*******************************************************************/
  /*
  std::cout << "name\teu\tfj_avg\tbc_avg" << std::endl;
  for(int i=1; i<N; i++)
  {
    std::vector<std::tuple<double, double, std::string>> member_eu_sim;
    
    for(auto &[eu_name, members] : eu_member)
    {
      if(members.size() <= 2) continue;
      double fj_eu = 0.;
      double bc_eu = 0.;
      
      std::multiset<int> eu_fingerprint = sampledFingerprint[members[0]];
      
      for(size_t j=1; j<members.size(); j++)
        eu_fingerprint = multiUnion(eu_fingerprint, sampledFingerprint[members[j]]);
      
      auto sim_A_B = get_sim(sampledFingerprint[i], eu_fingerprint);
      
      fj_eu = sim_A_B[0];
      bc_eu = sim_A_B[1];
      member_eu_sim.emplace_back(-fj_eu, -bc_eu, eu_name);
    }
    
    std::sort(member_eu_sim.begin(), member_eu_sim.end());
    
    for(size_t j=0; j<member_eu_sim.size(); j++)
    {
      if(j == 20) break;
      std::cout << name[i] << "\t";
      std::cout << std::get<2>(member_eu_sim[j]) << "\t";
      std::cout << -std::get<0>(member_eu_sim[j]) << "\t";
      std::cout << -std::get<1>(member_eu_sim[j]) << std::endl;
    }
  }
  */
  
  /*******************************************************************/
  /* EXP 4: D1 : ogni partito -> gli altri partiti (anche qui top 20) */
  /*******************************************************************/
  /*
  std::cout << "party_1\tparty_2\tfj_avg\tbc_avg\tfj_real\tbc_real" << std::endl;
  for(auto &[party_name1, members1] : party_member)
  {
    if(members1.size() <= 2) continue;

    std::vector<std::tuple<double, double, double, double, std::string>> party_party_sim;

    for(auto &[party_name2, members2] : party_member)
    {
      if(members2.size() <= 2) continue;
       
      // SAMPLED
      std::multiset<int> party1_finger = sampledFingerprint[members1[0]];
      for(size_t j=1; j<members1.size(); j++)
        party1_finger = multiUnion(party1_finger, sampledFingerprint[members1[j]]);
       
      std::multiset<int> party2_finger = sampledFingerprint[members2[0]];
      for(size_t j=1; j<members2.size(); j++)
        party2_finger = multiUnion(party2_finger, sampledFingerprint[members2[j]]);

      // REAL
      //std::multiset<int> real_party1_finger = realFingerprint[members1[0]];
      //for(size_t j=1; j<members1.size(); j++)
      //  real_party1_finger = multiUnion(real_party1_finger, realFingerprint[members1[j]]);
       
      //std::multiset<int> real_party2_finger = sampledFingerprint[members2[0]];
      //for(size_t j=1; j<members2.size(); j++)
      //  real_party2_finger = multiUnion(real_party2_finger, realFingerprint[members2[j]]);

      auto sim_A_B = get_sim(party1_finger, party2_finger);
      auto sim_A_B_real = std::array<double, 2>{0.0, 0.0}; // get_sim(real_party1_finger, real_party2_finger);

      party_party_sim.emplace_back(-sim_A_B[0], -sim_A_B[1], -sim_A_B_real[0], -sim_A_B_real[1], party_name2);
 
    }

    std::sort(party_party_sim.begin(), party_party_sim.end());
    
    for(size_t j=0; j<party_party_sim.size(); j++)
    {
      if(party_party_sim.size() >= 22 && 11 <= j && j <= party_party_sim.size() - 11) continue;
      std::cout << party_name1 << "\t";
      std::cout << std::get<4>(party_party_sim[j]) << "\t";
      std::cout << -std::get<0>(party_party_sim[j]) << "\t";
      std::cout << -std::get<1>(party_party_sim[j]) << "\t";
      std::cout << -std::get<2>(party_party_sim[j]) << "\t";
      std::cout << -std::get<3>(party_party_sim[j]) << std::endl;
    }
    
  }
  */
  
  /*******************************************************************/
  /* EXP 4: D2 : ogni partito -> gruppi eu (anche qui top 20)        */
  /*******************************************************************/
  
  /*
  int validi = 0;
  int tot_validi = 0;

  std::cout << "party_1\teu_2\tfj_avg\tbc_avg\tfj_real\tbc_real" << std::endl;
  for(auto &[party_name1, members1] : party_member)
  {
    if(members1.size() <= 2) continue;

    std::vector<std::tuple<double, double, double, double, std::string>> party_party_sim;

    for(auto &[eu_name2, members2] : eu_member)
    {
      if(members2.size() <= 2) continue;
       
      // SAMPLED
      std::multiset<int> party1_finger = sampledFingerprint[members1[0]];
      for(size_t j=1; j<members1.size(); j++)
        party1_finger = multiUnion(party1_finger, sampledFingerprint[members1[j]]);
       
      std::multiset<int> party2_finger = sampledFingerprint[members2[0]];
      for(size_t j=1; j<members2.size(); j++)
        party2_finger = multiUnion(party2_finger, sampledFingerprint[members2[j]]);

      // REAL
      //std::multiset<int> real_party1_finger = realFingerprint[members1[0]];
      //for(size_t j=1; j<members1.size(); j++)
      //  real_party1_finger = multiUnion(real_party1_finger, realFingerprint[members1[j]]);
       
      //std::multiset<int> real_party2_finger = sampledFingerprint[members2[0]];
      //for(size_t j=1; j<members2.size(); j++)
      //  real_party2_finger = multiUnion(real_party2_finger, realFingerprint[members2[j]]);

      
      auto sim_A_B = get_sim(party1_finger, party2_finger);
      auto sim_A_B_real = std::array<double,2>{0.0, 0.0}; // get_sim(real_party1_finger, real_party2_finger);

      party_party_sim.emplace_back(-sim_A_B[0], -sim_A_B[1], -sim_A_B_real[0], -sim_A_B_real[1], eu_name2);
 
    }

    std::sort(party_party_sim.begin(), party_party_sim.end());
    
    for(size_t j=0; j<party_party_sim.size(); j++)
    {
      if(j == 20) break;
      std::cout << party_name1 << "\t";
      std::cout << std::get<4>(party_party_sim[j]) << "\t";
      std::cout << -std::get<0>(party_party_sim[j]) << "\t";
      std::cout << -std::get<1>(party_party_sim[j]) << "\t";
      std::cout << -std::get<2>(party_party_sim[j]) << "\t";
      std::cout << -std::get<3>(party_party_sim[j]) << std::endl;
    }
    
    if(party_name1 == "-1") continue;    
    if(party_name1 == "random") continue;    
    if(party_name1 == "random-large") continue;    
    if(members1.size() <= 2) continue;
    
    tot_validi++;    
    
         if(std::get<4>(party_party_sim[0]) == party_to_eu[party_name1]) validi++;
    else if(std::get<4>(party_party_sim[1]) == party_to_eu[party_name1]) validi++;
    else if(std::get<4>(party_party_sim[2]) == party_to_eu[party_name1]) validi++;
    //else if(std::get<4>(party_party_sim[3]) == party_to_eu[party_name1]) validi++;
    //else if(std::get<4>(party_party_sim[4]) == party_to_eu[party_name1]) validi++;
    //else if(std::get<2>(member_party_sim[5]) == party_to_name[party_id[i]]) validi++;
    //else if(std::get<2>(member_party_sim[6]) == party_to_name[party_id[i]]) validi++;
    //else if(std::get<2>(member_party_sim[7]) == party_to_name[party_id[i]]) validi++;
    //else if(std::get<2>(member_party_sim[8]) == party_to_name[party_id[i]]) validi++;
    //else if(std::get<2>(member_party_sim[9]) == party_to_name[party_id[i]]) validi++;
    
  }
  std::cout << "VALIDI: " << validi << "/" << tot_validi << std::endl;
  */
  
  /*
  VALIDI: 34/46
  $ ./nSimGram-BloomFilter-multilabeled --bruteforce --baseline --fsample -e 10 -r 100 -t 4 -Q 4 -Z 32 -H 24 --verbose
  */
  
  // Flush output buffers
  std::cout.flush();
  std::cerr.flush();

  return 0;
}
