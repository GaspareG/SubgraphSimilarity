#include <bits/stdc++.h>
#include <omp.h>
#include <parallel/algorithm>

using namespace std;

#define COLORSET uint32_t
#define MAXQ 31
#define DEBUG(x) fprintf(stderr, x)
#define DEBUG1(x, y) fprintf(stderr, (x), (y))
#define DEBUG2(x, y, z) fprintf(stderr, (x), (y), (z))
#define MAXN 3000000

typedef uint_fast32_t fuint;
typedef uint_fast64_t ull;

typedef struct {
  pair<fuint,fuint> pos;
  string seq;
  string label;
} read;

fuint N=0; // number of reads
fuint k=0; // Size of overlap
fuint q=0; // number of color/length of paths
fuint e=0; // number of experiments

read reads[MAXN];
// vector<read> reads;
vector< pair<fuint, fuint> > permutations;
vector<fuint> lcp;
unordered_set<fuint> Gr[MAXN];
vector<fuint> G[MAXN];
vector<fuint> sortedPos;
map<COLORSET, ull> M[MAXQ + 1][MAXN];
COLORSET *color;

// Random generator
mt19937_64 eng;
uniform_int_distribution<ull> distr;

// Get pos-th bit in n
bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }

// Set pos-th bit in n
COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }

// Reset pos-th bit in n
COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }

// Complementary set of a COLORSET
COLORSET getCompl(COLORSET n) { return ((1 << q) - 1) & (~n); }

/*************************************/
// Base functions
void readInput()
{
  ios::sync_with_stdio(false);
  cin >> N;
  // reads.resize(N);
  read tmp;
  for(fuint i=0; i<N; ++i)
  {
    cin >> tmp.pos.first >> tmp.pos.second;
    cin >> tmp.seq;
    tmp.seq += "$";
    reads[i] = tmp;
  }
}

void buildGraph()
{

  // Gr = new unordered_set<fuint>[N];
  // G = new vector<fuint>[N];
  permutations.reserve(N*150);

  DEBUG("Building permutation...\n");
  for(fuint i=0; i<N; i++)
  {
    auto oSize = reads[i].seq.size();
    if( oSize-1 < k ) continue;
    for(fuint j=0; j<oSize-k-1; j++)
      permutations.push_back( make_pair(j, i) );
  }

  auto comp = [&](pair<fuint,fuint> a, pair<fuint,fuint> b){
    return strcmp( &(reads[a.second].seq.c_str()[a.first]), &(reads[b.second].seq.c_str()[b.first])) < 0;
  };

  DEBUG("Sort permutation...\n");
  // in the remote future...
  // sort(parallel::par_unseq, permutations.begin(), permutations.end(), comp);
  __gnu_parallel::sort(permutations.begin(), permutations.end(), comp);

  lcp.reserve(permutations.size());
  lcp.push_back(0);

  DEBUG("Calculated LCP Array...\n");
  #pragma omp parallel for schedule(guided)
  for(size_t i=1; i<permutations.size(); i++)
  {
    auto a = &(reads[permutations[i].second].seq.c_str()[permutations[i].first]);
    auto b = &(reads[permutations[i-1].second].seq.c_str()[permutations[i-1].first]);
    fuint pref = 0;
    while( a[pref] != '$' && b[pref] != '$' && a[pref] == b[pref] ) ++pref;
    lcp[i] = pref;
  }

  DEBUG("Building graph edges...\n");
  #pragma omp parallel for schedule(guided)
  for(size_t i=0; i<permutations.size()-1; ++i)
  {
    if( permutations[i].first != 0 ) continue;
    size_t j=i+1;
    fuint minPref = lcp[j];
    while( minPref >= k && j < permutations.size() )
    {
      if( permutations[i].second != permutations[j].second )
        Gr[permutations[i].second].insert(permutations[j].second);
      minPref = min(lcp[++j], minPref);
    }
  }

  DEBUG("Reverting edges...\n");
  for(fuint i=0; i<N; i++)
  {
    for(fuint j : Gr[i])
      G[j].push_back(i);
    Gr[i].clear();
  }
}

void randomColor() {
  color = new COLORSET[N];
  for (fuint i = 0; i < N; i++) color[i] = eng() % q;
}

void processDP()
{
  //#pragma omp parallel for schedule(guided)
  //for(fuint i=0; i<=q; i++) M[i] = new map<COLORSET, ull>[N];

  #pragma omp parallel for schedule(guided)
  for (fuint u = 0; u < N; u++) M[1][u][setBit(0, color[u])] = 1ll;

  for (fuint i = 2; i <= q; i++) {
    #pragma omp parallel for schedule(guided)
    for (fuint u = 0; u < N; u++) {
      for (int v : G[u]) {
        for (auto d : M[i - 1][v]) {
          if (getBit(d.first, color[u])) continue;
          M[i][u][setBit(d.first, color[u])] += d.second;
        }
      }
    }
  }
}

/*************************************/
// Sampling functions
string L(vector<fuint> P) {
  string l = "";
  for (size_t i = 0; i < P.size(); i++) l += reads[P[i]].label;
  return l;
}

vector<fuint> randomPathTo(fuint u) {
  list<fuint> P;
  P.push_front(u);
  COLORSET D = getCompl(setBit(0l, color[u]));
  for (fuint i = q - 1; i > 0; i--) {
    vector<ull> freq;
    for (fuint v : G[u]) freq.push_back(M[i][v][D]);
    discrete_distribution<fuint> distribution(freq.begin(), freq.end());
    u = G[u][distribution(eng)];
    P.push_front(u);
    D = clearBit(D, color[u]);
  }
  vector<fuint> ret;
  ret.clear();
  ret = vector<fuint>(begin(P), end(P));
  return ret;
}

set<string> randomColorfulSample(vector<fuint> X, fuint r) {
  set<string> W;
  set<vector<fuint>> R;
  vector<ull> freqX;
  for (fuint x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<fuint> distribution(freqX.begin(), freqX.end());
  while (R.size() < (size_t)r) {
    fuint u = X[distribution(eng)];
    vector<fuint> P = randomPathTo(u);
    if (R.find(P) == R.end()) R.insert(P);
  }
  for (auto r : R) {
    reverse(r.begin(), r.end());
    W.insert(L(r));
  }
  return W;
}

map<pair<fuint, string>, ull> randomColorfulSamplePlus(vector<fuint> X, fuint r) {
  map<pair<fuint, string>, ull> W;
  set<vector<fuint>> R;
  vector<ull> freqX;
  for (fuint x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while( R.size() < (size_t)r)
  {
    fuint rem = r - R.size();
    for(fuint i=0; i<rem; i++)
    {
      fuint u = X[distribution(eng)];
      vector<fuint> P = randomPathTo(u);
      R.insert(P);
    }
  }
  for (auto r : R) {
    reverse(r.begin(), r.end());
    W[make_pair(*r.begin(), L(r))]++;
  }
  return W;
}

/*************************************/
// Similarity functions
double bcw(set<string> W, map<string, ull> freqA, map<string, ull> freqB) {
  double ret = 0.;
  for (string x : W) {
    ull fax = freqA[x];
    ull fbx = freqB[x];
    ret += (double) min(fax, fbx) / ( fax + fbx );
  }
  return ((double) 2 / W.size()) * ret;
}

double fjw(set<string> W, map<string, ull> freqA, map<string, ull> freqB, ull R) {
  ull num = 0ll;
  for (string x : W)
    num += min(freqA[x], freqB[x]);
  return (double) num / (double) R;
}

/*************************************/
int main(int argc, char **argv)
{

  if( argc < 4 )
  {
    printf("Usage: %s K Q E\n", argv[0]);
    printf("\tK = length of matching suffixes-prefixes\n");
    printf("\tQ = length of paths\n");
    printf("\tE = number of experiments\n");
    return 0;
  }

  k = atoi(argv[1]);
  q = atoi(argv[2]);
  e = atoi(argv[3]);

  DEBUG("Read input...\n");
  readInput();

  DEBUG("Building graph...\n");
  buildGraph();

  DEBUG("Random coloring graph...\n");
  randomColor();

  DEBUG("Processing DP Table...\n");
  processDP();

  DEBUG("Sorting positions...\n");
  sortedPos.resize(N);
  iota(sortedPos.begin(),sortedPos.end(), 0);

  __gnu_parallel::sort(sortedPos.begin(), sortedPos.end(), [&](fuint a, fuint b){
    return reads[a].pos < reads[b].pos;
  });

  DEBUG("Similarity tests...\n");

  for(fuint i=0; i<e; i++)
  {
    DEBUG1("Experiment #%lu\n", i);
  }

  return 0;
}
