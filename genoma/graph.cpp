#include <bits/stdc++.h>
#include <omp.h>
#include <parallel/algorithm>

using namespace std;

#define COLORSET uint32_t
#define MAXQ 31
#define DEBUG(x) printf(x)

/*
1. parole originali
2. rotazioni cicle (dove la parte prima del $ e almeno k --> 150-k)
3. ordino
4. calcolo l'lcp array
5. per ogni inizio vado avanti finche lcp(i, j) < threshold
   i -> parola originale, j -> qualsiasi
*/

typedef long long int ll;
typedef struct {
  pair<int,int> pos;
  string seq;
  string label;
} read;

unsigned int N=0; // number of reads
unsigned int k=0; // Size of overlap
unsigned int q=0; // number of color/length of paths

vector<read> reads;
vector< pair<int, int> > permutations;
vector<int> lcp;
set<int> *Gr;
vector<int> *G;
map<COLORSET, ll> *M[MAXQ + 1];
COLORSET *color;

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

/*************************************/
// Base functions
void readInput()
{
  ios::sync_with_stdio(false);
  cin >> N;
  for(int i=0; i<N; i++)
  {
    read tmp;
    cin >> tmp.pos.first >> tmp.pos.second;
    cin >> tmp.seq;
    tmp.seq += "$";;
    reads.push_back(tmp);
  }
}

void buildGraph()
{

  Gr = new set<int>[N];
  G = new vector<int>[N];

  DEBUG("Building permutation...\n");
  for(int i=0; i<N; i++)
  {
    auto oSize = reads[i].seq.size();
    if( oSize-1 < k ) continue;
    for(unsigned j=0; j<oSize-k-1; j++)
      permutations.push_back( make_pair(j, i) );
  }

  auto comp = [&](pair<int,int> a, pair<int,int> b){
    return strcmp( &(reads[a.second].seq.c_str()[a.first]), &(reads[b.second].seq.c_str()[b.first])) < 0;
  };

  DEBUG("Sort permutation...\n");
  __gnu_parallel::sort(permutations.begin(), permutations.end(), comp);

  lcp.reserve(permutations.size());
  lcp.push_back(0);

  DEBUG("Calculated LCP Array...\n");
  #pragma omp parallel for schedule(guided)
  for(int i=1; i<permutations.size(); i++)
  {
    auto a = &(reads[permutations[i].second].seq.c_str()[permutations[i].first]);
    auto b = &(reads[permutations[i-1].second].seq.c_str()[permutations[i-1].first]);
    auto sa = strlen(a);
    auto sb = strlen(b);
    int pref = 0;
    while( a[pref] != '$' && b[pref] != '$' && a[pref] == b[pref] ) pref++;
    lcp[i] = pref;
  }

  DEBUG("Building graph edges...\n");
  #pragma omp parallel for schedule(guided)
  for(int i=0; i<permutations.size()-1; ++i)
  {
    if( permutations[i].first != 0 ) continue;
    int j=i+1;
    int minPref = lcp[j];
    while( minPref >= k && j < permutations.size() )
    {
      if( permutations[i].second != permutations[j].second )
        Gr[permutations[i].second].insert(permutations[j].second);
      minPref = min(lcp[++j], minPref);
    }
  }

  DEBUG("Reverting edges...\n");
  for(int i=0; i<N; i++)
    for(int j : Gr[i])
      G[j].push_back(i);

}

void randomColor() {
  color = new COLORSET[N];
  for (unsigned int i = 0; i < N; i++) color[i] = eng() % q;
}

void processDP()
{
  #pragma omp parallel for schedule(guided)
  for(int i=0; i<=q; i++) M[i] = new map<COLORSET, ll>[N];

  #pragma omp parallel for schedule(guided)
  for (unsigned int u = 0; u < N; u++) M[1][u][setBit(0, color[u])] = 1ll;

  for (unsigned int i = 2; i <= q; i++) {
    #pragma omp parallel for schedule(guided)
    for (unsigned int u = 0; u < N; u++) {
      for (int v : G[u]) {
        for (auto d : M[i - 1][v]) {
          COLORSET s = d.first;
          if (getBit(s, color[u])) continue;
          M[i][u][setBit(s, color[u])] += d.second;
        }
      }
    }
  }
}

/*************************************/
// Sampling functions
string L(vector<int> P) {
  string l = "";
  for (size_t i = 0; i < P.size(); i++) l += reads[P[i]].label;
  return l;
}

vector<int> randomPathTo(int u) {
  list<int> P;
  P.push_front(u);
  COLORSET D = getCompl(setBit(0l, color[u]));
  for (int i = q - 1; i > 0; i--) {
    vector<ll> freq;
    for (int v : G[u]) freq.push_back(M[i][v][D]);
    discrete_distribution<int> distribution(freq.begin(), freq.end());
    u = G[u][distribution(eng)];
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
  for (int x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  while( R.size() < (size_t)r)
  {
    int rem = r - R.size();
    for(int i=0; i<rem; i++)
    {
      int u = X[distribution(eng)];
      vector<int> P = randomPathTo(u);
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
double bcw(set<string> W, map<string, ll> freqA, map<string, ll> freqB) {
  double ret = 0.;
  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    ret += (double) min(fax, fbx) / ( fax + fbx );
  }
  return ((double) 2 / W.size()) * ret;
}

double fjw(set<string> W, map<string, ll> freqA, map<string, ll> freqB, long long R) {
  ll num = 0ll;
  for (string x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    num += min(fax, fbx);
  }
  return (double)num / (double) R;
}

int main(int argc, char **argv)
{

  if( argc < 3 )
  {
    printf("Usage: ./%s K Q\n", argv[0]);
    return 1;
  }

  k = atoi(argv[1]);
  q = atoi(argv[2]);

  DEBUG("Read input...\n");
  readInput();

  DEBUG("Building graph...\n");
  buildGraph();

  DEBUG("Random coloring graph...\n");
  randomColor();

  DEBUG("Processing DP Table...\n");
  processDP();

  return 0;
}
