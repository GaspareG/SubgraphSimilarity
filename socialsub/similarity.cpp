#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

// CONSTANT DEFINE
#define MAXN 100000
#define MAXQ 15
#define TOPK 10
#define FACEBOOK_USER 4039
#define LINKEDIN_USER 29040
#define COLORSET uint32_t

// Flags for experiments
bool classmateFlag = false; // Test classmate x facebook
bool familyFlag = false;    // Test family x facebook
bool schoolFlag = true;     // Test school x linkedin
bool workFlag = false;      // Test work x linkedin

bool exact = false;         // true -> exact value || false -> color-coding + sampling
unsigned int q = 6;         // length of paths (n° of nodes) (FOR EXACT ONLY 3 OR 4)
unsigned int ncolor = 6;    // number of colors used for color-codign
unsigned int w = 1000;      // size of the samples
unsigned int seed = 42;     // random seed

// Subgraph categories
string classmateStrings[] = {"user", "education;concentration;id",
                             "education;school;id", "education;degree;id"};
string familyStrings[] = {"user", "hometown;id", "last_name", "location;id",
                             "last_name"};
string schoolStrings[] = {"user", "college"};
string workStrings[] = {"user", "employer"};

// File names
char filename_facebook_graph_node[] = "facebook/graph.node";
char filename_facebook_graph_edge[] = "facebook/graph.edge";

char filename_linkedin_graph_node[] = "linkedin/graph.node";
char filename_linkedin_graph_edge[] = "linkedin/graph.edge";

char filename_classmate_query_in[] = "classmate-query-node";
char filename_classmate_query_out[] = "classmate-query-node.out";

char filename_family_query_in[] = "family-query-node";
char filename_family_query_out[] = "family-query-node.out";

char filename_school_query_in[] = "school-query-node";
char filename_school_query_out[] = "school-query-node.out";

char filename_work_query_in[] = "work-query-node";
char filename_work_query_out[] = "work-query-node.out";

// exact path definition
typedef long long ll;
typedef tuple<int, int, int> label3;
typedef tuple<int, int, int, int> label4;
typedef label4 labelE;

// File pointers
FILE *node, *edge;
FILE *classmate, *classmate_out;
FILE *family, *family_out;
FILE *work, *work_out;
FILE *school, *school_out;

// Graph definition
unsigned int N, E;    // number of nodes, edges
int label[MAXN];          // id node -> id label
int color[MAXN];      // id node -> color
bool in[MAXN];        // id node -> picked in induced subgraph
vector<int> G[MAXN];  // id node -> { id node neighbour }
map<int, pair<string, int>> userCatId;  // id node -> <category, id cat>

// COLORSET functions
bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }
COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }
COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }
COLORSET getCompl(COLORSET n) { return ((1 << q) - 1) & (~n); }
inline void randomColor() {
  mt19937_64 eng = mt19937_64(seed*MAXN);
  for (unsigned int i = 0; i < MAXN; i++) color[i] = eng() % ncolor;
}

// Path label
vector<int> L(vector<int> P) {
  vector<int> l;
  for (size_t i = 0; i < P.size(); i++) l.push_back(label[P[i]]);
  return l;
}

// Color-Coding Dynamic Programming preprocessing
map<COLORSET, ll> M[MAXQ + 1][MAXN];

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

// Definition of Y(x) for facebook' users
set<int> Yfacebook(int x, string cat) {
  set<int> Guser;
  set<int> myCat;

  for(int y : G[x])
    if( userCatId[y].first == cat )
      myCat.insert(y);

  for(int y : G[x])
    if( userCatId[y].first == "user" )
      for(int z : G[y])
        if( myCat.find(z) != myCat.end() )
        {
          Guser.insert(y);
          break;
        }

  return Guser;
}

// Definition of Y(x) for linkedin' users
set<int> Ylinkedin(int x)
{
  set<int> ys;
  for (int y : G[x])
    if (userCatId[x].first == "user") ys.insert(y);
  return ys;
}

// Enumerate all the label sequences from a nodes
vector<labelE> labelsFrom(int x) {
  vector<labelE> out;

  assert( q >= 3 && q <= 4 );

  if( q == 3 )
  {
    for (int y : G[x]) {
      for (int z : G[y]) {
        if (z == x) continue;
        out.push_back(make_tuple(label[x], label[y], label[z], 0));
      }
    }
  }
  else
  {
    for (int y : G[x]) {
      for (int z : G[y]) {
        if (z == x) continue;
        for(int w : G[z])
        {
          if(w == x) continue;
          if(w == y) continue;
          out.push_back(make_tuple(label[x], label[y], label[z], label[w]));
        }
      }
    }
  }

  return out;
}
// bray-curtis weighted
double BCW(set<vector<int>> W, map<vector<int>, ll> freqA, map<vector<int>, ll> freqB) {
  ll num = 0ll;
  ll den = 0ll;
  for (auto x : W) {
    ll fax = freqA[x];
    ll fbx = freqB[x];
    num += 2 * min(fax, fbx);
    den += fax + fbx;
  }
  return (double)num / (double)den;
}

vector<int> randomPathTo(int u) {
  vector<int> P;
  P.push_back(u);
  COLORSET cs = getCompl(setBit(0ll, color[u]));
  mt19937_64 eng = mt19937_64(seed*u);
  for (int qi = q - 1; qi > 0; qi--) {
    vector<ll> freq;
    for (int v : G[u])
    {
      int add=0;
      if( M[qi][v].find(cs) != M[qi][v].end() )
        add = M[qi][v][cs];
      freq.push_back(add);
    }
    discrete_distribution<int> distribution(freq.begin(), freq.end());
    u = G[u][distribution(eng)];
    P.push_back(u);
    cs = clearBit(cs, color[u]);
  }
  reverse(P.begin(), P.end());
  return P;
}

map<pair<int, vector<int>>, ll> randomColorfulSamplePlus(vector<int> X, int r) {
  map<pair<int, vector<int>>, ll> dict;
  set<vector<int>> R;
  vector<ll> freqX;
  for (int x : X) freqX.push_back(M[q][x][getCompl(0ll)]);
  discrete_distribution<int> distribution(freqX.begin(), freqX.end());
  ll generated = 0;
  mt19937_64 eng = mt19937_64(seed*X[0]);

  while( R.size() < (size_t)r && generated < r*4 ) // avoid infinite loop
  {
    int u = X[distribution(eng)];
    vector<int> P = randomPathTo(u);
    R.insert(P);
    generated++;
  }
  for (auto r : R) {
    reverse(r.begin(), r.end());
    dict[make_pair(*r.begin(), L(r))]++;
  }
  return dict;
}

// bray-curtis via fsample approach
double fsampleBrayCurtis(int x, int y)
{
  map<vector<int>, ll> freqA, freqB;
  set<vector<int>> W;
  vector<int> X;

  X.push_back(x);
  X.push_back(y);

  for (auto dw : randomColorfulSamplePlus(X, w)) {
    int u = dw.first.first;
    W.insert(dw.first.second);
    if (u == x) freqA[dw.first.second] += dw.second;
    if (u == y) freqB[dw.first.second] += dw.second;
  }

  return BCW(W, freqA, freqB);
}

// Exact bray-curtis
double exactBrayCurtis(int x, int y)
{
  vector<labelE> Lx = labelsFrom(x);
  vector<labelE> Ly = labelsFrom(y);
  vector<labelE> Lxy(Lx.size() + Ly.size());
  auto it =
      set_intersection(Lx.begin(), Lx.end(), Ly.begin(), Ly.end(), Lxy.begin());
  Lxy.resize(it - Lxy.begin());
  double out = (double)((double)2. * Lxy.size()) / (Lx.size() + Ly.size());
  return out;
}

// bray-curtis computation
double brayCurtis(int x, int y) {
  return exact ? exactBrayCurtis(x,y) : fsampleBrayCurtis(x,y);
}

// Exact computation of top-10 similar nodes in a facebook node
int query_count=0;
vector<int> queryFacebook(int x, string cat) {
  priority_queue<pair<double, int>> Q;
  printf("%3d) X = [%d] \n", query_count++, x);
  auto yset = Yfacebook(x, cat);
  vector<int> ys(yset.begin(), yset.end());
  sort(ys.begin(), ys.end());
  #pragma omp parallel for schedule(guided)
  for (unsigned int i=0; i<ys.size(); i++)
  {
    int y = ys[i];
    // printf("\tBC(%d, %d) [%d]\n", x, y, omp_get_thread_num());
    double bcvalue = brayCurtis(x, y);
    #pragma omp critical
    {
      Q.push(make_pair(bcvalue, y));
    }
  }
  printf("\n");
  vector<int> out;
  while (!Q.empty() && out.size() < 10) {
    out.push_back(Q.top().second);
    Q.pop();
  }
  return out;
}

// Exact computation of top-10 similar nodes in a linkedin node
vector<int> queryLinkedin(int x) {
  priority_queue<pair<double, int>> Q;
  printf("%3d) X = [%d] \n", query_count++, x);
  auto yset = Ylinkedin(x);
  vector<int> ys(yset.begin(), yset.end());
  sort(ys.begin(), ys.end());
  #pragma omp parallel for schedule(guided)
  for (unsigned int i=0; i<ys.size(); i++)
  {
    int y = ys[i];
    // printf("\tBC(%d, %d) [%d]\n", x, y, omp_get_thread_num());
    double bcvalue = brayCurtis(x, y);
    #pragma omp critical
    {
      Q.push(make_pair(bcvalue, y));
    }
  }
  vector<int> out;
  while (!Q.empty() && out.size() < 10) {
    out.push_back(Q.top().second);
    Q.pop();
  }
  return out;
}

vector<int> classmateQuery(int x) { return queryFacebook(x, "education;school;id"); }
vector<int> familyQuery(int x) { return queryFacebook(x, "last_name"); }
vector<int> schoolQuery(int x) { return queryLinkedin(x); }
vector<int> workQuery(int x) { return queryLinkedin(x); }

int main() {

  // Opening files
  if( classmateFlag || familyFlag )
  {
    node = fopen(filename_facebook_graph_node, "r");
    edge = fopen(filename_facebook_graph_edge, "r");
  }
  else
  {
    node = fopen(filename_linkedin_graph_node, "r");
    edge = fopen(filename_linkedin_graph_edge, "r");
  }

  if (classmateFlag) {
    classmate = fopen(filename_classmate_query_in, "r");
    classmate_out = fopen(filename_classmate_query_out, "w");
  }

  if (familyFlag) {
    family = fopen(filename_family_query_in, "r");
    family_out = fopen(filename_family_query_out, "w");
  }

  if (schoolFlag) {
    school = fopen(filename_school_query_in, "r");
    school_out = fopen(filename_school_query_out, "w");
  }

  if (workFlag) {
    work = fopen(filename_work_query_in, "r");
    work_out = fopen(filename_work_query_out, "w");
  }

  assert(NULL != node);
  assert(NULL != edge);

  // Random generator
  srand(seed);

  // induces subgraph categories
  set<string> classmateCat(classmateStrings, classmateStrings + 4);
  set<string> familyCat(familyStrings, familyStrings + 5);
  set<string> schoolCat(schoolStrings, schoolStrings + 2);
  set<string> workCat(workStrings, workStrings + 2);

  // Users number for facebook/linkedin
  int MAXU = (classmateFlag || familyFlag  ) ? FACEBOOK_USER : LINKEDIN_USER;

  // Reading nodes
  N=0;
  while (!feof(node)) {
    int uid, labelN = 0;  // label 0 for users
    char line[1024];
    char cat[512];

    if ( NULL == fgets(line, 1023, node) )
      break;

    assert(1 == sscanf(line, "%d\t", &uid));

    if (uid < MAXU) {
      strcpy(cat, "user");
    } else {
      assert(1 == sscanf(line, "%*d\t%s", cat));
      labelN = uid;
    }

    // Filter node by category
    in[uid] = false;
    if (classmateFlag && classmateCat.find(string(cat)) == classmateCat.end())
      continue;
    else if (familyFlag && familyCat.find(string(cat)) == familyCat.end())
      continue;
    else if (schoolFlag && schoolCat.find(string(cat)) == schoolCat.end())
      continue;
    else if (workFlag && workCat.find(string(cat)) == workCat.end())
      continue;
    in[uid] = true;
    N++;

    // catUserId[string(cat)].push_back(make_pair(uid, 0));
    userCatId[uid] = make_pair(string(cat), 0);
    label[uid] = labelN;
  }

  // Reading edges
  E=0;
  while (!feof(edge)) {
    int x, y;
    assert(2 == fscanf(edge, "%d %d\n", &x, &y));
    if (!in[x] || !in[y]) continue;  // Only edges in induced subgraph
    G[x].push_back(y);
    E++;
  }

  // Preprocessing for approximated queries
  if( !exact )
  {
    // Create and process DP Table
//    for (unsigned int i = 0; i <= q + 1; i++) M[i] = new map<COLORSET, ll>[MAXN + 1];
    printf("Random coloring...\n");
    randomColor();
    printf("Process DP table...\n");
    processDP();
    printf("Finish Preprocessing...\n");
  }

  // Query processing
  if (classmateFlag) {
    while (!feof(classmate)) {
      int x;
      assert(1 == fscanf(classmate, "%d\n", &x));
      auto bcClassmate = classmateQuery(x);
      fprintf(classmate_out, "%d", x);
      for (int y : bcClassmate) fprintf(classmate_out, "\t%d", y);
      fprintf(classmate_out, "\n");
      fflush(classmate_out);
    }
  } else if (familyFlag) {
    while (!feof(family)) {
      int x;
      assert(1 == fscanf(family, "%d\n", &x));
      auto bcFamily = familyQuery(x);
      fprintf(family_out, "%d", x);
      for (int y : bcFamily) fprintf(family_out, "\t%d", y);
      fprintf(family_out, "\n");
      fflush(family_out);
    }
  } else if (schoolFlag) {
    while (!feof(school)) {
      int x;
      assert(1 == fscanf(school, "%d\n", &x));
      auto bcSchool = schoolQuery(x);
      fprintf(school_out, "%d", x);
      for (int y : bcSchool) fprintf(school_out, "\t%d", y);
      fprintf(school_out, "\n");
      fflush(school_out);
    }
  } else if (workFlag) {
    while (!feof(work)) {
      int x;
      assert(1 == fscanf(work, "%d\n", &x));
      auto bcWork = workQuery(x);
      fprintf(work_out, "%d", x);
      for (int y : bcWork) fprintf(work_out, "\t%d", y);
      fprintf(work_out, "\n");
      fflush(work_out);
    }
  }

  // Closing fiels
  fclose(node);
  fclose(edge);
  if (classmateFlag) fclose(classmate);
  if (classmateFlag) fclose(classmate_out);
  if (familyFlag) fclose(family);
  if (familyFlag) fclose(family_out);
  if (workFlag) fclose(work);
  if (workFlag) fclose(work_out);
  if (schoolFlag) fclose(school);
  if (schoolFlag) fclose(school_out);
  return 0;
}
