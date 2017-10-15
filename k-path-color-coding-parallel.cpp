/*
  Autore: Gaspare Ferraro
  Conta i k-path colorful in un grafo non orientato
  con la tecnica del color-coding
*/
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <omp.h>
#include <getopt.h>
#include "graph_read.hpp"

#ifndef MAXK
  #define MAXK 32
#endif

#if MAXK <= 8
  #define COLORSET uint8_t
#elif MAXK <= 16
  #define COLORSET uint16_t
#elif MAXK <= 32
  #define COLORSET uint32_t
#elif MAXK <= 64
  #define COLORSET uint64_t
#else
  #define COLORSET uint32_t
#endif
  
using namespace std;
typedef long long ll;
/*
  N = nodi del grafo
  M = archi del grafo
  k = lunghezza dei path da cercare
  kp = numero di colori da usare (>= k)
*/
unsigned int N, M, k = 0, kp = 0, thread_count = 0;
int *color;
vector<int> *G;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Ritorna il pos-esimo bit di n
bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }

// Setta il pos-esimo bit di n
COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }

// Resetta il pos-esimo bit di n
COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }

// Colora casualmente il grafo
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = rand() % kp;
}

// Hash di pair<int,int>
struct PairHash {
  size_t operator()(const std::pair<int, int> &p) const {
    return p.first ^ p.second;
  }
};

// Link
unordered_map<pair<int, COLORSET>, vector<int>, PairHash> linkGlobal;
inline void addLink(unordered_map<pair<int, COLORSET>, vector<int>, PairHash> link, int x, COLORSET C, int j) {
  auto key = make_pair(x, C);
  if (link.find(key) == link.end()) link[key] = vector<int>();
  link[key].push_back(j);
}

// Oracolo
vector<int> H(int x, COLORSET C) { return linkGlobal[make_pair(x, C)]; }

ll cont = 0;

/*
Lista e stampa i path trovati
void list_k_path(vector<int> ps, ll cs, int x)
{
  vector<int> N = H(x, cs);
  for(int v : N)
  {
    if( (ps.size() + 2) == k )
    {
      ps.push_back(v);
      for(int j : ps) printf("%d ", j);
      printf("\n");
      ps.pop_back();
    }
    else
    {
      ps.push_back(v);
      list_k_path(ps, setBit(cs, color[v]), v);
      ps.pop_back();
    }
  }
}*/

// Conta i path trovati
void list_k_path_c(vector<int> ps, COLORSET cs, int x, int kp) {
  vector<int> N = H(x, cs);
  if (kp + 2ull == k)
    cont += N.size();
  else
  {
    for (int v : N) list_k_path_c(ps, setBit(cs, color[v]), v, kp + 1);
  }
}

unordered_set<COLORSET> *DP[MAXK+1];

void processDP() {
  DP[1][N].insert(setBit(0ll, color[N]));

  for (unsigned int l = 2; l <= k; l++)
  {
      // Parallelizzare  
      #pragma omp parallel for
      for (unsigned int j = 0; j <= N; j++)
      {
        for (int x : G[j])
          for (COLORSET C : DP[l - 1][x])
            if (!getBit(C, color[j])) DP[l][j].insert(setBit(C, color[j]));  
      }  
  }
}

void backProp() {
  for (int i = k - 1; i >= 0; i--) {
    // Parallelizzare
    #pragma omp parallel for shared(linkGlobal)
    for (unsigned int x = 0; x <= N; x++) {
     // printf("%d %u\n", i, x);
      vector<ll> toDel;
      for (COLORSET C : DP[i][x]) {
        bool find = false;
        
        for (int j : G[x]) {
          if (getBit(C, color[j])) continue;

          if (DP[i + 1][j].find(setBit(C, color[j])) != DP[i + 1][j].end()) {
            find = true;
            addLink(linkGlobal, x, C, j);
          }
        }
        if (!find) toDel.push_back(C);
      }
      for (COLORSET C : toDel) DP[i][x].erase(C);
    }
  }
}

void print_usage(char *filename)
{
  printf("Usage: ./%s arguments\n",filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph (default stdin)\n");

  printf("-f, --format formatname\n");
  printf("\tFormat of input file (snap, nde, gasp)\n");

  printf("-t, --tableout filename\n");
  printf("\tOutput DP table (default stdout)\n");

  printf("-T, --tablein filename\n");
  printf("\tImport DP table (default stdin)\n");

  printf("-l, --list filename\n");
  printf("\tList k-path (default stdout)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("-h, --help\n");
  printf("\tDisplay help text and exit.\n");

  printf("-v, --verbose\n");
  printf("\tDisplay help text and exit.\n");

  
}

static int verbose_flag, help_flag;

char *input_graph = NULL;
char *table_in = NULL;
char *table_out = NULL;
char *list_path = NULL;
char *format_name = NULL;

int main(int argc, char **argv)
{

  static struct option long_options[] =
  {
    {"path",    required_argument,             0, 'k'},
    {"color",   required_argument,             0, 'K'},
    {"input",   required_argument,             0, 'g'},
    {"format",  required_argument,             0, 'f'},
    {"tableout",required_argument,             0, 't'},
    {"tablein", required_argument,             0, 'T'},
    {"list",    required_argument,             0, 'l'},
    {"parallel",required_argument,             0, 'p'},
    {"help",          no_argument, &   help_flag,  1 },
    {"verbose",       no_argument, &verbose_flag,  1 },
    {0, 0, 0, 0}
  };

  int option_index = 0;
  int c;
  while(1)
  {
    c = getopt_long (argc, argv, "k:K:g:f:t:T:l:hv", long_options, &option_index);
    
    if( c == -1 )
      break;

    switch( c )
    {
      case 'k':
        if( optarg != NULL ) k = atoi(optarg);
      break;
      case 'K':
        if( optarg != NULL ) kp = atoi(optarg);
      break;
      case 'g':
        if( optarg != NULL ) input_graph = optarg;
      break;
      case 'f':
        if( optarg != NULL ) format_name = optarg;
      break;
      case 't':
        if( optarg != NULL ) table_in = optarg;
      break;
      case 'T':
        if( optarg != NULL ) table_out = optarg;
      break;
      case 'l':
        if( optarg != NULL ) list_path = optarg;
      break;
      case 'p':
        if( optarg != NULL ) thread_count = atoi(optarg);
      break;
    }
  }


  char *input_graph = NULL;
  char *table_in = NULL;
  char *table_out = NULL;
  char *list_path = NULL;
  char *format_name = NULL;

  if( k == 0 )
    printf("Invalid or missing path length value.\n");

  if( help_flag || k == 0 )
  {
    print_usage(argv[0]);
    return 0;
  }

  if( kp == 0 ) kp = k;

  if( verbose_flag )
  {
    printf("Options:\n");
    printf("k = %d\n", k);
    printf("kp = %d\n", kp);
    printf("thread = %d\n", thread_count);

    printf("input_graph = %s\n", input_graph != NULL ? input_graph : "stdin");
    printf("format_name = %s\n", format_name != NULL ? format_name : "stdin");
    printf("table_in    = %s\n", table_in    != NULL ? table_in    : "stdin");
    printf("table_out   = %s\n", table_out   != NULL ? table_out   : "stdin");
    printf("list_path   = %s\n", list_path   != NULL ? list_path   : "stdin");
  }

  N = nextInt();
  M = nextInt();

  k = atol(argv[1]);
  kp = atol(argv[2]);

  color = new int[N + 1];
  G = new vector<int>[N + 1];

  assert(k <= MAXK);

  for(unsigned int i=0; i <= k+1; i++)
    DP[i] = new unordered_set<COLORSET>[N+1];

  for (unsigned int i = 0; i < M; i++) {
    int a = nextInt();
    int b = nextInt();
    G[a].push_back(b);
    G[b].push_back(a);
  }

  // Creo un nodo fittizio N che collegato a tutti nodi di G
  for (unsigned int i = 0; i < N; i++) {
    G[N].push_back(i);
    G[i].push_back(N);
  }

  // Coloro il grafo casualmente
  randomColor();
  color[N] = kp;
  k++;

  // Riempie la tabella di programmazione dinamica
  processDP();

  // Backward-propagation della tabella di programmazione dinamica
  backProp();

  // Conto i k-path colorful
  // list_k_path_c(vector<int>(), setBit(0ll, color[N]), N, 0);

  for(unsigned int i= 0 ; i <= k ; i++ )
    for(unsigned int j = 0 ; j <= N; j++)
      cont += DP[i][j].size();

  printf("%llu\n", cont);
  return 0;
}


