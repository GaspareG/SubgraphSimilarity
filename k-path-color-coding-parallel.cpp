/*
  Author: Gaspare Ferraro
  Count and find simple k colorful-path in a graph
  using the color-coding technique (parallel version)
*/
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <omp.h>
#include <getopt.h>
#include <unistd.h>

#ifdef K_8
  #define MAXK 8
  #define COLORSET uint8_t
#elif K_16
  #define MAXK 16
  #define COLORSET uint16_t
#elif K_64
  #define MAXK 64
  #define COLORSET uint64_t
#else
  #define MAXK 32
  #define COLORSET uint32_t
#endif

using namespace std;
typedef long long ll;

unsigned int N, M;
unsigned k = 0, kp = 0;
unsigned thread_count = 0;
static int verbose_flag, help_flag;


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
unordered_map< pair<int, COLORSET>, vector<int>, PairHash> links;

inline void addLink(int x, COLORSET C, int j) {
  auto key = make_pair(x, C);
  if (links.find(key) == links.end()) links[key] = vector<int>();
  links[key].push_back(j);
}

// Oracolo
vector<int> H(int x, COLORSET C) { return links[make_pair(x, C)]; }
ll cont = 0;

void list_k_path(FILE *out, vector<int> ps, COLORSET cs, int x)
{
  vector<int> oracle = H(x, cs);

  if( ps.size()+1 == k )
  {
    cont++;
    for(int j : ps) fprintf(out, "%d ", j);
    fprintf(out, "\n");
  }
  else
    for(int v : oracle)
    {
      ps.push_back(v);
      list_k_path(out, ps, setBit(cs, color[v]), v);
      ps.pop_back();
    }
}

unordered_set<COLORSET> *DP[MAXK+2];

void processDP() {
  DP[1][N].insert(setBit(0ll, color[N]));

  for (unsigned int l = 2; l <= k; l++)
  {
      if( verbose_flag ) printf("K = %d\n", l);
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
  for (int i = k - 1; i > 0; i--) {
    if( verbose_flag ) printf("K = %d\n", i);
    #pragma omp parallel for
    for (unsigned int x = 0; x <= N; x++) {
      vector<ll> toDel;
      for (COLORSET C : DP[i][x]) {
        bool find = false;

        for (int j : G[x]) {
          if (getBit(C, color[j])) continue;

          if (DP[i + 1][j].find(setBit(C, color[j])) != DP[i + 1][j].end()) {
            find = true;
            #pragma omp critical
            {
              addLink(x, C, j);
            }
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
//  printf("Usage: ./%s -k length -K number -g filename -f format -t filename -T filename -p threadcount -h -v\n",filename);
  printf("Usage: ./%s -k length -K number -g filename -f format -p threadcount --help --verbose\n",filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph (default stdin)\n");

  printf("-f, --format format\n");
  printf("\tFormat of input file (snap, nde, nme)\n");
/*
  printf("-t, --tableout filename\n");
  printf("\tOutput DP table (default stdout)\n");

  printf("-T, --tablein filename\n");
  printf("\tImport DP table (default stdin)\n");
*/
  printf("-l, --list filename\n");
  printf("\tList k-path (default stdout)\n");

  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\tDisplay help text and exit.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;
/*
bool table_in_flag = false;
char *table_in = NULL;

bool table_out_flag = false;
char *table_out = NULL;
*/
bool list_path_flag = false;
char *list_path = NULL;

bool format_name_flag = false;
char *format_name = NULL;

int main(int argc, char **argv)
{

  static struct option long_options[] =
  {
    {"path",    required_argument,             0, 'k'},
    {"color",   required_argument,             0, 'K'},
    {"input",   required_argument,             0, 'g'},
    {"format",  required_argument,             0, 'f'},
    // {"tableout",required_argument,             0, 't'},
    // {"tablein", required_argument,             0, 'T'},
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
    // c = getopt_long (argc, argv, "k:K:g:f:t:T:l:p:", long_options, &option_index);
    c = getopt_long (argc, argv, "k:K:g:f:l:p:", long_options, &option_index);

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
        input_graph_flag = true;
        if( optarg != NULL ) input_graph = optarg;
      break;
      case 'f':
        format_name_flag = true;
        if( optarg != NULL ) format_name = optarg;
      break;
/*      case 't':
        table_in_flag = true;
        if( optarg != NULL ) table_in = optarg;
      break;
      case 'T':
        table_out_flag = true;
        if( optarg != NULL ) table_out = optarg;
      break;*/
      case 'l':
        list_path_flag = true;
        if( optarg != NULL ) list_path = optarg;
      break;
      case 'p':
        if( optarg != NULL ) thread_count = atoi(optarg);
      break;
    }
  }

  if( help_flag || argc == 1 )
  {
    print_usage(argv[0]);
    return 0;
  }

  if( k == 0 )
  {
    printf("Invalid or missing path length value.\n");
    return 1;
  }

  if( k > MAXK || kp > MAXK )
  {
    printf("k or kp to high! (max value: %d\n", MAXK);
    return 1;
  }

  if( kp == 0 ) kp = k;

  if( thread_count > 0 && (int) thread_count < omp_get_max_threads() )
  {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
  }

  if( verbose_flag )
  {
    printf("Options:\n");
    printf("k = %d\n", k);
    printf("kp = %d\n", kp);
    printf("thread = %d\n", thread_count);
    printf("input_graph = %s\n", input_graph != NULL ? input_graph : "stdin");
    printf("format_name = %s\n", format_name != NULL ? format_name : "stdin");
    // printf("table_in    = %s\n", table_in    != NULL ? table_in    : "stdin");
    // printf("table_out   = %s\n", table_out   != NULL ? table_out   : "stdin");
    printf("list_path   = %s\n", list_path   != NULL ? list_path   : "stdin");
  }
  if( verbose_flag ) printf("Reading graph...\n");

  if( input_graph_flag )
  {
    if( input_graph == NULL )
    {
      printf("Input file name missing!\n");
      return 1;
    }
    if( !format_name_flag || format_name == NULL )
    {
      printf("Input format missing!\n");
      return 1;
    }
    if( strcmp(format_name, "snap") == 0)
    {
      FILE *input_fd = fopen(input_graph, "r");
      if( input_fd == NULL )
      {
        perror("Error opening input file");
        return 1;
      }
      set< pair<int,int> > edge;
      char *buffer;
      size_t n;
      int ab[2];
      N = 0;
      do {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if( line_length == 0 ) continue;
        if( buffer[0] == '#' ) continue;
        sscanf(buffer, "%d %d", ab, ab+1);
        edge.insert(make_pair(ab[0], ab[1]));
        N = (unsigned) ab[0] > N ? ab[0] : N;
        N = (unsigned) ab[1] > N ? ab[1] : N;
      }
      while( !feof(input_fd) );

      M = edge.size();
      color = new int[N + 1];
      G = new vector<int>[N + 1];
      for(auto e : edge)
      {
        G[ e.first ].push_back( e.second );
        G[ e.second ].push_back( e.first );
      }
    }
    else if( strcmp(format_name, "nde") == 0 )
    {
      FILE *input_fd = fopen(input_graph, "r");
      if( input_fd == NULL )
      {
        perror("Error opening input file");
        return 1;
      }
      char *buffer;
      size_t n;

      do {      
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if( line_length == 0 ) continue;
        if( buffer[0] == '#' ) continue;
        break;
      }
      while( !feof(input_fd) );
      sscanf(buffer, "%u",&N);
      int Nrim = N;

      M = 0;
      int ab[2];
      while( Nrim > 0 && !feof(input_fd) ) {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if( line_length == 0 ) continue;
        if( buffer[0] == '#' ) continue;
        sscanf(buffer, "%d %d", ab, ab+1);
        Nrim--; 
        M += ab[1];
        N = (unsigned) ab[0] > N ? (unsigned) ab[0] : N;
      }
      
      color = new int[N + 1];
      G = new vector<int>[N + 1];

      int Mrim = M;
      while( Mrim > 0 && !feof(input_fd) ) {
        free(buffer);
        buffer = NULL;
        n = 0;
        int line_length = ::getline(&buffer, &n, input_fd);
        if( line_length == 0 ) continue;
        if( buffer[0] == '#' ) continue;
        sscanf(buffer, "%d %d", ab, ab+1);
        G[ ab[0] ].push_back( ab[1] );
        G[ ab[1] ].push_back( ab[0] );
      }

    }
    else if( strcmp(format_name, "nme") == 0)
    {
      int input_fd = open(input_graph, O_RDONLY, 0);
      if( input_fd == -1 )
      {
        perror("Error opening input file");
        return 1;
      }
      read(input_fd, &N, sizeof(int));
      read(input_fd, &M, sizeof(int));

      color = new int[N + 1];
      G = new vector<int>[N + 1];
      int ab[2];
      for(unsigned int i=0; i<M; i++)
      {
        read(input_fd, ab, 2*sizeof(int));
        G[ab[0]].push_back(ab[1]);
        G[ab[1]].push_back(ab[0]);
      }
    }
    else
    {
      printf("Wrong input format (only 'snap', 'nde' or 'nme'\n");
      return 1;
    }
  }
  else
  {
    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();

    color = new int[N + 1];
    G = new vector<int>[N + 1];

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  if( verbose_flag ) printf("N = %d | M = %d\n", N, M);

  for (unsigned int i = 0; i < N; i++) {
    G[N].push_back(i);
    G[i].push_back(N);
  }

  // Create DP Table
  for(unsigned int i=0; i <= k+1; i++)
    DP[i] = new unordered_set<COLORSET>[N+1];

  // Random color graph
  if( verbose_flag ) printf("Random coloring graph...\n");
  randomColor();
  color[N] = kp;
  k++;

  // Fill dynamic programming table
  if( verbose_flag ) printf("Processing DP table...\n");
  processDP();
  if( verbose_flag ) printf("End processing DP table...\n");

  // Backward-propagation of DP table
  if( verbose_flag ) printf("Backward propagation...\n");
  backProp();
  if( verbose_flag ) printf("End backward propagation...\n");

  // Count ad list k-colorful path
  if( list_path_flag )
  {
    FILE *list_fd = stdout;
    if( list_path != NULL )
    {
      list_fd = fopen(list_path, "w");
      if( list_fd == NULL )
      {
        perror("Error opening list file");
        return 1;
      }
    }
    if( verbose_flag ) printf("Listing k-path...\n");
    list_k_path(list_fd, vector<int>(), setBit(0, color[N]), N);
    if( verbose_flag ) printf("%llu k-path found!\n", cont);
    fclose(list_fd);
  }

  cont = 0;
  for(unsigned int i= 0 ; i <= k ; i++ )
    for(unsigned int j = 0 ; j <= N; j++)
      cont += DP[i][j].size();
  if( verbose_flag ) printf("DP elements: %llu\n", cont);
  
  cont = 0;
  for(auto l : links)
    cont += l.second.size();
  if( verbose_flag ) printf("Oracle links: %llu\n", cont);
  
  return 0;
}
