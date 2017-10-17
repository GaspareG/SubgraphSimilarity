/*
  Author: Gaspare Ferraro
  Count the frequency of labels sequency of length k
  start/finish in every node using color-coding technique (parallel version)
  (by Andrea Marino)

  h: V -> [n]    // etichetta
  path P = v1, ..., vr
  label H(P) = h(v1), ..., h(vr)

  t € [n]^k (stringa casuale)

  Create empty dictionaries dict[i][u] and freq[u]

  freq[N] path -> freq
  dict[k][N] path -> <colore, h, label>

  for u in V:
    dict[0][u].put(<u, color(u), H(u), label(u)>)

    < path, color(path), H(path), label(P) >


  for i in [k]:
  	for u in [V]:
  		for v in N(u):
  			for path in dict[i-1][v] such that color(u) NOT IN color(path) and H(path+u) <= t:
  				f =freq[v].lookup(label(path))
  				f'=freq[u].lookup(label(path+.u))
  				freq[u].put(label(path+u), f+f')
  				dict[i][u].put(path+u)

//Note that \leq_L is the lexicographic order and a.b means the concatenation between a and b, \chi(u) is the color of u

Sketches[u]=freq[u].set() //Skecthes[u] are the pairs <\pi,f> in freq[u]
*/

#include <vector>
#include <set>
#include <map>
#include <string>
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
unsigned int k = 0, kp = 0;
int x = -1;
unsigned int thread_count = 0;
static int verbose_flag, help_flag;

ll cont = 0;
char *labels;
int *color;
vector<int> *G;

inline int nextInt() {
  int r;
  scanf("%d", &r);
  return r;
}

// Get pos-th bit in n
inline bool getBit(COLORSET n, int pos) { return ((n >> pos) & 1) == 1; }

// Set pos-th bit in n
inline COLORSET setBit(COLORSET n, int pos) { return n |= 1 << pos; }

// Reset pos-th bit in n
inline COLORSET clearBit(COLORSET n, int pos) { return n &= ~(1 << pos); }

// Random coloring graph using kp color
inline void randomColor() {
  for (unsigned int i = 0; i < N; i++) color[i] = rand() % kp;
}

// Hashing function
int p1 = 23;
int p2 = 29;
inline int h(int x)
{
  return (p1*x+p2) % N;
}

vector<int> H(vector<int> X)
{
  vector<int> ret;
  for(auto x : X) ret.push_back(h(x));
  return ret;
}

// Dynamic programming processing
map< vector<int>, tuple<COLORSET, vector<int>, string> > *DP[MAXK+1];
map< string, ll > *freq;

void processDP() {
  if( verbose_flag ) printf("K = %d\n", 0);
  for(unsigned int u=0; u<N; u++)
  {
    vector<int> vu;
    vu.push_back(u);
    DP[0][u][vu] = make_tuple( setBit(0, color[u]), H(vu), string(&labels[u],1));
    freq[u][string(&labels[u],1)] = 1;
  }

  for(unsigned int i=1; i<=k; i++)
  {
    if( verbose_flag ) printf("K = %d\n", i);
    #pragma omp parallel for
    for(unsigned int u=0; u<N; u++)
    {
      for(auto v : G[u])
      {
        for(auto e : DP[i-1][v])
        {
          vector<int> p = e.first;
          COLORSET c = get<0>(e.second);
          vector<int> h = get<1>(e.second);
          string l = get<2>(e.second);

          if( getBit(c, color[u]) ) continue;

          vector<int> pu = vector<int>(p);
          pu.push_back(u);
          COLORSET cu = setBit(c, color[u]);
          vector<int> hu = H(pu);
          string lu = string(l);
          lu += labels[u];

          ll f = freq[v][l];
          ll fp = freq[u][lu];
          freq[u][lu] = f+fp;
          DP[i][u][pu] = make_tuple(cu, hu, lu);
        }
      }
    }
  }
}


void print_usage(char *filename) {
  //  printf("Usage: ./%s -k length -K number -g filename -f format -t filename
  //  -T filename -p threadcount -h -v\n",filename);
  printf(
      "Usage: ./%s -k length -K number -g filename -p threadcount -n index"
      "--help --verbose\n",
      filename);
  printf("Valid arguments:\n");

  printf("-k, --path length\n");
  printf("\tLength of the path.\n");

  printf("-K, --color number\n");
  printf("\tNumber of colors to use (default path length).\n");

  printf("-g, --input filename\n");
  printf("\tInput file of graph in .nme.bin format (default stdin)\n");

  /*printf("-f, --format format\n");
  printf("\tFormat of input file (snap, nde, nme)\n");

    printf("-t, --tableout filename\n");
    printf("\tOutput DP table (default stdout)\n");

    printf("-T, --tablein filename\n");
    printf("\tImport DP table (default stdin)\n");

  printf("-l, --list filename\n");
  printf("\tList k-path (default stdout)\n");
  */
  printf("-p, --parallel threadcount\n");
  printf("\tNumber of threads to use (default maximum thread avaiable)\n");

  printf("-n, --node index\n");
  printf("\tIndex of the node to print sequence and frequency\n");

  printf("--help\n");
  printf("\tDisplay help text and exit.\n");

  printf("--verbose\n");
  printf("\tPrint status messages.\n");
}

bool input_graph_flag = false;
char *input_graph = NULL;
/*
bool table_in_flag = false;
char *table_in = NULL;

bool table_out_flag = false;
char *table_out = NULL;

bool list_path_flag = false;
char *list_path = NULL;

bool format_name_flag = false;
char *format_name = NULL;
*/

int main(int argc, char **argv) {
  static struct option long_options[] = {
      {"path", required_argument, 0, 'k'},
      {"color", required_argument, 0, 'K'},
      {"input", required_argument, 0, 'g'},
      //{"format", required_argument, 0, 'f'},
      // {"tableout",required_argument,             0, 't'},
      // {"tablein", required_argument,             0, 'T'},
      //{"list", required_argument, 0, 'l'},
      {"parallel", required_argument, 0, 'p'},
      {"node", required_argument, 0, 'n'},
      {"help", no_argument, &help_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  while (1) {
    // c = getopt_long (argc, argv, "k:K:g:f:t:T:l:p:", long_options,
    // &option_index);
    //c = getopt_long(argc, argv, "k:K:g:f:l:p:", long_options, &option_index);
    c = getopt_long(argc, argv, "k:K:g:p:n:", long_options, &option_index);

    if (c == -1) break;

    switch (c) {
      case 'k':
        if (optarg != NULL) k = atoi(optarg);
        break;
      case 'K':
        if (optarg != NULL) kp = atoi(optarg);
        break;
      case 'n':
        if (optarg != NULL) x = atoi(optarg);
        break;
      case 'g':
        input_graph_flag = true;
        if (optarg != NULL) input_graph = optarg;
        break;
      /*case 'f':
        format_name_flag = true;
        if (optarg != NULL) format_name = optarg;
        break;
           case 't':
              table_in_flag = true;
              if( optarg != NULL ) table_in = optarg;
            break;
            case 'T':
              table_out_flag = true;
              if( optarg != NULL ) table_out = optarg;
            break;
      case 'l':
        list_path_flag = true;
        if (optarg != NULL) list_path = optarg;
        break;*/
      case 'p':
        if (optarg != NULL) thread_count = atoi(optarg);
        break;
    }
  }

  if (help_flag || argc == 1) {
    print_usage(argv[0]);
    return 0;
  }

  if (k == 0) {
    printf("Invalid or missing path length value.\n");
    return 1;
  }

  if (k > MAXK || kp > MAXK) {
    printf("k or kp to high! (max value: %d\n", MAXK);
    return 1;
  }

  if (kp == 0) kp = k;

  if (thread_count > 0 && (int)thread_count < omp_get_max_threads()) {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
  }

  if (verbose_flag) {
    printf("Options:\n");
    printf("k = %d\n", k);
    printf("kp = %d\n", kp);
    printf("thread = %d\n", thread_count);
    printf("input_graph = %s\n", input_graph != NULL ? input_graph : "stdin");
    //printf("format_name = %s\n", format_name != NULL ? format_name : "stdin");
    // printf("table_in    = %s\n", table_in    != NULL ? table_in    :
    // "stdin");
    // printf("table_out   = %s\n", table_out   != NULL ? table_out   :
    // "stdin");
    //printf("list_path   = %s\n", list_path != NULL ? list_path : "stdin");
  }
  if (verbose_flag) printf("Reading graph...\n");

  if (input_graph_flag) {
    if (input_graph == NULL) {
      printf("Input file name missing!\n");
      return 1;
    }

    int input_fd = open(input_graph, O_RDONLY, 0);
    if (input_fd == -1) {
      perror("Error opening input file");
      return 1;
    }
    read(input_fd, &N, sizeof(int));
    read(input_fd, &M, sizeof(int));

    color = new int[N];
    labels = new char[N];
    G = new vector<int>[N];

    read(input_fd, labels, N*sizeof(char));

    for(unsigned int i=0; i<N; i++)
      labels[i] += 'A';

    int ab[2];
    for (unsigned int i = 0; i < M; i++) {
      read(input_fd, ab, 2 * sizeof(int));
      G[ab[0]].push_back(ab[1]);
      G[ab[1]].push_back(ab[0]);
    }

  } else {
    // Read from stdin, nme format
    N = nextInt();
    M = nextInt();

    color = new int[N];
    labels = new char[N];
    G = new vector<int>[N];

    for(unsigned int i = 0; i<N; i++)
      labels[i] = 'A'+nextInt();

    for (unsigned int i = 0; i < M; i++) {
      int a = nextInt();
      int b = nextInt();
      G[a].push_back(b);
      G[b].push_back(a);
    }
  }

  if (verbose_flag) printf("N = %d | M = %d\n", N, M);

  if( x >= (int) N )
  {
    x = -1;
    if( verbose_flag ) printf("Invalid node index, it will be ignored!\n");
  }

  // Create DP Table
  for (unsigned int i = 0; i <= k; i++)
    DP[i] = new map< vector<int>, tuple<COLORSET, vector<int>, string> >[N];

  freq = new map<string, ll>[N];

  // Random color graph
  if (verbose_flag) printf("Random coloring graph...\n");
  randomColor();

  // Processing DP
  if (verbose_flag) printf("Processing DP table...\n");
  processDP();
  if (verbose_flag) printf("Finish processing DP table...\n");

/*
  if (verbose_flag)
  {
    for(unsigned int u = 0; u<N ; u++)
    {
      printf("[%u]: ", u);
      for(auto v : DP[k][u])
        printf("(%s, %llu) ", v.first.second.c_str(), v.second);
      printf("\n");
    }
  }
*/
  if( x != -1 )
  {
    printf("Pair<String, Frequency> from node %d:\n", x);
    for(auto v : freq[x])
    {
      printf("\tS=%s F=%llu\n", v.first.c_str(), v.second);
    }
  }
  return 0;
}
