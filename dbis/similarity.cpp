#include <bits/stdc++.h>
#include <omp.h>
#include <iostream>
#include <random>
#include <string>
#include <iterator>
#include <algorithm>

using namespace std;

// CONSTANT DEFINE
#define MAXN 1080000
#define MAXQ 15
#define TOPK 10
#define COLORSET uint32_t

// Flags for experiments
bool classmateFlag = false; // Test classmate x facebook
bool familyFlag = false;    // Test family x facebook
bool schoolFlag = true;     // Test school x linkedin
bool workFlag = false;      // Test work x linkedin

bool exact = true;         // true -> exact value || false -> color-coding + sampling
unsigned int w = 2000;      // size of the samples
unsigned int seed = 42;     // random seed

// File names
char filename_author[] = "author.txt";
char filename_paper[] = "paper.txt";
char filename_conf[] = "conf.txt";

char filename_paper_author[] = "paper_author.txt";
char filename_paper_conf[] = "paper_conf.txt";

// exact path definition
typedef long long ll;

// File pointers
FILE *author, *paper, *conf;
FILE *conf_paper, *paper_author;

// Graph definition
unsigned int N, E;    // number of nodes, edges
int label[MAXN];      // id node -> id label
vector<int> CP[MAXN]; // conference -> { paper }
vector<int> PA[MAXN]; // paper -> { authoro }

double BC(vector<int> Lx, vector<int> Ly)
{

  sort(Lx.begin(), Lx.end());
  sort(Ly.begin(), Ly.end());

  vector<int> Lxy(Lx.size() + Ly.size());
  auto it =
      set_intersection(Lx.begin(), Lx.end(), Ly.begin(), Ly.end(), Lxy.begin());
  Lxy.resize(it - Lxy.begin());
  double out = (double)((double)2. * Lxy.size()) / (Lx.size() + Ly.size());
  return out;
}

int main() {

  // printf("Opening files\n");

  assert(NULL != (author = fopen(filename_author,"r")));
  assert(NULL != (paper = fopen(filename_paper,"r")));
  assert(NULL != (conf = fopen(filename_conf,"r")));

  assert(NULL != (conf_paper = fopen(filename_paper_conf,"r")));
  assert(NULL != (paper_author = fopen(filename_paper_author,"r")));

  // Random generator
  srand(seed);

  vector<int> author_id, conf_id, paper_id;
  vector<int> conf_author[MAXN];
  vector<int> conf_author_samp[MAXN];
  map<int, string> map_author, map_conf;

  // printf("Read author\n");
  while( !feof(author) )
  {
    int x;
    fscanf(author, "%d\t", &x);
    char line[1024];
    if ( NULL == fgets(line, 1023, author) )
      break;
    line[strlen(line)-2] = '\0';
//    printf("%d -> %s\n",x, line);
    author_id.push_back(x);
    map_author[x] = string(line);
  }

  // printf("Read conference\n");
  while( !feof(conf) )
  {
    int x;
    fscanf(conf, "%d\t", &x);
    char line[2024];
    if ( NULL == fgets(line, 2023, conf) )
      break;
    line[strlen(line)-2] = '\0';
    // printf("%d -> %s\n",x, line);
    conf_id.push_back(x);
    map_conf[x] = string(line);
    // printf("[%d] [%s]\n", x, line);
  }

  // printf("Read paper\n");
  while( !feof(paper) )
  {
    char line[1024];
    if ( NULL == fgets(line, 1023, paper) )
      break;
    int x;
    assert(1 == sscanf(line, "%d",&x));
    paper_id.push_back(x);
  }

  // printf("Read conf-paper\n");
  while( !feof(conf_paper) )
  {
    char line[1024];
    if ( NULL == fgets(line, 1023, conf_paper) )
      break;
    int x, y;
    assert(2 == sscanf(line, "%d\t%d",&x,&y));
    CP[y].push_back(x);
  }

  // printf("Read paper-author\n");
  while( !feof(paper_author) )
  {
    char line[1024];
    if ( NULL == fgets(line, 1023, paper_author) )
      break;
    int x, y;
    assert(2 == sscanf(line, "%d\t%d",&x,&y));
    // printf("%d\n",x);
    PA[x].push_back(y);
  }

  // scanf("%d", &conf_node);
  /*
    SIGMOD 3329
    VLDB 3594
    ICDE 1798
    PODS 3027
    EDBT 1234
    DASFAA 983
    KDD 2504
    ICDM 1801
    PKDD 3011
    SDM 3230
    PAKDD 2934
    WWW 3771
    SIGIR 3318
    TREC 3510
    APWeb 294
  */
  int conf_nodes[] = {
    3011
  /*    3329,
      3594,
      1798,
      3027,
      1234,
       983,
      2504,
      1801,
      3011,
      3230,
      2934,
      3771,
      3318,
      3510,
       294*/
  };

  ios::sync_with_stdio(true);

  // printf("Creation of multiset\n");
  for(int conf_node : conf_nodes)
  {

    for(int c : conf_id)
      for(int p : CP[c])
        for(int a : PA[p])
            conf_author[c].push_back(a);

    for(int c : conf_id)
      sample(conf_author[c].begin(), conf_author[c].end(), back_inserter(conf_author_samp[c]), w, mt19937{random_device{}()});
    //    conf_author_samp[c] = conf_author[c];

    printf("%s\n", map_conf[conf_node].c_str());
    fflush(stdout);
    vector< tuple<double, int, int> > top;

    for(int c : conf_id)
      top.push_back( make_tuple( BC(conf_author_samp[conf_node], conf_author_samp[c]), conf_node, c) );

    sort(top.begin(), top.end());
    reverse(top.begin(), top.end());

    for(int i=0; i<20; i++)
    {
      // cout << get<0>(top[i]) << map_conf[get<1>(top[i])] << " " << map_conf[get<2>(top[i])] << "\n
      printf("%.8lf\t%s\t%s\n", get<0>(top[i]), map_conf[get<1>(top[i])].c_str(), map_conf[get<2>(top[i])].c_str() );
      fflush(stdout);
    }

    printf("\n");
  }
  //
  //
  // printf("\t");
  // for(int x : conf_id)
  //   printf("%d\t", x);
  // printf("\n");
  //
  // for(int x : conf_id)
  // {
  //   printf("%d\t", x);
  //   for(int y : conf_id)
  //     printf("%.6f\t", BC(conf_author_samp[x], conf_author_samp[y]));
  //   printf("\n");
  // }

  // Closing fiels
  fclose(author);
  fclose(paper);
  fclose(conf);

  fclose(conf_paper);
  fclose(paper_author);

  return 0;
}
