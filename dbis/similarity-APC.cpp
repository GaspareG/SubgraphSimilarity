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
unsigned int w = 5000;      // size of the samples
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
vector<int> PC[MAXN]; // paper -> { conference }
vector<int> PA[MAXN]; // paper -> { author }PPP
vector<int> AP[MAXN]; // author -> { paper }

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
  vector<int> author_conf[MAXN];
  vector<int> author_conf_samp[MAXN];


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
    char line[1024];
    if ( NULL == fgets(line, 1023, conf) )
      break;
    int x;
    assert(1 == sscanf(line, "%d",&x));
    conf_id.push_back(x);
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
    PC[x].push_back(y);
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
    AP[y].push_back(x);
  }

  int auth_node = 113755;

  // printf("Creation of multiset\n");
  for(int a : author_id)
    for(int p : AP[a])
      for(int c : PC[p])
      author_conf[a].push_back(c);

  printf("Sampling\n");
  for(int a : author_id)
    sample(author_conf[a].begin(), author_conf[a].end(), back_inserter(author_conf_samp[a]), w, mt19937{random_device{}()});


  printf("%s\n", map_author[auth_node].c_str());
  fflush(stdout);
  vector< tuple<double, int, int> > top;


  for(int c : author_id)
    top.push_back( make_tuple( BC(author_conf_samp[auth_node], author_conf_samp[c]), auth_node, c) );

  sort(top.begin(), top.end());
  reverse(top.begin(), top.end());

  for(int i=0; i<20; i++)
  {
    printf("%.8lf\t%s\t%s\n", get<0>(top[i]), map_author[get<1>(top[i])].c_str(), map_author[get<2>(top[i])].c_str() );
    fflush(stdout);
  }

  printf("\n");
  // Closing fiels
  fclose(author);
  fclose(paper);
  fclose(conf);

  fclose(conf_paper);
  fclose(paper_author);

  return 0;
}
