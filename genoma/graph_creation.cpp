#include <bits/stdc++.h>
#include <omp.h>
#include <parallel/algorithm>

using namespace std;

/*
1. parole originali
2. rotazioni cicle (dove la parte prima del $ e almeno k --> 150-k)
3. ordino
4. calcolo l'lcp array
5. per ogni inizio vado avanti finche lcp(i, j) < threshold
   i -> parola originale, j -> qualsiasi
*/

vector<string> reads;
vector< pair<int, int> > permutations;
vector<int> lcp;

set<int> *G;

int main(int argc, char **argv)
{
  if( argc < 2 )
  {
    printf("Usage: ./%s K\n", argv[0]);
    return 1;
  }

  ios::sync_with_stdio(false);

  unsigned k = atoi(argv[1]);

  int N;
  cin >> N;

  G = new set<int>[N];

  for(int i=0; i<N; i++)
  {
    string tmp;
    cin >> tmp;
    int oSize = tmp.size();

    tmp += "$" ; // + tmp;
    reads.push_back(tmp);

    if( oSize < k ) continue;
    for(unsigned j=0; j<oSize-k; j++)
      permutations.push_back( make_pair(j, i) );
  }

  cout << permutations.size() << endl;

  cout << "start sort" << endl;

  auto comp = [&](pair<int,int> a, pair<int,int> b){
    return strcmp( &(reads[a.second].c_str()[a.first]), &(reads[b.second].c_str()[b.first])) < 0 ;
  };

  __gnu_parallel::sort(permutations.begin(), permutations.end(), comp);

  cout << "end sort" << endl;

  lcp.reserve(permutations.size());
  lcp.push_back(0);

  cout << "start build lcp" << endl;

  #pragma omp parallel for schedule(guided)
  for(int i=1; i<permutations.size(); i++)
  {
    auto a = &(reads[permutations[i].second].c_str()[permutations[i].first]);
    auto b = &(reads[permutations[i-1].second].c_str()[permutations[i-1].first]);
    auto sa = strlen(a);
    auto sb = strlen(b);
    int pref = 0;
    while( a[pref] != '$' && b[pref] != '$' && a[pref] == b[pref] ) pref++;
    // cout << " I = " << i << " SA " << sa << " SB " << sb << " PREF " << pref << endl;
    // printf("\t[%s]\n\t[%s]\n\t%d\n\n", a, b, pref);
    lcp[i] = pref;
  }

  cout << "end build lcp" << endl;

  cout << "start building graph" << endl;

  long long deg = 0;

  #pragma omp parallel for schedule(guided)
  for(int i=0; i<permutations.size()-1; i++)
  {
    if( permutations[i].first != 0 ) continue;
    int j=i+1;
    int minPref = lcp[j];
    while( minPref >= k && j < permutations.size() )
    {
      if( permutations[i].second != permutations[j].second )
      {
        /*printf("MATCH [%d, %d][%d, %d] \n\t%s\n\t%s\n\n",
          permutations[i].second,
          permutations[i].first,
          permutations[j].second,
          permutations[j].first,
          &(reads[permutations[i].second].c_str()[0]),
          &(reads[permutations[j].second].c_str()[0])
        );*/
        // printf("%d -> %d\n", j, i);
        G[permutations[i].second].insert(permutations[j].second);
      }
      j++;
      minPref = min(lcp[j], minPref);
    }
    deg += G[permutations[i].second].size();
    G[permutations[i].second].clear();
  }

  cout << "end building graph" << endl;

//  for(int i=0; i<N; i++)
//    deg += G[i].size();

  cout << "K = " << k << " |V| = " << N << " |E| = " << deg << endl;
  return 0;
}
