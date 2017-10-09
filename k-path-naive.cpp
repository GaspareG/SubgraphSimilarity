#include <bits/stdc++.h>

using namespace std;

int N,M,k;
vector<int> *G;

int nextInt()
{
  int tmp;
  scanf("%d", &tmp);
  return tmp;
}


bool *in;
vector<int> path;

long long cont = 0;
void dfs(int a)
{
  path.push_back(a);
  in[a] = true;
  if( path.size() == (size_t)k )
  {
   cont++;
   // for(int i : path)
   //  printf("%d ", i);
   // printf("\n");
  }
  else
  {
    for(int b : G[a])
    {
      if( in[b] ) continue; 
      dfs(b);
    }
  }
  in[a] = false;
  path.pop_back();
}

int main(int argc, char **argv)
{
  if( argc < 2 ){
    printf("Usage: ./%s k\n", argv[0]);
  }

  k = atol(argv[1]);
  N = nextInt();
  M = nextInt();

  G = new vector<int>[N];
  in = new bool[N];

  for(int i=0; i<M; i++) 
  {
    int a = nextInt();
    int b = nextInt();
    G[a].push_back(b);
    G[b].push_back(a);
  }

  for(int i=0; i<N; i++) dfs(i);
  printf("%llu\n", cont);
  return 0;
}
