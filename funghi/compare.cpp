#include <bits/stdc++.h>

typedef uint64_t ull;

double jaccard(std::set<std::string>& W, std::map<std::string, ull>& fX, std::map<std::string, ull>& fY)
{
  ull num = 0;
  ull den = 0;
  for(auto w : W)
  {
    num += std::min(fX[w], fY[w]);
    den += std::max(fX[w], fY[w]);
  }
  if(den==ull(0)) return 0.;
  return static_cast<double>(num) / static_cast<double>(den);
}

double brayCurtis(std::set<std::string>& W, std::map<std::string, ull>& fX, std::map<std::string, ull>& fY)
{
  ull num = 0;
  ull den = 0;
  for(auto w : W)
  {
    num += 2*std::min(fX[w], fY[w]);
    den += fX[w] + fY[w];
  }
  if(den==ull(0)) return 0.;
  return static_cast<double>(num) / static_cast<double>(den);
}

int main(int argc, char **argv)
{

  if( argc < 3 )
  {
    return 1;
  }

  std::ifstream f1(argv[2]);
  std::ifstream f2(argv[3]);

  std::set<std::string> W;
  std::map<std::string, ull> fX, fY;

  while(f1)
  {
    std::string s;
    ull freq;
    f1 >> s >> freq;
    W.insert(s);
    fX[s] = freq;
  }

  while(f2)
  {
    std::string s;
    ull freq;
    f2 >> s >> freq;
    W.insert(s);
    fY[s] = freq;
  }

  double j = jaccard(W, fX, fY);
  double bc = brayCurtis(W, fX, fY);

  std::cout << argv[2] << " " << argv[3] << std::endl;
  std::cout << "Jaccard: " << j << std::endl;
  std::cout << "brayCurtis: " << bc << std::endl;

  return 0;
}
