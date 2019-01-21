#include <bits/stdc++.h>

using ull = uint64_t;

std::optional<std::string> getRead(std::istream &in) 
{
  std::string data;
  while (in) 
  {
    std::getline(in, data);
    if( !in || data[0] == '@' ) 
        continue;
    std::getline(in, data);
    if (data.back() == '\n')
        data.pop_back(); 
    return make_optional(data);
  }
  return std::nullopt;
}

std::vector<std::string> reads;
std::multiset<std::string_view, ull> kmers;

int main(int argc, char **argv)
{
    if( argc < 5 )
    {
        std::cout << "Usage: " << argv[0] << " k q s file1 file2 ..." << std::endl;
        return 1;
    }

    int k = atoi(argv[1]);
    int q = atoi(argv[2]);

    std::optional<std::string> tmp;

    for(int i=3; i<argc; i++)
    {
        std::ifstream fi(argv[i]);
        while( (tmp=getRead(fi)).has_value() )
            reads.push_back(tmp.value());
    }

    for(auto read : reads)
        for(size_t i=0; i+k<read.size(); i++)
            kmers.insert( std::string_view(read.data() + i, k) );

    return 0;
}
