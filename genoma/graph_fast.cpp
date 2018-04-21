#include <bits/stdc++.h>
#include "graph.hpp"

int main()
{
    const unsigned k = 8;
    const unsigned seed = 42;

    std::vector<std::string> reads;
    std::multiset<std::string> kmers;
    std::map<std::string, std::set<int>> kmersMap;
    std::map<std::string, std::size_t> idMap;

    gaspare::graph<gaspare::vertex_t, gaspare::edge_t, k, seed> G;

    int N;
    std::cin >> N;
    reads.reserve(N);

    for(auto i=0; i<N; i++)
    {
        std::string seq, revcseq;
        std::cin >> seq;
        revcseq = gaspare::genoma::getRevComp(seq);
        if( revcseq < seq ) seq = revcseq;
        reads.push_back(seq);
    }

    for(auto read : reads)
        for(auto i=0; i+k-1<read.size(); i++)
            kmers.insert(read.substr(i, k));

    for(auto kmer : kmers)
    {
        gaspare::vertex_t v;
        v.seq = kmer;
        v.freq = kmers.count(kmer);
        auto idx = G.addVertex(v);
        idMap[kmer] = idx;
        kmersMap[kmer.substr(0, k-1)].insert(idx);
    }

    for(auto u : kmers)
        for(auto idv : kmersMap[u.substr(1, k-1)] )
        {
            auto v = G.getVertex()[idv].seq;
            gaspare::edge_t e;
            e.freq = kmers.count(u)*kmers.count(v);
            G.addEdge(idMap[u], idMap[v], e);
        }

    G.processDP();

    auto path = G.randomColorfulPath(0, gaspare::sampling::pathChooser);
    std::cout << gaspare::genoma::getSequence(path) << std::endl;

    return 0;
}
