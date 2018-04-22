#include <bits/stdc++.h>
#include "graph.hpp"

int main()
{
    const unsigned k = 50;
    const unsigned seed = 42;

    std::vector<std::string> reads;
    std::multiset<std::string> kmers;
    std::map<std::string, std::set<int>> kmersMap;
    std::map<std::string, std::size_t> idMap;

    gaspare::graph<gaspare::vertex_t, gaspare::edge_t, k, seed> G;

    int N;
    std::cin >> N;
    reads.reserve(N);

    std::cout << "Reading reads" << std::endl;
    for(auto i=0; i<N; i++)
    {
        std::string seq, revcseq;
        std::cin >> seq;
        revcseq = gaspare::genoma::getRevComp(seq);
        if( revcseq < seq ) seq = revcseq;
        reads.push_back(seq);
    }

    std::cout << "Creating k-mers" << std::endl;
    for(auto read : reads)
    {
        #pragma omp parallel for schedule(guided)
        for(auto i=0; i+k-1<read.size(); i++)
        {
            auto kmer = read.substr(i, k);
            #pragma omp critical
            {
              kmers.insert(kmer);
            }
        }
    }

    std::cout << "Create (k-1)-mers mapping" << std::endl;
    for(auto kmer : kmers)
    {
        gaspare::vertex_t v;
        v.seq = kmer;
        v.freq = kmers.count(kmer);
        auto idx = G.addVertex(v);
        idMap[kmer] = idx;
        kmersMap[kmer.substr(0, k-1)].insert(idx);
    }

    std::cout << "Create graph' edges" << std::endl;
    for(auto u : kmers)
        for(auto idv : kmersMap[u.substr(1, k-1)] )
        {
            auto v = G.getVertex()[idv].seq;
            gaspare::edge_t e;
            e.freq = kmers.count(u)*kmers.count(v);
            G.addEdge(idMap[u], idMap[v], e);
        }


    std::cout << "Processing DP table" << std::endl;
    G.processDP();

    std::cout << "Sampling random path" << std::endl;
    std::mt19937 rng(seed);
    for(int i=0; i<10; i++)
    {
        auto path = G.randomColorfulPath(rng() % G.getVertex().size(), gaspare::sampling::pathChooser);
        std::cout << gaspare::genoma::getSequence(path) << std::endl;
    }
    return 0;
}
