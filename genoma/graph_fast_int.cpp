#include <bits/stdc++.h>
#include "graph.hpp"
#include <omp.h>



int main()
{
    const unsigned k = 32;
    const unsigned seed = 42;

    std::vector<std::string> reads;
    std::multiset<gaspare::genoma::read_t> kmers;
    std::map<gaspare::read_t, std::set<int>> kmersMap;
    std::map<gaspare::read_t, std::size_t> idMap;

    gaspare::graph<gaspare::vertex_i_t, gaspare::edge_t, k, seed> G;

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
    auto kmersT = new std::multiset<gaspare::genoma::read_t>[omp_get_max_threads()];
    #pragma omp parallel for schedule(guided)
    for(size_t j=0; j<reads.size(); j++)
    {
        int t = omp_get_thread_num();
        auto read = reads[j];
        auto kmer = gaspare::genoma::getRead(read, 0, k);
        kmersT[t].insert(kmer);
        for(size_t i=k; i+k-1<read.size(); i++)
        {
            kmer = gaspare::genoma::nextKmer(kmer, read[i]);
            kmersT[t].insert(kmer);
        }
    }    

    for(int i=0; i<omp_get_max_threads(); i++)
    {
        std::cout << "map " << i << " size: " << kmersT[i].size() << std::endl;
        //kmers.merge(kmersT[i]);
    }
    
    std::cout << "Create (k-1)-mers mapping" << std::endl;
    for(auto kmer : kmers)
    {
        gaspare::vertex_i_t v;
        v.seq = kmer;
        v.freq = kmers.count(kmer);
        auto idx = G.addVertex(v);
        idMap[kmer] = idx;
        kmersMap[gaspare::genoma::substr(kmer, 1, k-1)].insert(idx);
    }

    std::cout << "Create graph' edges" << std::endl;
    for(auto u : kmers)
        for(auto idv : kmersMap[gaspare::genoma::substr(u, 1, k-1)] )
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
