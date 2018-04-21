
// Refactor with libraries
#include <bits/stdc++.h>

namespace gaspare 
{

    typedef unsigned long long int ull;

    typedef struct vertex_t {
        std::string seq;
        unsigned long long freq;
    } vertex_t;

    typedef struct edge_t {
        unsigned long long freq;
    } edge_t;

    // Colors handler functions
    typedef uint_fast32_t color_t;
    bool getBit(color_t n, int pos) { return ((n >> pos) & 1) == 1; }
    color_t setBit(color_t n, int pos) { return n |= 1 << pos; }
    color_t clearBit(color_t n, int pos) { return n &= ~(1 << pos); }
    color_t getCompl(color_t n, int q) { return ((1 << q) - 1) & (~n); }

    // V = typename of vertex
    // E = typename of vertex
    // Q = number of colors (default 8)
    // S = seed of rng (default 42)
    // D = direct graph (default true)
    template<typename V, typename E, unsigned int Q = 8, unsigned int S = 42, bool D = true>
    class graph 
    {
        std::vector<V> vertex;
        std::vector<std::vector<std::pair<std::size_t, E>>> edge;
        std::vector<color_t> color; 
        std::mt19937 rng;
        std::vector<std::map<color_t, ull>> M[1+Q]; // Color-Coding DP table

        color_t getRandomColor()
        {
            return rng() % Q;
        }

        public:

        graph() 
        {
            rng = std::mt19937(S);
        }

        std::size_t addVertex(V v)
        {
            std::size_t idv = vertex.size();
            vertex.push_back(v);
            edge.push_back( std::vector<std::pair<std::size_t, E>>() );
            color.push_back( getRandomColor() );
            for(unsigned i=0; i<=Q; i++)
                M[i].push_back(std::map<color_t, ull>());
            return idv;
        }

        void addEdge(std::size_t idu, std::size_t idv, E e)
        {
            assert( 0 <= idu && idu < vertex.size() );
            assert( 0 <= idv && idv < vertex.size() );
            edge[idu].push_back(std::pair<std::size_t, E>(idv, e));
            if( !D )
                edge[idv].push_back(std::pair<std::size_t, E>(idu, e));
        }

        const std::vector<V>& getVertex()
        {
            return vertex;
        }

        std::vector<std::vector<E>> getEdges()
        {
            return edge;
        }

        void processDP()
        {
            size_t N = vertex.size();

            #pragma omp parallel for schedule(guided)
            for (unsigned u = 0; u < N; u++) M[1][u][setBit(0, color[u])] = 1ll;
 
            for (unsigned i = 2; i <= Q; i++) {
                #pragma omp parallel for schedule(guided)
                for (unsigned u = 0; u < N; u++) 
                {
                    for (auto e : edge[u]) 
                    {
                        auto v = e.first;
                        for (auto d : M[i - 1][v]) 
                        {
                            if (getBit(d.first, color[u])) continue;
                            M[i][u][setBit(d.first, color[u])] += d.second;
                        }
                    }
                }
            }
        }

        std::vector<V> randomColorfulPath(ssize_t u, std::function<ssize_t (const std::vector<std::pair<std::size_t, E>>&)> chooser)
        {
            color_t cs = 0;
            std::vector<V> path;
            if( u == -1 ) return path;
            path.reserve(Q);
            path.push_back(vertex[(size_t)u]);
            cs = gaspare::setBit(cs, color[u]); 
            while(path.size() < Q)
            {
                std::vector<std::pair<std::size_t, E>> uedge;
                uedge.resize(edge[u].size());
                for( auto e : edge[u] )
                {
                    if( !getBit(cs, color[ e.first ]) )
                        uedge.push_back(e);
                }
                u = chooser(uedge);
                if( u == -1 ) break;
                cs = gaspare::setBit(cs, color[ (size_t) u]);
                path.push_back(vertex[u]);
            }
            return path;
        }
    };

    namespace genoma 
    {
        std::string getReverse(const std::string& seq)
        {
            std::string rev(seq);
            std::reverse(rev.begin(), rev.end());
            return rev;
        }

        char getComp(const char& b)
        {
            if( b == 'A' || b == 'T' ) return (b^('A'^'T'));
            if( b == 'C' || b == 'G' ) return (b^('C'^'G'));
            return b;
        }

        std::string getComp(const std::string& seq)
        {
            std::string comp;
            comp.reserve(seq.size());
            for(auto it=seq.begin(); it != seq.end(); ++it)
                comp += gaspare::genoma::getComp(*it);
            return comp;
        }

        std::string getRevComp(const std::string& seq)
        {
            return gaspare::genoma::getComp(gaspare::genoma::getReverse(seq));
        }

        std::string getSequence(const std::vector<gaspare::vertex_t>& path)
        {
            if( path.size() == 0 ) return "";
            std::string seq;
            seq.reserve(path[0].seq.size() + path.size() - 1);
            seq += path[0].seq;
            for(size_t i=1; i<path.size(); ++i)
                seq += path[i].seq[path[i].seq.size()-1];
            return seq;
        }
    }
    
    namespace sampling 
    {
        
        double jaccard(const std::set<std::string>& W, std::map<std::string, ull>& fA, std::map<std::string, ull>& fB, ull R)
        {
            ull num = 0;
            for(auto w : W)
                num += std::min(fA[w], fB[w]);
            return static_cast<double>(num) / static_cast<double>(R);
        }

        double brayCurtis(const std::set<std::string>& W, std::map<std::string, ull>& fA, std::map<std::string, ull>& fB)
        {
            ull num = 0, den = 0;
            for(auto w : W)
            {
                num += std::min(fA[w], fB[w]);
                den += fA[w] + fB[w];
            }
            return static_cast<double>(num) / static_cast<double>(den);
        }

        ssize_t pathChooser(const std::vector<std::pair<std::size_t, gaspare::edge_t>>& edges)
        {
            if( edges.size() == 0 ) return -1;
            std::vector<ull> freq;
            freq.resize(edges.size());
            for(auto e : edges)
                freq.push_back(e.second.freq);
            std::default_random_engine generator;
            std::discrete_distribution<ull> distr(freq.begin(), freq.end());
            return edges[ distr(generator) ].first;
        }

    }

}
