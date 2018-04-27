
// Refactor with libraries
#include <bits/stdc++.h>

namespace gaspare
{

    typedef uint64_t ull;

    typedef std::pair<size_t, gaspare::ull> read_t;
        
    typedef struct vertex_t {
        std::string seq;
        ull freq;
    } vertex_t;

    typedef struct edge_t {
        ull freq;
    } edge_t;

    typedef struct vertex_i_t {
        gaspare::read_t seq;
        ull freq;
    } vertex_i_t;

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

        unsigned int getQ()
        {
            return Q;
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

        std::vector<V> randomColorfulPath(ssize_t u, std::function<ssize_t (const std::vector<std::pair<std::size_t, E>>&, std::mt19937&)> chooser)
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
                u = chooser(uedge, rng);
                if( u == -1 ) break;
                cs = gaspare::setBit(cs, color[ (size_t) u]);
                path.push_back(vertex[u]);
            }
            return path;
        }
    };

    namespace genoma
    {

        typedef std::pair<size_t, gaspare::ull> read_t;

         gaspare::ull getBaseCode(const char& b)
         {
             return ((b == 'G' || b == 'T')<<1) | (b == 'C' || b == 'T');
         }

         char getBaseChar(const ull& x)
         {
             if( (x&3) == 0 ) return 'A';
             if( (x&3) == 1 ) return 'C';
             if( (x&3) == 2 ) return 'G';
             if( (x&3) == 3 ) return 'T';
             return 'X';
         }

        std::string getStringRead(const gaspare::read_t& r)
        {
            std::string seq;
            seq.reserve(r.first);
            for(size_t i = 0; i<r.first; i += 2)
                seq += gaspare::genoma::getBaseChar(r.second >> (i<<1));
            return seq;
        }

        gaspare::read_t substr(const gaspare::read_t& r, int i, int j)
        {
            gaspare::read_t sub;
            sub.first = (j-i);
            sub.second = r.second >> (2*( r.first - j )); 
            sub.second &= (1<<(2*sub.first))-1;
            return sub;
        }

        gaspare::read_t getRead(const std::string& seq, size_t i, size_t j)
        {
            gaspare::read_t read;
            read.first = (j-i);
            read.second = 0;
            for(; i < j; i++)
                read.second = (read.second<<2) | gaspare::genoma::getBaseCode(seq[i]);
            return read;
        }

        gaspare::read_t getRead(const std::string& seq)
        {
            return gaspare::genoma::getRead(seq, 0, seq.size());
        }

        gaspare::read_t nextKmer(const gaspare::read_t& seq, const char&b)
        {
            //std::cout << "next [" << seq.second << "]" << b << " " <<  std::endl;
            gaspare::read_t ret;
            ret.first = seq.first;
            ret.second = (seq.second << 2) | gaspare::genoma::getBaseCode(b);
            // ret.second &= (1<<(2*ret.first))-1;
            return ret;
        }
                
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

        std::string getSequence(const gaspare::read_t& read)
        {
            std::string seq;
            seq.reserve(read.first);
            gaspare::ull tmp = read.second;
            for(size_t i = 0; i<read.first; i++)
            {
                seq += gaspare::genoma::getBaseChar(tmp);
                tmp >>= 2;
            }
            std::reverse(seq.begin(), seq.end());
            return seq;
        }

        std::string getSequence(const gaspare::vertex_i_t& vertex)
        {
            return gaspare::genoma::getSequence(vertex.seq);
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

        std::string getSequence(const std::vector<gaspare::vertex_i_t>& path)
        {
            if( path.size() == 0 ) return "";
            std::string ret;
            ret.resize(path[0].seq.first + path.size() - 1);
            ret += gaspare::genoma::getSequence(path[0]);
            // TODO 
            return ret;
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

        ssize_t pathChooser(const std::vector<std::pair<std::size_t, gaspare::edge_t>>& edges, std::mt19937& rng)
        {
            if( edges.size() == 0 ) return -1;
            std::vector<ull> freq;
            freq.resize(edges.size());
            for(auto e : edges)
                freq.push_back(e.second.freq);
            std::discrete_distribution<ull> distr(freq.begin(), freq.end());
            return edges[ distr(rng) ].first;
        }

    }

}
