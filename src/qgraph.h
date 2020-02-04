#ifndef QGRAPH_H
#define QGRAPH_H

#include <unordered_map>
#include <vector>
#include <tuple>
#include <cassert>
#include <utility>

namespace DET {
    class ABStrings {
        public:
            ABStrings(int nel_, int norb_, bool verbose_ = false) : nel(nel_), norb(norb_), verbose(verbose_) {
                nstrings = calculate_vertex_weights();
                calculate_edge_weights();
            }
            ABStrings() : nel(1), norb(1), verbose(false) {
                nstrings = calculate_vertex_weights();
                calculate_edge_weights();
            }
            ABStrings(const ABStrings &other) : 
                nel(other.nel), norb(other.norb), 
                verbose(other.verbose), nstrings(other.nstrings),  
                vertex_weights(other.vertex_weights), edge_weights(other.edge_weights) { }
            ABStrings(ABStrings &&other) :
                nel(other.nel), norb(other.norb), 
                verbose(other.verbose), nstrings(other.nstrings),  
                vertex_weights(std::move(other.vertex_weights)), edge_weights(std::move(other.edge_weights)) { }
            size_t nstrings; // Can be used within the Basis class
            std::vector<size_t> address2str(int id) {
                return a2s(std::make_tuple((size_t)norb, (size_t)nel), id + 1);
            }
            int str2address(std::vector<size_t> &s) {
                assert (s.size() == (size_t)nel);
                int address = 1;
                for (size_t i = 0; i < (size_t)nel; i++) {
                    auto o = s[i];
                    address += (int)edge_weights[pair2idx(std::make_tuple(o + 1, i + 1))];
                }
                return address - 1;
            }

        private:
            int nel, norb; // Note that we should always have  norb >= nel
            bool verbose;
            std::unordered_map<size_t, size_t> vertex_weights, edge_weights;
            // The following functions are meant to be called from the construnctor 
            // to populate the maps; Once the maps are set, we can use public
            // interface to perform orbitals string encoding/decoding
            size_t calculate_vertex_weights(); // Returns the total number of strings described by the graph
            void calculate_edge_weights();
            // Some auxiliary functions handling vertex encoding and graph navigation
            std::tuple<size_t, size_t> idx2pair(size_t idx) {
                // Will this work correctly if nel == norb?
                size_t m = idx % ((size_t)nel + 1);
                size_t k = (idx - m) / ((size_t)nel + 1);
                return std::make_tuple(k, m);
            }
            size_t pair2idx(const std::tuple<size_t, size_t> &p) {
                auto& [k, m] = p;
                return k * ((size_t)nel + 1) + m;
            }
            std::vector<size_t> get_parents(const std::tuple<size_t, size_t> &pvertex); // Input: vertex as a (k,m) pair; output: array of encoded parent vertices 
            // Recursive function calculating vertex weights
            // Could I remove the second argument by using "this" pointer?
            void vertex_weight(size_t vertex, std::unordered_map<size_t, size_t> &vweights); 
            // Private recursive function converting address to string; should be called from inside the address2str
            std::vector<size_t> a2s(const std::tuple<size_t, size_t> &pvertex, int remaining_weight);
    };

// The following class is memory hungry and should mostly be used for testing purposes
    class ABStrings_simple {
        public:
            ABStrings_simple(int nel_, int norb_, bool verbose_);
            ABStrings_simple() : nel(1), norb(1), verbose(false) {
                nstrings = 1;
                det_str.push_back(std::vector<size_t>({0}));
            }
            ABStrings_simple(const ABStrings_simple &other) : 
                nel(other.nel), norb(other.norb), 
                verbose(other.verbose), nstrings(other.nstrings),  
                det_str(other.det_str) {
            }
            ABStrings_simple& operator=(const ABStrings_simple &other) { 
                nel = other.nel; 
                norb = other.norb;
                verbose = other.verbose; 
                nstrings = other.nstrings;
                if (det_str.size() > 0) det_str.resize(0);
                det_str.insert(det_str.end(), other.det_str.begin(), other.det_str.end()); 
                return *this;
            }
            ABStrings_simple(ABStrings_simple &&other) :
                nel(other.nel), norb(other.norb), 
                verbose(other.verbose), nstrings(other.nstrings),  
                det_str(std::move(other.det_str)) {
            }
            size_t nstrings; // Can be used within the Basis class
            std::vector<size_t> address2str(int id);
            int str2address(std::vector<size_t> &s);
        private:
            std::vector< std::vector<size_t> > det_str;
            int nel, norb; // Note that we should always have  norb >= nel
            bool verbose;
            size_t construct_strings();
    };
}


#endif
