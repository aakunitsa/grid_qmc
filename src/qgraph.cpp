#include "qgraph.h"
#include <cassert>
#include <cstdio>

// Helper functions first

std::vector<size_t> DET::ABStrings::get_parents(const std::tuple<size_t, size_t> &pvertex) {
    auto& [k, m] = pvertex;
    std::vector<size_t> parents(0);
    size_t i;
    size_t n_if_branches = 0;
    if (k > 0 && m > 0 && k != m) {
        n_if_branches++;
        // Diagonal neighbor
        i = pair2idx(std::make_tuple(k - 1, m - 1));
        parents.push_back(i);
        // Vertical neighbor
        i = pair2idx(std::make_tuple(k - 1, m));
        parents.push_back(i);
    } else if (m == 0 && k > 0) {
        n_if_branches++;
        // Vertical neighbor only
        i = pair2idx(std::make_tuple(k - 1, m));
        parents.push_back(i);
    } else if (m == k && m != 0) {
        n_if_branches++;
        // Diagonal neighbor only
        i = pair2idx(std::make_tuple(k - 1, m - 1));
        parents.push_back(i);
    }

    assert (n_if_branches == 0 || n_if_branches == 1);

    return parents;
}

// High level calculate_* functions

size_t DET::ABStrings::calculate_vertex_weights() {
    std::vector<size_t> vertex_labels;
    for (size_t m = 0; m < (size_t)nel + 1; m++) {
        for (size_t k = m; k < (size_t)norb - (size_t)nel + m + 1; k++) 
            vertex_labels.push_back(pair2idx(std::make_tuple(k, m)));
    }
    vertex_weights[0] = 1; // Very important!
    for (const auto &l : vertex_labels) vertex_weight(l, vertex_weights);

    if (verbose) {
        printf("Vertex weights for the graph (%2d, %4d)", nel, norb);
        for (auto p = vertex_weights.begin(); p != vertex_weights.end(); p++) {
            // Decode vertex label
            auto [k, m] = idx2pair(p->first);
            printf("%zu [%zu, %zu] --> %zu\n", p->first, k, m, p->second);
        }
    }
    assert (vertex_weights.find(pair2idx(std::make_tuple(norb, nel))) != vertex_weights.end());

    return vertex_weights[pair2idx(std::make_tuple(norb, nel))];
}

void DET::ABStrings::calculate_edge_weights() {

    int k, m;

    for (auto p = vertex_weights.begin(); p != vertex_weights.end(); p++) {
        std::tie(k, m) = idx2pair(p->first); // cast as int; should work fine since all the numbers are positive
        k--; m--;
        if (k >= 0 && m >= 0) {
            edge_weights[p->first] = p->second - vertex_weights[pair2idx(std::make_tuple((size_t)k, (size_t)m))];
        } else {
            edge_weights[p->first] = 0;
        }
    }

    if (verbose) {
        // Print all the weights
        printf("Edge weights for the graph (%2d, %4d)", nel, norb);
        for (auto p = edge_weights.begin(); p != edge_weights.end(); p++) {
            auto [k, m] = idx2pair(p->first);
            printf("%zu [%zu, %zu] --> %zu\n", p->first, k, m, p->second);
        }
    }
}

// Recursive worker functions 

void DET::ABStrings::vertex_weight(size_t label, std::unordered_map<size_t, size_t> &vweights) {
    if (vweights.find(label) == vweights.end()) {
        vweights[label] = 0;
        auto parents = get_parents(idx2pair(label));
        for (auto p : parents) {
            vertex_weight(p, vweights);
            vweights[label] += vweights[p];
        }
    }
}

std::vector<size_t> DET::ABStrings::a2s(const std::tuple<size_t, size_t> &pvertex, int remaining_weight) {

    const auto& [k, m] = pvertex;
    auto ivertex = pair2idx(pvertex);
    auto parents = get_parents(pvertex);
    // Result will be saved to orbitals
    std::vector<size_t> orbitals(0);
    for (const auto &p : parents) {
        const auto& [k1, m1] = idx2pair(p);
        bool diagonal = ((k - k1) == 1 && (m - m1) == 1);
        if (diagonal) {
            int new_remaining_weight = remaining_weight - edge_weights[ivertex];
            if (new_remaining_weight > 0 && new_remaining_weight <= vertex_weights[p]) {
                auto orbitals_ = a2s(std::make_tuple(k1, m1), new_remaining_weight);
                orbitals.insert(orbitals.end(), orbitals_.begin(), orbitals_.end());
                orbitals.push_back(k - 1);
            }
        } else {
            if (remaining_weight <= vertex_weights[p]) {
                auto orbitals_ = a2s(std::make_tuple(k1, m1), remaining_weight);
                orbitals.insert(orbitals.end(), orbitals_.begin(), orbitals_.end());
            }
        }
    }
    return orbitals;
}

