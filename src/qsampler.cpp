#include "qsampler.h"
#include <algorithm> // upper_bound, lower_bound
#include <cassert>
#include <numeric>

void Uniform::single_replace(std::vector<size_t> &det) {
    // det will be modified in this function!!!
    auto nocc = det.size();
    bool accepted = false;
    size_t n1porb = (size_t)bas.get_n1porb();
    std::uniform_int_distribution<size_t> g_occ(0, nocc - 1), g_offset(1, n1porb - 1);
    auto to_remove = std::next(det.begin(), g_occ(g));
    size_t removed_orbital = *to_remove;
    det.erase(to_remove);
    while(!accepted) {
        size_t new_orbital = (removed_orbital + g_offset(g)) % n1porb;
        // Special case
        if (det.size() == 0) {
            det.push_back(new_orbital);
            accepted = true;
        // General case
        } else {
            auto it_up = std::lower_bound(det.begin(), det.end(), new_orbital);
            if (it_up != det.end() && *it_up == new_orbital)
                continue;
            else {
                accepted = true;
                det.insert(it_up, new_orbital); // will insert the orbital ___before___ the element it_up
            }
        }
    }
}

void Uniform::double_replace_simple(std::vector<size_t> &det) {
    auto nocc = det.size();
    bool accepted = false;
    size_t n1porb = (size_t)bas.get_n1porb();
    std::uniform_int_distribution<size_t> g_occ(0, nocc - 1), g_offset_occ(1, nocc - 1);
    // Generate an occupied pair
    auto o1 = g_occ(g), o2 = (o1 + g_offset_occ(g)) % nocc; 
    if (o1 > o2) std::swap(o1, o2);
    size_t removed_orbital1 = *std::next(det.begin(), o1), removed_orbital2 = *std::next(det.begin(), o2);
    det.erase(std::next(det.begin(), o1));
    det.erase(std::next(det.begin(), o2 - 1));
    std::uniform_int_distribution<size_t> u(0, n1porb - 1);
    while (!accepted) {
        // Generate two orbitals
        auto v1 = u(g), v2 = u(g);
        if (v1 == v2 || v1 == removed_orbital1 || v1 == removed_orbital2 || v2 == removed_orbital1 || v2 == removed_orbital2 || 
                find(det.begin(), det.end(), v1) != det.end() ||  find(det.begin(), det.end(), v2) != det.end()) continue;
        det.push_back(v1);
        det.push_back(v2);
        accepted = true;
        std::sort(det.begin(), det.end());
    }
}

void Uniform::double_replace(std::vector<size_t> &det) {
    // It would be good to make sure that cases with no available orbitals are handled correctly as well
    auto nocc = det.size();
    auto [na, nb] = bas.get_ab();
    bool accepted = false;
    size_t n1porb = (size_t)bas.get_n1porb(), nvir = n1porb - nocc;
    const auto &triangular_code = (nocc == nb && !na_eq_nb ? triangular_code_beta : triangular_code_alpha); 
    std::uniform_int_distribution<size_t> g_occ(0, nocc - 1), g_offset_occ(1, nocc - 1);
    // Generate an occupied pair
    auto o1 = g_occ(g), o2 = (o1 + g_offset_occ(g)) % nocc; 
    // We impose a restriction o1 < o2:
    if (o1 > o2) std::swap(o1, o2);
    size_t removed_orbital1 = *std::next(det.begin(), o1), removed_orbital2 = *std::next(det.begin(), o2);
    assert (removed_orbital1 < removed_orbital2); // Should always be the case since orbital lists are sorted
    // The occupied orbitals will be removed after the new virtuals are generated
    std::vector<size_t> ospace(n1porb);  
    std::iota(ospace.begin(), ospace.end(), (size_t)0);
    size_t offset = 0;
    for (size_t i = 0; i < nocc; i++) {
        ospace.erase(std::next(ospace.begin(), det[i] - offset));
        offset++; 
    }
    det.erase(std::next(det.begin(), o1));
    det.erase(std::next(det.begin(), o2 - 1));
    std::uniform_int_distribution<size_t> pair(0, nvir * (nvir - 1) / 2 - 1);
    assert(triangular_code.size() ==  nvir * (nvir - 1) / 2);
    auto [v2, v1] = triangular_code[pair(g)];
    assert(v1 < v2);
    det.insert(det.end(), {ospace[v1], ospace[v2]});
    std::sort(det.begin(), det.end());
}


std::tuple<double, std::vector<size_t> > Uniform::excitation(size_t det_id, size_t n_samples) {
    auto [na, nb] = bas.get_ab();
    auto [nas, nbs] = bas.get_num_str(); // Strings are stored as follows b0a0, b0a1, b0a2, ..., b1a0, b1a1, ...
    auto [ia, ib] = bas.unpack_str_index(det_id);
    const auto &det_a = bas.a(ia);
    const auto &det_b = (nb > 0 ? bas.b(ib) : std::vector<size_t>());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<size_t> samples_list; // Result
    for (size_t isample = 0; isample < n_samples; isample++) {
        // Selection is driven by single/double choice
        if (u(g) < p_s) {
            // Do singles after deciding which string to use
            if (u(g) < p_a) {
                // We will be working with alpha string here
                // 1. Make a copy of det_a
                std::vector<size_t> det_a_(na);
                std::copy(det_a.begin(), det_a.end(), det_a_.begin());
                single_replace(det_a_);
                // Temporary
                //assert(ex_order(det_a, det_a_) == 1);
                //
                samples_list.push_back(ib * nas + bas.inv_a(det_a_)); 
            } else {
                // We will be working with beta string here
                // Just in case:
                if (nb == 0) assert(false);
                std::vector<size_t> det_b_(nb);
                std::copy(det_b.begin(), det_b.end(), det_b_.begin());
                single_replace(det_b_);
                // Temporary
                //assert(ex_order(det_b, det_b_) == 1);
                //
                samples_list.push_back(bas.inv_b(det_b_) * nas + ia); 
            }
        } else {
            // Do doubles; Choose from aa, ab, bb
            int r = double_ex_type(g);
            bool aa = (r == 0), ab = (r == 2), bb = (r == 1);
            assert (aa || ab || bb);
            if (aa) {
                std::vector<size_t> det_a_(na);
                std::copy(det_a.begin(), det_a.end(), det_a_.begin());
                double_replace(det_a_);
                //double_replace_simple(det_a_);
                // Temporary
                //assert(ex_order(det_a, det_a_) == 2);
                //
                /*
                if (ex_order(det_a_, det_a) != 2) {
                    std::cout << "Error in excitation function " << std::endl;
                    std::cout << "Final : " << std::endl;
                    print_vector(det_a_);
                    std::cout << "Source : " << std::endl;
                    print_vector(det_a);
                    assert(false);
                }
                */
                samples_list.push_back(ib * nas + bas.inv_a(det_a_));
            }
            if (ab) {
                std::vector<size_t> det_a_(na);
                std::copy(det_a.begin(), det_a.end(), det_a_.begin());
                std::vector<size_t> det_b_(nb);
                std::copy(det_b.begin(), det_b.end(), det_b_.begin());
                single_replace(det_a_);
                single_replace(det_b_);
                // Temporary
                // assert(ex_order(det_a, det_a_) == 1);
                // assert(ex_order(det_b, det_b_) == 1);
                //
                samples_list.push_back(bas.inv_b(det_b_) * nas + bas.inv_a(det_a_));
            }
            if (bb) {
                std::vector<size_t> det_b_(nb);
                std::copy(det_b.begin(), det_b.end(), det_b_.begin());
                double_replace(det_b_);
                //double_replace_simple(det_b_);
                // Temporary
                //assert(ex_order(det_b, det_b_) == 2);
                /*
                if (ex_order(det_b_, det_b) != 2) {
                    std::cout << "Error in excitation function " << std::endl;
                    std::cout << "Final : " << std::endl;
                    print_vector(det_b_);
                    std::cout << "Source : " << std::endl;
                    print_vector(det_b);
                    assert(false);
                }
                */
                samples_list.push_back(bas.inv_b(det_b_) * nas + ia); 
            }
        }
    }
    assert (samples_list.size() == n_samples);
    return std::make_tuple(p_u, samples_list);
}
