#include "qgraph.h"
#include <cassert>
#include <algorithm>

// This program is a simple test of the alpha, beta string encoding

using namespace DET;

int main(int argc, char** argv) {

    int nel = 4, norb = 6;
    ABStrings graph46(4, 6, true); // Verbose mode

    // Test #1: make sure that we have exactly 15 walks in the graph
    assert (graph46.nstrings == 15);
    // Test #2: loop over walks and generate strings from addresses
    for (int a = 0; a < 15; a++) {
        auto s = graph46.address2str(a);
        // Sanity check: all the orbital indeces are valid
        for (const auto &o : s) assert (o >= 0 && o < 6);
        // Another sanity check : make sure that the string is sorted
        assert (is_sorted(s.begin(), s.end()));
        // Encode the string and check if the code coincides with the original one
        int a_ = graph46.str2address(s);
        assert (a == a_);
    }
    // Test #3: similar but for ABStrings_simple
    ABStrings_simple graph46s(4, 6, true);
    assert (graph46s.nstrings == 15);
    for (int a = 0; a < 15; a++) {
        auto s = graph46s.address2str(a);
        // Sanity check: all the orbital indeces are valid
        for (const auto &o : s) assert (o >= 0 && o < 6);
        // Another sanity check : make sure that the string is sorted
        assert (is_sorted(s.begin(), s.end()));
        // Encode the string and check if the code coincides with the original one
        int a_ = graph46s.str2address(s);
        assert (a == a_);
    }
    // Test #4: edge cases
    ABStrings graph11(1,1,true);
    assert (graph11.nstrings == 1);
    // Test #5: move semantics
    auto graph46m = ABStrings(4, 6, false);
    assert (graph46m.nstrings == 15);
    // Test #6: move semantics for "simple" version of the class
    auto graph46ms = ABStrings_simple(4, 6, false);
    assert (graph46ms.nstrings == 15);

    return 0;

}
