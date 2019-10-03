#ifndef QMC1
#define QMC1

#include <tuple>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <random>

#include "qhamiltonian.h"
#include "qestimator.h"


/*
// Idea was suggested by Alex E. Doran

class Det_Lookup_ABC {
	public:
		Det_Lookup_ABC();
		virtual ~Det_Lookup_ABC() {
			// fill me in
		}
		virtual void add(const std::vector<int>&) = 0;
		virtual int lookup(const std::vector<int>&) = 0;
	private:
};

*/

class Hash_Tree_Det_Lookup {

	public:
		Hash_Tree_Det_Lookup() : value(-1), nmax(-1) { }
		Hash_Tree_Det_Lookup(int value_, int nmax_) : value(value_), nmax(nmax_) { }


		int lookup(const std::vector<int>& keys) {
			// Make sure that the vector is sorted
			if (!std::is_sorted(keys.cbegin(), keys.cend())) return -1;
			// ...................................
			auto [ nmax_changed, index ] = internal_lookup(keys, 0);
			return index;
		}

		void update_nmax(int nmax_) { nmax = nmax_; }
		int get_nmax() { return nmax; }

	private:
		std::tuple<bool, int> internal_lookup(const std::vector<int>& keys, int depth) {

			int i = keys[depth];

			if (depth == keys.size()) { 
				if (value == -1) {
					nmax++;
					value = nmax;
					return std::make_tuple(true, value); 
				} else {
					return std::make_tuple(false, value);
				}

			} else {

				if (children.find(i) == children.end()) {
					children[i] = Hash_Tree_Det_Lookup(-1, nmax); 
				} else {
					children[i].update_nmax(nmax);
				}

				auto [ nmax_changed, result ] = children[i].internal_lookup(keys, depth+1);
				if (nmax_changed) nmax = result;
				return std::make_tuple(nmax_changed, result);

			}
		}
		std::map<int, Hash_Tree_Det_Lookup> children;
		int value, nmax;

};


class FCIQMC_simple {

		typedef std::mt19937_64 random_engine; // Should be defined at compile time; will use this for now

    private:

		// Member variables 
		std::map<string, int> &par;
        vector<int> m_walker_ensemble, spawned; // weights now have signs; this should be properly accounted for

        int m_N_uniq, m_N;
        int m_steps_per_block, m_N_blocks, m_N_equil;
        double m_E_T, m_E_M;
		int init_guess_subspace;

		double dt; // imaginary time step

		// All the objects below have to be compatible with each other
        Hamiltonian &gh; // The Hamiltonian has to be compatible with the basis
		Basis &gb; 
		ProjEstimator &en_proj;
		// ...........................................................

		Hash_Tree_Det_Lookup det_index;
		random_engine *g; // Depends on the number of OpenMP threads; Each thread has its own engine

		// Member functions

        void OneTimeStep(bool equil = false);
        void initialize(bool uniform = true);
        // Update energy shift; Not used in the present code - is reserved for the production version
        void set_shift(double E_T) { m_E_T = E_T; };
		//Returns a list of determinant id-s 
		std::vector< size_t > sample_connected(const int &src, int n_samples); 

    public:

        FCIQMC_simple(std::map<string, int> &par, Hamiltonian &h, Basis &b, ProjEstimator &e);
        ~FCIQMC_simple();

        // Propagation and energy estimators

		// High level driver to handle 
		// propagation and statistics; 
		void run();

        double get_growth_estimator() { return m_E_T; }
        double get_mixed_estimator()  { return m_E_M; }

        // Probe ensemble size

        int get_num_uniq() {return m_N_uniq;}
        int get_num_total() {return m_N;}

        // Save ensemble to a text file
        // Will mute it for now 
        // void save_walkers(fstream &f);

};

#endif

