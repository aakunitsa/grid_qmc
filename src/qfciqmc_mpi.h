#ifndef QMC2
#define QMC2

#include <tuple>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <string>
#include <random>
#include <iostream>
#include <mpi.h>

#include "qhamiltonian_mpi.h"
#include "qbasis.h"
#include "qestimator.h"


struct Walker {
    Walker(int det_id_, int weight_) : det_id(det_id_), weight(weight_) {}
    Walker() : det_id(0), weight(0) {}
    int det_id;
    int weight;
};

inline bool comp_less (const Walker &w1, const Walker &w2) { return w1.det_id < w2.det_id; }

class FCIQMC_mpi {
#ifdef MT64
	typedef std::mt19937_64 random_engine; // Should be defined at compile time; will use this for now
#endif
#ifdef MT32
	typedef std::mt19937 random_engine; // Should be defined at compile time; will use this for now
#endif

    private:

	// Member variables 
        int me, size;
        int local_it_count = 0; // local iteration count for debugging purposes
	std::map<string, int> &par;
        std::unordered_map<int, int> m_walker_ensemble;
        //std::map<int, int> m_walker_ensemble;
        std::vector< std::vector<Walker> > local_spawned; // Stores the determinants asigned to all ranks
        std::vector<Walker> global_spawned; // Needed to store the walkers collected from all ranks and assigned to "me"
        std::vector<int> local_n_spawned, global_n_spawned;
        //std::vector<int> disp; // For MPI_Gatherv
        std::vector<int> sdisp, rdisp; // For MPI_Alltoallv


        int m_N_global, m_N;
        int m_N_uniq_global, m_N_uniq; 
        int m_steps_per_block, m_N_blocks, m_N_equil;
        double m_E_T, m_E_M;

	int init_guess_subspace;
	double B, dt; // imaginary time step

        bool debug; // Turns on the verbose mode for run_step; see below and the implementation

	// All the objects below have to be compatible with each other
        Hamiltonian_mpi &gh; // The Hamiltonian has to be compatible with the basis
	Basis &gb; 
	Estimator &en_proj;

	random_engine g; 

	// Member functions

        void run_step(bool verbose = false);
        void initialize(bool uniform = true);
        // Update energy shift; Not used in the present code - is reserved for the production version
        void set_shift(double E_T) { m_E_T = E_T; };
	//Returns a list of determinant id-s 
	std::vector< size_t > sample_connected(const int &src, int n_samples); 
        int fnv_hash(const int &src); // Determines the assignment of the determinant
        //std::vector<Walker> compress(std::vector<Walker> &v); // Compresses a list of walkers to minimize communication
        void update_walker_lists();
        void merge_walker(Walker &new_walker, std::vector<Walker> &v) {
            auto up = std::upper_bound(v.begin(), v.end(), new_walker, comp_less);
            if (up != v.begin()) {
                auto prev = std::prev(up);
                if (prev->det_id == new_walker.det_id) {
                    prev->weight += new_walker.weight;
                } else {
                    v.insert(up, new_walker);
                }
            } else {
                v.insert(up, new_walker);
            }
        }

    public:

        FCIQMC_mpi(std::map<string, int> &par, std::map<string, double> &dpar, Hamiltonian_mpi &h, Basis &b, Estimator &e);
        ~FCIQMC_mpi();

        // Propagation and energy estimators

	// High level driver to handle propagation, statistics, and MPI communication
	void run();

        double get_growth_estimator() { return m_E_T; }
        double get_mixed_estimator()  { return m_E_M; }

        // Probe ensemble size

        int get_num_uniq() {return m_N_uniq_global;}
        int get_num_total() {return m_N_global;}

        // Save ensemble to a text file
        void save_walkers(fstream &f);

};

#endif

