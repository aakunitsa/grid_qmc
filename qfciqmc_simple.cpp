#include "qfciqmc_simple.h"
#include "qrandom_seed.h"
#include <vector>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <chrono>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include "omp.h"


FCIQMC_simple::FCIQMC_simple(std::map<string, int> &p, std::map<string, double> &dp, Hamiltonian &h, Basis &b, ProjEstimator &e) : gh(h), gb(b), en_proj(e), par(p) {

	// Process input parameters and/or set defaults

    m_N = p["N"];  // set N to target number specified by user in the input file
    m_N_uniq = 0;
    m_steps_per_block = p["steps_per_block"];
	m_N_blocks = p["N_blocks"];
	m_N_equil = int(0.25 * m_N_blocks);

	dt = dp["dt"];
	B = dp["B"];

	// Retrieve information about the basis set & initialize big array to zeros
    m_walker_ensemble.resize(gb.get_basis_size());
	std::fill(m_walker_ensemble.begin(), m_walker_ensemble.end(), 0);

	// Determine the maximum number of OpenMP threads and set up
	// appropriate number of random engines
	int max_threads = omp_get_max_threads(); // Check what would happen by default!!
	g = new random_engine[max_threads];
	for (int iengine = 0; iengine < max_threads; iengine++) {
		g[iengine] = Rand_seed<random_engine>(gnu_fortran).setup_engine();
	}

	// Populate initial walker ensemble
	bool uniform_init = (p["fciqmc_guess_subspace"] > 0 ? false : true);
	// This will be used in initialize if uniform_init == false
	if(!uniform_init) init_guess_subspace = std::min(p["fciqmc_guess_subspace"], int(gb.get_basis_size())); 
    initialize(uniform_init);

	// Hash_Tree_Det_Lookup will be populated below
	int ndet = gb.get_basis_size(), n1porb = gb.get_n1porb();
	auto [na, nb] = gb.get_ab();

	// Loop over determinants and add them to the tree
/*	
	for (int idet = 0; idet< ndet; idet++) {
		auto [ ialpha, ibeta ] = gb.unpack_str_index(idet);
		if (nb == 0) {
			auto index = det_index.lookup(gb.a(ialpha)); // Lookup wants a vector of ints!!!!!!!!!!
			assert (index == idet);
		} else {
			std::vector<int> combined_str; // represents a complete alpha-beta string
			combined_str.insert(combined_str.end(), gb.a(ialpha).begin(), gb.a(ialpha).end());

            // Retrieve beta string and transform the orbital indeces so that they differ from 
			// those referring to the alpha orbitals

			std::vector<size_t> beta_tmp(gb.b(ibeta));
			std::transform(beta_tmp.begin(), beta_tmp.end(), beta_tmp.begin(), [&](size_t x) { return x+n1porb ; });

			combined_str.insert(combined_str.end(), beta_tmp.begin(), beta_tmp.end());
			auto index = det_index.lookup(combined_str);
			assert( index == idet);
		}
	}
*/
	// Prepare spawned 
	spawned.resize(gb.get_basis_size()); 
	// Temporary
	gh.save_matrix();

}

void FCIQMC_simple::initialize (bool uniform) {

    auto basis_size = gb.get_basis_size();

	int nwalkers = 0;

    double max_diag = std::numeric_limits<double>::min(), min_diag = std::numeric_limits<double>::max(); 
	if (uniform) {
		uniform_int_distribution<int> init_distr(0, basis_size - 1);
		#pragma omp parallel 
		{
			#pragma omp for schedule(static) 
			for (int n = 0; n < m_N; n++) {
				int tid = omp_get_thread_num();
				int idx = init_distr(g[tid]);
				#pragma omp critical
				{
					m_walker_ensemble[idx] += 1;
				}
			}
        
			#pragma omp barrier

			#pragma omp for schedule(static) reduction(max:max_diag) reduction(min:min_diag) reduction(+:nwalkers)
			for (int i = 0; i < basis_size; i++) {
				double diag_i = gh.matrix(i, i);
				max_diag = max(max_diag, diag_i);
				min_diag = min(min_diag, diag_i);
				nwalkers += m_walker_ensemble[i];
			}
		}

		assert (nwalkers == m_N);

		// I should setup a separate parameter reader for the floating point 
		// numbers

		m_E_T = min_diag - 0.1, m_E_M = min_diag;  // initial guess for the ground state energy
		m_N_uniq = 0;
		for (const auto &w : m_walker_ensemble)
			m_N_uniq += abs(w);
	} else {
		// 1. Construct a truncated basis
		auto H_diag = gh.build_diagonal();
		TruncatedBasis tr_gb(par, gb.get_n1porb(), init_guess_subspace, H_diag, gb);
		Hamiltonian tr_gh(gh.get_integral_factory(), tr_gb);
		// 2. Diagonalize Hamiltonian in truncated basis and obtain initial guess function and energy
		auto guess_en = tr_gh.diag(true);
		std::vector<double> guess_wfn = tr_gh.get_wfn();
		// 3. Use coefficients of the w.f. expansion to populate global array (use sampling!)
		double norm = 0.0;
		for (auto &c : guess_wfn) norm += abs(c);
		std::vector<double> init_distr1(tr_gb.get_basis_size(), 0.0);
		std::transform(guess_wfn.begin(), guess_wfn.end(), init_distr1.begin(), [&](double c) { return abs(c) / norm; });
		std::discrete_distribution<int> init_distr2(init_distr1.begin(), init_distr1.end());
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::default_random_engine generator (seed);
#pragma omp parallel for
		for (int i = 0; i < m_N; i++) {
			int tid = omp_get_thread_num();
			int tr_det = init_distr2(g[tid]);
			int sign = guess_wfn[tr_det] > 0 ? 1 : -1;
#pragma omp critical 
			{
				m_walker_ensemble[tr_gb.get_id(tr_det)] += sign;
			}
		}

		m_N_uniq = 0;
		for (const auto &w : m_walker_ensemble)
			m_N_uniq += abs(w);

		m_E_T = guess_en[0]; m_E_M = guess_en[0];
		max_diag = *std::max_element(H_diag.begin(), H_diag.end());
		min_diag = *std::min_element(H_diag.begin(), H_diag.end());
	}

	if (dt < 0) 
		dt = 0.1 *  1. / ( max_diag - min_diag); // See PRB 44 9410 (1991)

	std::cout << " Maximum H_ii = " << max_diag << std::endl;
	std::cout << " Minimum H_ii = " << min_diag << std::endl;
	std::cout << " Number of walkers = " << m_N << std::endl;
	std::cout << " Number of unique walkers = " << m_N_uniq << std::endl;
	printf(" Imaginary time-step = %10.6f\n", dt);
	printf(" Initial guess for the energy offset = %10.6f\n", m_E_T);

	// Print m_walker_ensemble
	//std::cout << " Printing the occupation vector" << std::endl;
	//for (auto p : m_walker_ensemble)
	//	std::cout << p << ' ';
    //std::cout << std::endl;

}

// High level driver to start the simulation

void FCIQMC_simple::run() {

	// Structure
	// 1. Equilibration
	// 2. Production run:
	// 2.1 Loop over blocks
	// 2.1.1 Loop over steps
	// 2.1.2 Update E_T
	// 2.1.3 Collect stats
	
	// Some parameters
	if (B < 0)
		B = 1.0; // default damping parameter for population control
	// Running statistics (GSL)
    gsl_rstat_workspace *rstat_m = gsl_rstat_alloc(), *rstat_g = gsl_rstat_alloc();

	printf( ">>>>>>>>> Running FCIQMC calculation of %d threads <<<<<<<<<<\n", omp_get_max_threads());

	// Starting equilibration run here (dt will not be adjusted)
    printf( "block #  total pop.  E_m (mixed)  E_g (growth)  <E_m>   <E_g>\n");
    printf( "=======  ==========  ===========  ============  =====   =====\n");
    printf( "---------------------- Equilibration run --------------------\n");
	std::cout.flush();

	for (size_t iblock = 0; iblock < m_N_equil; iblock++) {
		int N_after, N_before = m_N;
		for (size_t istep = 0; istep < m_steps_per_block; istep++) OneTimeStep(true);
		N_after = m_N; 
		if (N_after == 0) {
			std::cout << "All walkers died! Aborting..." << std::endl; 
			exit(EXIT_FAILURE);
		}
		m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before);
        printf( "%-7d %-10d %-13.6f %-13.6f %-13s %-13s\n", iblock, get_num_total(), m_E_M, m_E_T, "-", "-"); 
		std::cout.flush();
	}

	// Production run
	
    printf( "---------------------- Production run ------------------------\n");
	std::cout.flush();

	for (size_t iblock = 0; iblock < m_N_blocks; iblock++) {
		int N_after, N_before = m_N;
		for (size_t istep = 0; istep < m_steps_per_block; istep++) OneTimeStep(false);
		N_after = m_N; 
		if (N_after == 0) {
			std::cout << "All walkers died! Aborting..." << std::endl; 
			exit(EXIT_FAILURE);
		}
		m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before);
        gsl_rstat_add(m_E_M, rstat_m);
        gsl_rstat_add(m_E_T, rstat_g);
        printf( "%-7d %-10d %-13.6f %-13.6f %-13.6f %-13.6f\n", iblock, get_num_total(), m_E_M, m_E_T, gsl_rstat_mean(rstat_m), gsl_rstat_mean(rstat_g)); 
		std::cout.flush();
	}

	gsl_rstat_free(rstat_m);
	gsl_rstat_free(rstat_g);


}


void FCIQMC_simple::OneTimeStep(bool equil) {

	int basis_size = int(gb.get_basis_size());
	assert (m_walker_ensemble.size() == basis_size);
    const int N_0 = m_N;
    int total_spawned = 0, anti_creations = 0, total_killed = 0, total_cloned = 0, N = 0, N_pr = 0, N_uniq = 0;
    uniform_real_distribution<double> u(0.0, 1.0);

    double e_num = 0.0, e_denom = 0.0;

	static int steps_in_block = 0;

	std::fill(spawned.begin(), spawned.end(), 0); // resetting spawned array before making another iteration


#ifdef DEBUG
    clock_t t0 = clock();
	std::cout << "Starting loop..." << std::endl;
#endif

    #pragma omp parallel 
    {
        #pragma omp for schedule(static) reduction(+:anti_creations,total_spawned,total_cloned,total_killed, N_pr)
        for (int i = 0; i < basis_size; i++) {

            int tid = omp_get_thread_num();

            const int n_walkers = abs(m_walker_ensemble[i]);
			if (n_walkers == 0) continue; // This is important since we are working with a full population vector
            N_pr += n_walkers;
			auto neigh = gb.get_neigh(i);
			double sprob = 1. / (neigh.size() - 1); // -1 appears because determinant i is included in the list
            const int sign_ref = ( m_walker_ensemble[i] > 0 ? 1 : -1);
			//std::cout << "Doing spawning for walker " << i << std::endl;
			auto conn = sample_connected(i, n_walkers); 
            for (auto &j : conn) {
                
                // Spawning
                double h = gh.matrix(i, j);
                int sign = (h >= 0 ? 1 : -1) * (-1) * sign_ref;
                double ps = abs(h) * dt / sprob; 
                //if (ps > 2)
                //    cout << "ps = " << ps << "!" << endl; 

                int survivors = int(ps);
                if (ps - survivors > u(g[tid])) 
                    survivors++;

				// should I do pragma omp atomic?
                #pragma omp critical
                {
					// Add to spawned
					total_spawned += survivors;
					spawned[j] += survivors * sign;
                }
            }

            //std::cout << "Finished spawning for walker # " << i << endl;

            // Birth-death process
            double rate = dt * (gh.matrix(i, i) - m_E_T) * n_walkers;
            int nkill = int(abs(rate));

            if((abs(rate) - nkill) > u(g[tid])) 
                nkill++;

            if (rate >= 0) {
                int new_unsigned_weight = n_walkers - nkill;
                if (new_unsigned_weight < 0) {
                    anti_creations += abs(new_unsigned_weight);
                    total_killed += n_walkers;
                } else {
                    total_killed += nkill;
                }
                m_walker_ensemble[i] = new_unsigned_weight * sign_ref;

            } else {
				//std::cout << "Inside cloning branch" << endl;
                total_cloned += nkill;
                m_walker_ensemble[i] = (n_walkers + nkill) * sign_ref;
            }
            //std::cout << "Finished cloning/death for walker # " << i << endl;
        }

		// Should I put a BLAS call here? (i.e. calculate the range 
		// for each thread and write a proper BLAS call for it)

        #pragma omp for schedule(static) reduction(+:N,N_uniq)
        for (int i = 0; i < basis_size; i++) {
            m_walker_ensemble[i] += spawned[i];
            N += abs(m_walker_ensemble[i]);
            if (m_walker_ensemble[i] != 0) N_uniq++;
        }

		// Can be replaced with smth non-blocking (pragma omp master?)
		// but probably would not make much difference

        #pragma omp single
        {
            m_N = N;
            m_N_uniq = N_uniq;

#ifdef DEBUG
			std::cout << "Number of walkers processed in the loop = " << N_pr << std::endl;
			std::cout << "Number of unique walkers " << m_N_uniq << std::endl;
			std::cout << "Total number of walkers " << m_N << std::endl;
            int N_debug = 0, N_pos = 0, N_neg = 0;
            for (auto w : m_walker_ensemble) {
                N_debug += abs(w);
                if (w > 0)
                    N_pos += w;
                if (w < 0)
                    N_neg += abs(w);
            }

			std::cout << "Total number of walkers (debug) " << N_debug << std::endl;
			std::cout << "Total number of positive (debug) " << N_pos << std::endl;
			std::cout << "Total number of negative (debug) " << N_neg << std::endl;
			std::cout << "Number of anti-particles created = " << anti_creations << std::endl;
			std::cout << "Number of particles killed = " << total_killed << std::endl;
			std::cout << "Death rate = " << double(total_killed)/N_0 << std::endl;
			std::cout << "Number of particles cloned = " << total_cloned << std::endl;
			std::cout << "Average spawning rate = " << double(total_spawned)/N_0 << std::endl;
			std::cout << "====" << std::endl;
#endif

            if (m_N == 0) {
				std::cout << "All the walkers died! Bye-bye!\n";
                abort();
            }


            if (double(total_spawned)/N_0 > 1) {
                cout << "WARNING: the spawning rate is = " << double(total_spawned)/N_0 << std::endl
                     << "consider reducing the time step..." << std::endl;
            }

        }

        if (steps_in_block == m_steps_per_block - 1 && !equil) {

			// In order to perform the calculation of the mixed estimator
			// one has to have a trial function which in Booth et al 
			// paper was taken as a HF determinant; that will be 
			// implemented later - for now I will use truncated CI

			auto n_bf = gb.get_basis_size();
            #pragma omp for reduction(+:e_num,e_denom)
			for (size_t i = 0; i < n_bf; i++) {
                if (m_walker_ensemble[i] == 0) continue;
				auto [n_, d_] = en_proj.eval(i);
				e_num += n_ * m_walker_ensemble[i];
				e_denom += d_ * m_walker_ensemble[i];
			}
			#pragma omp single
			{
				assert (abs(e_denom) >= 1e-10); // !!! 
				m_E_M = e_num/e_denom;
			}
        }

    }


    if(!equil) {
        if (steps_in_block == m_steps_per_block - 1) {
            steps_in_block = 0;
        } else {
            steps_in_block++;
        }
    }
}


std::vector< size_t > FCIQMC_simple::sample_connected(const int &src, int n_samples) {

	std::vector< size_t > sample;
	auto neigh = gb.get_neigh(src);

	uniform_int_distribution<int> neigh_dist(0, neigh.size() - 1);

	int tid = omp_get_thread_num();

	// I will assume that orbital lists are arranged such that alpha orbitals
	// preceed beta ones; That ordering applies separately to occupied and 
	// virtual lists

	for (size_t is = 0; is < n_samples; is++) {
		bool accepted = false;
		while(!accepted) {
			auto index = neigh_dist(g[tid]);
			if (neigh[index] != src) {
				accepted = true;
				sample.push_back(neigh[index]);
			}
		}
	}

	assert (sample.size() == n_samples);

	return sample;
}


/*
void FCIQMC_simple::save_walkers(fstream &out) {
    assert (out.is_open());
    vector<Det> &det_basis = m_h.m_basis; // Need to make sure that it has been constructed!

    for (int i = 0; i < det_basis.size(); i++) {
        auto d_repr = det_basis[i].to_vec();
        int weight = m_walker_ensemble[i];

        for (auto i : d_repr)
            out << i << '\t';
        out << weight << endl;

    }
}
*/

FCIQMC_simple::~FCIQMC_simple() {
	// Deallocate g here!!!
	delete [] g;
    m_walker_ensemble.clear();
}
