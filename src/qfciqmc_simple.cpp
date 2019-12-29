#include "qfciqmc_simple.h"
#include "qrandom_seed.h"
#include <vector>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_cblas.h>
#include "omp.h"


FCIQMC_simple::FCIQMC_simple(std::map<string, int> &p, std::map<string, double> &dp, Hamiltonian &h, Basis &b, Estimator &e) : gh(h), gb(b), en_proj(e), par(p) {

	// Process input parameters and/or set defaults

    m_N = p["N"];  // set N to target number specified by user in the input file
    m_N_uniq = 0;
    m_steps_per_block = p["steps_per_block"];
    m_N_blocks = p["N_blocks"];
    m_N_equil = p["N_equil"];
    if (m_N_equil < 0) m_N_equil = int(0.25 * m_N_blocks); // 25 % of the number of blocks by default

    dt = dp["dt"];
    B = dp["B"];

    // Retrieve information about the basis set & initialize big array to zeros
    m_walker_ensemble.resize(gb.get_basis_size());
    std::fill(m_walker_ensemble.begin(), m_walker_ensemble.end(), 0);

	// Determine the maximum number of OpenMP threads and set up
	// appropriate number of random engines
#ifdef _OPENMP
	int max_threads = omp_get_max_threads(); // Check what would happen by default!!
#else
	int max_threads = 1;
#endif
	g = new random_engine[max_threads];
        auto seeding_algorithm = p["seeding_algorithm"];
	for (int iengine = 0; iengine < max_threads; iengine++) {
            if (seeding_algorithm == 0) {
		g[iengine] = Rand_seed<random_engine>(simple).setup_engine();
            } else if (seeding_algorithm == 1) {
		g[iengine] = Rand_seed<random_engine>(gnu_fortran).setup_engine();
            } else if (seeding_algorithm == 2){
		g[iengine] = Rand_seed<random_engine>(sequence, iengine).setup_engine();
            } else if (seeding_algorithm >= 1000) {
                // This has been added for debugging; allows to reproduce the sequence of the random numbers used in the run
                g[iengine] = random_engine(seeding_algorithm * (iengine + 1));
                std::cout << " Seed for thread # " << iengine  << " is " << seeding_algorithm * (iengine + 1) << std::endl;
            } else {
                int seed = chrono::system_clock::now().time_since_epoch().count();
                g[iengine] = random_engine(seed);
                std::cout << " Seed for thread # " << iengine << " is " << seed << std::endl;
            }
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
        if (p["save_hamiltonian"] > 0) gh.save_matrix();
        debug = p["fciqmc_debug_mode"] > 0 ? true : false;
        power_method = p["fciqmc_power_method"] > 0 ? true : false;
        std::cout << " Power method variable is set " << std::endl;

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
#ifdef _OPENMP
				int tid = omp_get_thread_num();
				int idx = init_distr(g[tid]);
#else
				int idx = init_distr(g[0]);
#endif
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
#ifdef _OPENMP
			int tid = omp_get_thread_num();
			int tr_det = init_distr2(g[tid]);
#else
			int tr_det = init_distr2(g[0]);
                        //std::cout << tr_det << std::endl;
#endif
			int sign = guess_wfn[tr_det] > 0 ? 1 : -1;
#pragma omp critical 
			{
				m_walker_ensemble[tr_gb.get_id(tr_det)] += sign;
			}
		}

		m_N_uniq = 0; 
		for (const auto &w : m_walker_ensemble)
			m_N_uniq += abs(w); // Check if this is correct; in the present form it seems m_N_uniq == m_N

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
	printf(" Imaginary time-step = %20.10f\n", dt);
	printf(" Initial guess for the energy offset = %10.6f\n", m_E_T);

        // Printing the initial walker ensemble (this will be commented out later
        //std::cout << " Initial walker ensemble for reference purposes " << std::endl;
        //for (size_t jb = 0; jb < gb.get_basis_size(); jb++) std::cout << " On b.f. # " << jb << " : " << m_walker_ensemble[jb] << std::endl;
        //std::cout << " Control random number (to test the state of generator) " << g[0]() << std::endl;
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

        if (power_method) {

            std::cout << "Setting up a power iteration run " << std::endl;

            // Calculate the matrix and run power iteration with FCIQMC projector 
            size_t basis_size = gb.get_basis_size();
            std::vector<double> P_mat(basis_size * basis_size, 0.0), H_mat(basis_size * basis_size, 0.0);
            for (size_t i = 0; i < basis_size; i++)
                for (size_t j = i; j < basis_size; j++) {
                    double Hij = gh.matrix(i, j);
                    if (i == j) {
                        H_mat[i * basis_size + i] = Hij;
                        P_mat[i * basis_size + i] = 1. + dt * (m_E_T - Hij);
                    } else {
                        H_mat[i * basis_size + j] = Hij;
                        H_mat[j * basis_size + i] = Hij;
                        P_mat[i * basis_size + j] = -dt * Hij;
                        P_mat[j * basis_size + i] = -dt * Hij;
                    }

                }

            std::vector<double> v_old (m_walker_ensemble.begin(), m_walker_ensemble.end()),
                                v_new (basis_size, 0.0);

            // Normalize v_old
            {
                double norm = cblas_dnrm2((int)basis_size, v_old.data(), 1);
                for (auto &c : v_old) c /= norm;
            }

            printf( ">>>>>>>>> Running FCIQMC power iteration <<<<<<<<<<\n");
            printf( "block #   E_m (mixed #1) E_m (mixed #2)  E_r (Reileigh)\n");
            printf( "=======   ============== ==============  ==============\n");
            for (size_t iblock = 0; iblock < m_N_blocks; iblock++) {
                // Iterate
                cblas_dsymv(CblasRowMajor, CblasUpper, (int)basis_size, 1.0, P_mat.data(), (int)basis_size, v_old.data(), 1, 0.0, v_new.data(), 1);
                // Normalize
                double norm = cblas_dnrm2((int)basis_size, v_new.data(), 1);
                for (auto &c : v_new) c /= norm;
                // Report
                // First method
                auto [e1, e2] = en_proj.eval(v_new);
                assert (abs(e2) >= 1e-10);
                double E_first = e1/ e2, E_second = 0.0, num = 0.0, denom =0.0;
                // Second method
                for (size_t i = 0; i< basis_size; i++) {
                    auto [e1_, e2_] = en_proj.eval(i);
                    num += v_new[i] * e1_;
                    denom += v_new[i] * e2_;
                }
                assert (abs(denom) >= 1e-10);
                E_second = num/ denom;
                num = 0.0; denom = 0.0;
                double E_r = 0.0;
                for (size_t i = 0; i < basis_size; i++) {
                    denom += v_new[i] * v_new[i];
                    for (size_t j = 0; j < basis_size; j++) {
                        num += v_new[i] * H_mat[i * basis_size + j] * v_new[j];
                    }
                }
                assert (abs(denom) >= 1e-10);
                E_r = num/denom;
                printf( "%-7zu %-13.6f %-13.6f %-13.6f\n", iblock, E_first, E_second, E_r); 
                // Copy
                std::copy(v_new.begin(), v_new.end(), v_old.begin());
            }

        } else {
            // Running statistics (GSL)
            gsl_rstat_workspace *rstat_m = gsl_rstat_alloc(), *rstat_g = gsl_rstat_alloc();
			int nthreads = 1;
#ifdef _OPENMP
			nthreads = omp_get_max_threads();
#endif
            printf( ">>>>>>>>> Running FCIQMC calculation of %d threads <<<<<<<<<<\n", nthreads);

            // Starting equilibration run here (dt will not be adjusted)
            printf( "block #  total pop.  E_m (mixed)  E_g (growth)  <E_m>   <E_g>\n");
            printf( "=======  ==========  ===========  ============  =====   =====\n");
            printf( "---------------------- Equilibration run --------------------\n");
            std::cout.flush();

            for (size_t iblock = 0; iblock < m_N_equil; iblock++) {
		int N_after, N_before = m_N;
		run_block(m_steps_per_block, true, debug); // The second parameters indicates equilibration
		N_after = m_N; 
		if (N_after == 0) {
			std::cout << "All walkers died! Aborting..." << std::endl; 
			exit(EXIT_FAILURE);
		}
                m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before);
                printf( "%-7zu %-10d %-13.6f %-13.6f %-13s %-13s\n", iblock, get_num_total(), m_E_M, m_E_T, "-", "-"); 
                std::cout.flush();
            }

	// Production run
	
            printf( "---------------------- Production run ------------------------\n");
            std::cout.flush();

            for (size_t iblock = 0; iblock < m_N_blocks; iblock++) {
		int N_after, N_before = m_N;
		run_block(m_steps_per_block, false, debug);
		N_after = m_N; 
		if (N_after == 0) {
			std::cout << "All walkers died! Aborting..." << std::endl; 
			exit(EXIT_FAILURE);
		}
		m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before);
            gsl_rstat_add(m_E_M, rstat_m);
            gsl_rstat_add(m_E_T, rstat_g);
            printf( "%-7zu %-10d %-13.6f %-13.6f %-13.6f %-13.6f\n", iblock, get_num_total(), m_E_M, m_E_T, gsl_rstat_mean(rstat_m), gsl_rstat_mean(rstat_g)); 
            std::cout.flush();
            }

            gsl_rstat_free(rstat_m);
            gsl_rstat_free(rstat_g);
    }


}


void FCIQMC_simple::run_block(size_t nsteps, bool equil, bool debug_mode) {

    int basis_size = int(gb.get_basis_size());
    assert (m_walker_ensemble.size() == basis_size);
    const int N_0 = m_N;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<int> disp_walker(1, basis_size - 1); // This should be shared among threads
    int total_spawned = 0, anti_creations = 0, total_killed = 0, total_cloned = 0, N = 0, N_pr = 0, N_uniq = 0;
    double e_num = 0.0, e_denom = 0.0;
    //std::cout << " Control rn (from mt) : " << g[0]() << std::endl;
    //for (size_t dummy = 0; dummy < 1000; dummy++) std::cout << " Control rn : " << disp_walker(g[0]) << std::endl;
    //for (size_t dummy = 0; dummy < 1000; dummy++) std::cout << " Control rn : " << u(g[0]) << std::endl;

    std::stringstream int_rand_nums, double_rand_nums, double_rand_nums1;
    double_rand_nums << std::setw(20) << std::setprecision(10);
    double_rand_nums1 << std::setw(20) << std::setprecision(10);

    if (debug_mode) {
        int_rand_nums << " Determinant occupation numbers: " << std::endl;
        for (size_t i = 0; i < basis_size; i++) {
            int_rand_nums << m_walker_ensemble[i] << std::endl;
        }
        int_rand_nums << " *** End of determinant occupation numbers *** " << std::endl;
    }

#ifdef DEBUG
    clock_t t0 = clock();
	std::cout << "Starting loop..." << std::endl;
#endif

    #pragma omp parallel 
    {

        for (size_t step_in_block = 0; step_in_block < nsteps; step_in_block++) {

            std::fill(spawned.begin(), spawned.end(), 0); // resetting spawned array before making another iteration
            total_spawned = 0, anti_creations = 0, total_killed = 0, total_cloned = 0, N = 0, N_pr = 0, N_uniq = 0;
            e_num =0.0; e_denom = 0.0;

        #pragma omp for schedule(dynamic) reduction(+:anti_creations,total_spawned,total_cloned,total_killed, N_pr)
        for (int i = 0; i < basis_size; i++) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
			int tid = 0;
#endif

            const int n_walkers = abs(m_walker_ensemble[i]);
            const int sign_ref = ( m_walker_ensemble[i] > 0 ? 1 : -1);
            if (n_walkers == 0 && !debug_mode) continue; // This is important since we are working with a full population vector
            //std::cout << " On det. # " << i << std::endl;
            //N_pr += n_walkers;
            if (debug_mode) {
                for (size_t w = 0; w < n_walkers; w++) {
                    N_pr++;
                    auto di = disp_walker(g[tid]);
                    int_rand_nums << di << std::endl;
                    size_t j = size_t ((i + di) % basis_size);
                    double h = gh.matrix(i, j);
                    int sign = (h >= 0 ? 1 : -1) * (-1) * sign_ref;
                    double ps = abs(h) * dt * (basis_size - 1); 
                    int survivors = int(ps);
                    double rn = u(g[tid]);
                    double_rand_nums << rn << std::endl;
                    if (ps - survivors > rn) 
                        survivors++;
                    // Add to spawned
                    #pragma omp atomic 
                    total_spawned += survivors;
                    #pragma omp atomic
                    spawned[j] += survivors * sign;
                }

            } else {
                for (size_t w = 0; w < n_walkers; w++) {
                    N_pr++;
                    auto di = disp_walker(g[tid]);
                    //std::cout << di << std::endl;
                    size_t j = size_t ((i + di) % basis_size);
                    //int di;
                    //#pragma omp critical
                    //di = disp_walker(g[0]);
                    //size_t j = size_t ((i + di) % basis_size);
                    double h = gh.matrix(i, j);
                    int sign = (h >= 0 ? 1 : -1) * (-1) * sign_ref;
                    double ps = abs(h) * dt * (basis_size - 1); 
                    int survivors = int(ps);
                    if (ps - survivors > u(g[tid])) survivors++;
                    // Add to spawned
                    #pragma omp atomic 
                    total_spawned += survivors;
                    #pragma omp atomic
                    spawned[j] += survivors * sign;
                }
            }
            //std::cout << "Finished spawning for walker # " << i << endl;

            // Birth-death process
            double rate = dt * (gh.matrix(i, i) - m_E_T) * n_walkers;
            int nkill = int(abs(rate));

            if (debug_mode) {
                double rn = u(g[tid]);
                double_rand_nums1 << rn << std::endl;
                if((abs(rate) - nkill) > rn) nkill++;
            } else {
                if((abs(rate) - nkill) > u(g[tid])) 
    //            double rn;
    //#pragma omp critical
    //            rn = u(g[0]);
    //            if((abs(rate) - nkill) > rn) 
                    nkill++;
            }
            #pragma omp critical 
            {
                if (abs(rate) / n_walkers  > 2) {
                    std::cout << "Death rate per walker is " << abs(rate) / n_walkers << std::endl;
                    std::cout << "Consider decreasing the time step! " << std::endl;
                }
            }

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
        // -----------------------------------------------------------------------------------
        // The following code block is added for debugging purposes and will be commented out 
        // in the final version of the code
#pragma omp single 
        {
            /*
            std::cout << " >> Spawned walker array (only non-zero elements) << " << std::endl;
            for (size_t jb = 0; jb < gb.get_basis_size(); jb++) 
                if (spawned[jb] != 0) std::cout << " On b.f. # " << jb << " : " << spawned[jb] << std::endl;
            std::cout << " >> Main walker array (only non-zero elements) << " << std::endl;
            for (size_t ib = 0; ib < gb.get_basis_size() ; ib++)
                if (m_walker_ensemble[ib] != 0) std::cout << " On b.f. # " << ib << " : " << m_walker_ensemble[ib] << std::endl; 
            */
        }
        // -----------------------------------------------------------------------------------

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

#pragma omp barrier

        if (step_in_block == nsteps - 1 && !equil) {

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
    }
    if (debug_mode) {
        std::cout << " Lists of random numbers generated when running the block (valid only if run with one thread) " << std::endl;
        std::cout << int_rand_nums.str() << std::endl;
        std::cout << double_rand_nums.str() << std::endl;
        std::cout << "Second batch of doubles (for diagonal step)" << std::endl;
        std::cout << double_rand_nums1.str() << std::endl;
        std::cout << " ******************************************************************************************** " << std::endl;
    }
}

void FCIQMC_simple::save_walkers(fstream &out) {
    assert (out.is_open());
    auto n_bf = gb.get_basis_size();
    for (int i = 0; i < n_bf; i++) {
        out << m_walker_ensemble[i] << std::endl;;
    }
}

FCIQMC_simple::~FCIQMC_simple() {
	// Deallocate g here!!!
	delete [] g;
    m_walker_ensemble.clear();
}
