#include "qfciqmc_mpi.h"
#include "qrandom_seed.h"
#include <vector>
#include <algorithm>
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

MPI_Datatype Walker_dt;

FCIQMC_mpi::FCIQMC_mpi(std::map<string, int> &p, std::map<string, double> &dp, Hamiltonian_mpi &h, Basis &b, Estimator &e) : gh(h), gb(b), en_proj(e), par(p) {

    local_it_count = 0;
    // Process input parameters and/or set defaults
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Prepare arrays for MPI
    local_n_spawned.resize(size);
    local_spawned.resize(size);
    std::fill(local_n_spawned.begin(), local_n_spawned.end(), 0);
    global_n_spawned.resize(size);
    std::fill(global_n_spawned.begin(), global_n_spawned.end(), 0);
    //disp.resize(size);
    sdisp.resize(size);
    rdisp.resize(size);
    //std::fill(disp.begin(), disp.end(), 0);
    std::fill(sdisp.begin(), sdisp.end(), 0);
    std::fill(rdisp.begin(), rdisp.end(), 0);
    // Register the Walker structure with MPI
    MPI_Aint displacements[2]  = {offsetof(Walker, det_id), offsetof(Walker, weight)};
    int block_lengths[2]  = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    //MPI_Datatype Walker_dt;
    MPI_Type_create_struct(2, block_lengths, displacements, types, &Walker_dt);
    MPI_Type_commit(&Walker_dt);

    m_N = p["N"] / size;  // set N to target number specified by user in the input file
    m_N_uniq = 0;
    m_steps_per_block = p["steps_per_block"];
    m_N_blocks = p["N_blocks"];
    m_N_equil = p["N_equil"];
    debug = p["fciqmc_debug_mode"] > 0 ? true : false;
    if (m_N_equil < 0) m_N_equil = int(0.25 * m_N_blocks); // 25 % of the number of blocks by default

    dt = dp["dt"];
    B = dp["B"];

    auto seeding_algorithm = p["seeding_algorithm"];
    if (seeding_algorithm == 0) {
	g = Rand_seed<random_engine>(simple).setup_engine();
    } else if (seeding_algorithm == 1) {
	g = Rand_seed<random_engine>(gnu_fortran).setup_engine();
    } else if (seeding_algorithm == 2){
	g = Rand_seed<random_engine>(sequence, me).setup_engine();
    } else if (seeding_algorithm >= 1000) {
        // This has been added for debugging; allows to reproduce the sequence of the random numbers used in the run
        g = random_engine(seeding_algorithm * (me + 1));
        std::cout << " Seed for process # " << me  << " is " << seeding_algorithm * (me + 1) << std::endl;
    } else {
        int seed = chrono::system_clock::now().time_since_epoch().count() * (me + 1);
        g = random_engine(seed);
        std::cout << " Seed for process # " << me  << " is " << seed << std::endl;
    }
    // Populate initial walker ensemble
    bool uniform_init = (p["fciqmc_guess_subspace"] > 0 ? false : true);
    // This will be used in initialize if uniform_init == false
    if(!uniform_init) init_guess_subspace = std::min(p["fciqmc_guess_subspace"], int(gb.get_basis_size())); 
    initialize(uniform_init);

    if (me == 0) {
        std::vector<double> ndet_per_proc;
        ndet_per_proc.resize(size, 0.0);
        //std::cout << " Determinant mapping " << std::endl;
        for (size_t b = 0; b < gb.get_basis_size(); b++) {
        //    std::cout << " Det # " << b << " => " fnv_hash(b) << std::endl;
            ndet_per_proc[fnv_hash(b)] += 1.0;
        }
        if (debug) {
            std::cout << " Basis function counts per MPI rank " << std::endl;
            for (size_t iproc = 0; iproc < size; iproc++) 
                std::cout << "Hash function assigns " << ndet_per_proc[iproc] << " to rank # " << iproc << std::endl; 
        }
    }
}

/*
void FCIQMC_mpi::update_walker_lists() {
    // The function handles communication
    // of the spawned walkers to the assigned processes;
    // requires local_n_spawned and local_spawned to be 
    // set up correctly

    assert (local_n_spawned.size() == local_spawned.size());
    int total = 0, total_ = 0;
    for (const auto &l : local_spawned) {
        total += l.size();
    }

    total_ = std::accumulate(local_n_spawned.begin(), local_n_spawned.end(), 0);
    assert (total == total_);

    if (debug) {
        for (int iproc = 0; iproc < size; iproc++) {
            if (me != iproc) {
                MPI_Barrier(MPI_COMM_WORLD);
            } else {
                std::cout << " Local spawned on process # " << iproc << " at iteration " << local_it_count << std::endl;
                for (int jproc = 0; jproc < size; jproc++) {
                    std::cout << "To be sent to # " << jproc << std::endl;
                    for (const auto &w : local_spawned[jproc] )
                        std::cout << "(" << w.det_id << ", " << w.weight << ")" << std::endl;
                }
                std::cout.flush();
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }


    for (int i = 0; i < size; i++) {
        std::fill(global_n_spawned.begin(), global_n_spawned.end(), 0);
        std::fill(disp.begin(), disp.end(), 0);
        MPI_Gather(&local_n_spawned[i], 1, MPI_INT, global_n_spawned.data(), 1, MPI_INT, i, MPI_COMM_WORLD);
        if (i == me) {
            // Sum up all the elements of local_n_spawned and allocate 
            // an appropriate amount of memory
            int total_n_spawned_on_me = std::accumulate(global_n_spawned.begin(), global_n_spawned.end(), 0);
            if(debug) std::cout << " Total spawned on # " << i << " is " << total_n_spawned_on_me << " at iteration " << local_it_count << std::endl;
            global_spawned.resize(total_n_spawned_on_me);
            for (size_t iproc = 0; iproc < size - 1; iproc++) disp[iproc + 1] = disp[iproc] + global_n_spawned[iproc];
            if (debug) {
                std::cout << "Process " << i << " recieved the following walker counts " << " at iteration " << local_it_count << std::endl;
                for (int iproc = 0; iproc < size; iproc++) 
                    std::cout << " From rank # " << iproc << ": " << global_n_spawned[iproc] <<std::endl;
            }

            std::cout.flush();
        }
        MPI_Barrier(MPI_COMM_WORLD); // Do I need it here?
        MPI_Gatherv(local_spawned[i].data(), (int)local_spawned[i].size(), Walker_dt, global_spawned.data(), global_n_spawned.data(), disp.data(), Walker_dt, i, MPI_COMM_WORLD);
        // Merge the walkers to the main array
        if (i == me) {
            for (const auto &w : global_spawned) {
                if (m_walker_ensemble.find(w.det_id) != m_walker_ensemble.end()) {
                    m_walker_ensemble[w.det_id] += w.weight;
                } else {
                    m_walker_ensemble[w.det_id] = w.weight;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // This is needed in order to be able to calculate the local walker counts correctly

}
*/

void FCIQMC_mpi::update_walker_lists() {
    // The function handles communication
    // of the spawned walkers to the assigned processes;
    // requires local_n_spawned and local_spawned to be 
    // set up correctly

    assert (local_n_spawned.size() == local_spawned.size());
    int total = 0, total_ = 0;
    for (const auto &l : local_spawned) {
        total += l.size();
    }

    total_ = std::accumulate(local_n_spawned.begin(), local_n_spawned.end(), 0);
    assert (total == total_);

    // Obtain walker counts for me in a single AlltoAll call
    std::fill(global_n_spawned.begin(), global_n_spawned.end(), 0);
    std::fill(sdisp.begin(), sdisp.end(), 0);
    std::fill(rdisp.begin(), rdisp.end(), 0);
    MPI_Alltoall(local_n_spawned.data(), 1, MPI_INT, global_n_spawned.data(), 1, MPI_INT, MPI_COMM_WORLD);
    // Determine how many walkers the current process receives; resize global_spawned and calculate displacements
    int total_n_spawned_on_me = std::accumulate(global_n_spawned.begin(), global_n_spawned.end(), 0);
    global_spawned.resize(total_n_spawned_on_me);
    for (size_t iproc = 0; iproc < size - 1; iproc++) {
        rdisp[iproc + 1] = rdisp[iproc] + global_n_spawned[iproc];
        sdisp[iproc + 1] = sdisp[iproc] + local_n_spawned[iproc];
    }
    // Resize local_spawned and communicate walkers using MPI_Alltoallv
    std::vector<Walker> t_local_spawned;
    for (size_t iproc = 0; iproc < size; iproc++) 
        t_local_spawned.insert(t_local_spawned.end(), local_spawned[iproc].begin(), local_spawned[iproc].end());
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoallv(t_local_spawned.data(), local_n_spawned.data(), sdisp.data(), Walker_dt, global_spawned.data(), global_n_spawned.data(), rdisp.data(), Walker_dt, MPI_COMM_WORLD);
    // Merge all the walkers assigned to me
    for (const auto &w : global_spawned) {
        if (m_walker_ensemble.find(w.det_id) != m_walker_ensemble.end()) {
            m_walker_ensemble[w.det_id] += w.weight;
        } else {
            m_walker_ensemble[w.det_id] = w.weight;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // This is needed in order to be able to calculate the local walker counts correctly

}


int FCIQMC_mpi::fnv_hash(const int &src) {

    int hash = 0;
    constexpr size_t p = 1099511628211; // Mol. Phys. 112, 1855 (2014)
    auto [ na, nb ] = gb.get_ab();
    auto [ ia, ib ] = gb.unpack_str_index((size_t)src);
    auto src_a = gb.a(ia);
    for (size_t i = 0; i < na; i++) hash += (p*hash + i * src_a[i]);
    if (nb > 0) {
        auto src_b = gb.b(ib);
        for (size_t i = 0; i < nb; i++) hash += (p*hash + (i + na) * 2 * src_b[i]); // This needs to be updated so that indeces of beta orbitals are offset correctly!
    }

    return abs(hash % size);

}

void FCIQMC_mpi::initialize (bool uniform) {

    auto basis_size = gb.get_basis_size();
    double max_diag = std::numeric_limits<double>::min(), min_diag = std::numeric_limits<double>::max(); 
    double max_diag_global, min_diag_global;

    if (uniform) {
        assert (local_spawned.size() == size && local_n_spawned.size() == size);
	uniform_int_distribution<int> init_distr(0, basis_size - 1);
        for (int n = 0; n < m_N; n++) {
            int idx = init_distr(g);
            int proc = fnv_hash(idx);
            Walker new_walker {idx, 1};
            merge_walker(new_walker, local_spawned[proc]);
            assert (std::is_sorted(local_spawned[proc].begin(), local_spawned[proc].end(), comp_less));
	}

        for (size_t iproc = 0; iproc < size; iproc++)
            local_n_spawned[iproc] = local_spawned[iproc].size();

        MPI_Barrier(MPI_COMM_WORLD);
        update_walker_lists();
        
	for (int i = 0; i < basis_size; i++) {
            if (me == fnv_hash(i)) {
                double diag_i = gh.matrix(i, i);
                max_diag = max(max_diag, diag_i);
                min_diag = min(min_diag, diag_i);
            }
	}

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&max_diag, &max_diag_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&min_diag, &min_diag_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	m_E_T = min_diag_global, m_E_M = min_diag_global;  // initial guess for the ground state energy
	m_N = 0; m_N_uniq = 0;
	for (const auto &w : m_walker_ensemble) {
            m_N += abs(w.second); 
            m_N_uniq++;
        }

        max_diag = max_diag_global;
        min_diag = min_diag_global;

    } else {
        //assert (size == 1);
        assert (local_spawned.size() == size && local_n_spawned.size() == size);
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
        for (int i = 0; i < m_N; i++) {
            int tr_det = init_distr2(g);
            //std::cout << tr_det << std::endl;
            int sign = guess_wfn[tr_det] > 0 ? 1 : -1;
            int proc = fnv_hash(tr_gb.get_id(tr_det));
            Walker new_walker {(int)tr_gb.get_id(tr_det), sign};
            merge_walker(new_walker, local_spawned[proc]);
            assert (std::is_sorted(local_spawned[proc].begin(), local_spawned[proc].end(), comp_less));
        }

        for (size_t iproc = 0; iproc < size; iproc++) 
            local_n_spawned[iproc] = local_spawned[iproc].size();

        MPI_Barrier(MPI_COMM_WORLD);
        update_walker_lists();
        m_N = 0; // We need to redefine the local number of walkers on the current process
        m_N_uniq = 0;
        for (const auto &w : m_walker_ensemble) {
            m_N += abs(w.second);
            m_N_uniq++;
        }

        m_E_T = guess_en[0]; m_E_M = guess_en[0];
        max_diag = *std::max_element(H_diag.begin(), H_diag.end());
        min_diag = *std::min_element(H_diag.begin(), H_diag.end());
    }

    if (dt < 0) 
	dt = 0.1 *  1. / ( max_diag - min_diag); // See PRB 44 9410 (1991)

    // Collect the total number of walkers from all the ranks
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&m_N, &m_N_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&m_N_uniq, &m_N_uniq_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (me == 0) {
        std::cout << " Maximum H_ii = " << max_diag << std::endl;
        std::cout << " Minimum H_ii = " << min_diag << std::endl;
        std::cout << " Number of walkers = " << m_N_global << std::endl;
        std::cout << " Number of unique walkers = " << m_N_uniq_global << std::endl;
        printf(" Imaginary time-step = %20.10f\n", dt);
        printf(" Initial guess for the energy offset = %10.6f\n", m_E_T);

        // This is added temporary for debugging purposes
        // Printing the initial walker ensemble on the master process
        /*
        std::vector<int> extended_ensemble (gb.get_basis_size(), 0);
        for (const auto &w : m_walker_ensemble) extended_ensemble[w.first] += w.second;
        std::cout << " Initial walker ensemble for reference purposes " << std::endl;
        for (size_t jb = 0; jb < gb.get_basis_size(); jb++) std::cout << " On b.f. # " << jb << " : " << extended_ensemble[jb] << std::endl;
        std::cout << " Control random number (to test the state of generator) " << g() << std::endl;
        */
    }
}

// High level driver to start the simulation

void FCIQMC_mpi::run() {

	// Structure
	// 1. Equilibration
	// 2. Production run:
	// 2.1 Loop over steps
	// 2.1.1 Update E_T as soon as the block ends
	// 2.1.2 Collect stats
	
	// Some parameters
	if (B < 0)
            B = 1.0; // default damping parameter for population control

        double t_equil, t_prod;
        bool equil_mode = true;

        // Running statistics (GSL)
        gsl_rstat_workspace *rstat_m, *rstat_g;
        if (me == 0) {
            rstat_m = gsl_rstat_alloc(); rstat_g = gsl_rstat_alloc();
            printf( ">>>>>>>>> Running FCIQMC calculation on %d processes <<<<<<<<\n", size);

            // Starting equilibration run here (dt will not be adjusted)
            printf( "block #  total pop.  E_m (mixed)  E_g (growth)  <E_m>   <E_g>\n");
            printf( "=======  ==========  ===========  ============  =====   =====\n");
            printf( "---------------------- Equilibration run --------------------\n");
            std::cout.flush();
        }

        t_equil = -MPI_Wtime();

        for (size_t iblock = 0; iblock < m_N_equil; iblock++) {
	    int N_after, N_before = m_N_global; // Global has to be up-to-date
            for (size_t istep = 0; istep < m_steps_per_block; istep++) {
                local_it_count++;
                run_step(debug, equil_mode);
                MPI_Barrier(MPI_COMM_WORLD);
                update_walker_lists();
                m_N = 0; m_N_uniq = 0;
                for (const auto &w : m_walker_ensemble) {
                    m_N += abs(w.second); 
                    m_N_uniq++;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Allreduce(&m_N, &m_N_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Allreduce(&m_N_uniq, &m_N_uniq_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
                if (m_N_global == 0) {
                    std::cout << "All walkers died! Aborting..." << std::endl; 
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            N_after = m_N_global; // !!
            m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before); 
            if (me == 0) {
                printf( "%-7d %-10d %-13.6f %-13.6f %-13s %-13s\n", iblock, get_num_total(), m_E_M, m_E_T, "-", "-"); 
                std::cout.flush();
            }
        }

        t_equil += MPI_Wtime();

        //MPI_Barrier(MPI_COMM_WORLD);

	// Production run

        if (me == 0) {
            printf( "---------------------- Production run ------------------------\n");
            std::cout.flush();
        }

        t_prod = -MPI_Wtime();
        equil_mode = false; // We will now record diagonal elements of the Hamiltonian

        for (size_t iblock = 0; iblock < m_N_blocks; iblock++) {
	    int N_after, N_before = m_N_global;
            for (size_t istep = 0; istep < m_steps_per_block; istep++) {
                local_it_count++;
                run_step(debug, equil_mode);
                MPI_Barrier(MPI_COMM_WORLD);
                update_walker_lists();
                m_N = 0; m_N_uniq = 0;
                for (const auto &w : m_walker_ensemble) {
                    m_N += abs(w.second); 
                    m_N_uniq++;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Allreduce(&m_N, &m_N_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Allreduce(&m_N_uniq, &m_N_uniq_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
                if (m_N_global == 0) {
                    std::cout << "All walkers died! Aborting..." << std::endl; 
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            N_after = m_N_global;
            m_E_T -= B / (m_steps_per_block * dt) * log (double(N_after) / N_before); 
            // Projected estimator will be calculated below
            double e_num = 0.0, e_denom = 0.0;
            double e_num_global = 0.0, e_denom_global = 0.0;
            for (const auto &w : m_walker_ensemble) {
                if (w.second == 0) continue;
                if (mb_ints.find(w.first) == mb_ints.end()) {
                    auto [n_, d_] = en_proj.eval(w.first);
                    mb_ints[w.first] = std::make_tuple(n_, d_);
                }
                auto & [num, denom] = mb_ints[w.first];
                e_num += num * w.second;
                e_denom += denom * w.second;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&e_num, &e_num_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&e_denom, &e_denom_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            // Only root process will have the correct m_E_M 
            if (me == 0) {
                assert (abs(e_denom_global) >= 1e-10);
                m_E_M = e_num_global / e_denom_global;
                gsl_rstat_add(m_E_M, rstat_m);
                gsl_rstat_add(m_E_T, rstat_g);
                printf( "%-7d %-10d %-13.6f %-13.6f %-13.6f %-13.6f\n", iblock, get_num_total(), m_E_M, m_E_T, gsl_rstat_mean(rstat_m), gsl_rstat_mean(rstat_g)); 
                std::cout.flush();
            }
        }

        t_prod += MPI_Wtime();

        if (me == 0) {
            gsl_rstat_free(rstat_m);
            gsl_rstat_free(rstat_g);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double t_prod_, t_equil_; // recv buffers
        MPI_Reduce(&t_equil, &t_equil_, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_prod, &t_prod_, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (me == 0) {
            printf("(FCIQMC_mpi) Maximum equilibration run time %20.2f s\n", t_equil_);
            printf("(FCIQMC_mpi) Maximum production run time %20.2f s\n", t_prod_);
        }
        
}

void FCIQMC_mpi::run_step(bool verbose, bool equil) {

    int basis_size = int(gb.get_basis_size());
    const int N_0 = m_N; // Local walker count; N_pr should coincide with this number at the end of the loop
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<int> disp_walker(1, basis_size - 1); // This should be shared among threads
    int total_spawned = 0, anti_creations = 0, total_killed = 0, total_cloned = 0, N_pr = 0;
    //std::cout << " Control rn (from mt) : " << g() << std::endl;
    //for (size_t dummy = 0; dummy < 1000; dummy++) std::cout << " Control rn : " << disp_walker(g) << std::endl;
    //for (size_t dummy = 0; dummy < 1000; dummy++) std::cout << " Control rn : " << u(g) << std::endl;

    //const random_engine &r_g = g; // This is not allowed :(

    std::fill(local_n_spawned.begin(), local_n_spawned.end(), 0); 
    // For added safety?
    for (size_t pid = 0; pid < size; pid++) {
        local_spawned[pid].resize(0);
        assert (local_spawned[pid].size() == 0);
    }

    for (auto &w : m_walker_ensemble) {
        //const size_t i = w.first;
        const int i = w.first;
        const int n_walkers = abs(w.second);
        const int sign_ref = ( w.second > 0 ? 1 : -1);
        if (n_walkers == 0) continue; // This is important since we are working with a full population vector
        //std::cout << " On det # " << i << std::endl;
        for (size_t iw = 0; iw < n_walkers; iw++) {
            N_pr++;
            auto di = disp_walker(g);
            //std::cout << di << std::endl;
            size_t j = size_t ((i + di) % basis_size);
            double h = gh.matrix(i, j);
            int sign = (h >= 0 ? 1 : -1) * (-1) * sign_ref;
            double ps = abs(h) * dt * (basis_size - 1); 
            int survivors = int(ps);
            if (ps - survivors > u(g)) survivors++;
            if (verbose) std::cout << "Spawned " << survivors << " on det # " << j << std::endl;
            int proc = fnv_hash(j);
            total_spawned += survivors;
            Walker new_walkers {(int)j, sign * survivors};
            merge_walker(new_walkers, local_spawned[proc]);
            assert (std::is_sorted(local_spawned[proc].begin(), local_spawned[proc].end(), comp_less));
        }
        //std::cout << "Finished spawning for walker # " << i << endl;
        // Birth-death process
        double h_ii;
        if (equil) {
            h_ii = gh.matrix(i, i);
        } else {
            // Attemp to find the element in the lookup table; calculate and store if it's not there
            if (h_diag.find(i) == h_diag.end()) h_diag[i] = gh.matrix(i, i); 
            h_ii = h_diag[i];
        }
        double rate = dt * (h_ii - m_E_T) * n_walkers;
        int nkill = int(abs(rate));
        if((abs(rate) - nkill) > u(g)) nkill++;

        if (abs(rate) / n_walkers  > 2) {
            std::cout << "Death rate per walker is " << abs(rate) / n_walkers << " on rank # " << me << std::endl;
            std::cout << "Consider decreasing the time step! " << std::endl;
        }

        if (rate >= 0) {
            int new_unsigned_weight = n_walkers - nkill;
            if (new_unsigned_weight < 0) {
                anti_creations += abs(new_unsigned_weight);
                total_killed += n_walkers;
            } else {
                total_killed += nkill;
            }
            w.second = new_unsigned_weight * sign_ref;
        } else {
            //std::cout << "Inside cloning branch" << endl;
            total_cloned += nkill;
            w.second = (n_walkers + nkill) * sign_ref;
        }
    }

    // Since the walker counts cannot be defined during spawning process itself 
    // we have to calculate them afterward based on the content of the local_spawned
    for (size_t iproc = 0; iproc < size; iproc++) 
        local_n_spawned[iproc] = local_spawned[iproc].size();

    // Perform some obvious checks
    assert (N_pr == N_0);
    // -----------------------------------------------------------------------------------
    // The following code block is added for debugging purposes and will be commented out 
    // in the final version of the code
    /*
    std::vector<int> spawned_combined(gb.get_basis_size(), 0);
    for (int iproc = 0; iproc < size; iproc++) 
        for (const auto &w : local_spawned[iproc])
            spawned_combined[w.det_id] += w.weight;

    std::cout << " >> Spawned walker array (only non-zero elements) << " << std::endl;
    for (size_t jb = 0; jb < gb.get_basis_size(); jb++) 
        if (spawned_combined[jb] != 0) std::cout << " On b.f. # " << jb << " : " << spawned_combined[jb] << std::endl;

    std::vector<int> extended_ensemble (gb.get_basis_size(), 0);
    for (const auto &w : m_walker_ensemble) extended_ensemble[w.first] += w.second;

    std::cout << " >> Main walker array (only non-zero elements) << " << std::endl;
    for (size_t ib = 0; ib < gb.get_basis_size() ; ib++)
        if (extended_ensemble[ib] != 0) std::cout << " On b.f. # " << ib << " : " << extended_ensemble[ib] << std::endl; 
    */
    // -----------------------------------------------------------------------------------

    if (verbose) {

        std::cout << "At iteration " << local_it_count << std::endl;
        std::cout << "Number of walkers processed in the loop = " << N_pr << " on process # " << me << std::endl;
	std::cout << "Number of anti-particles created = " << anti_creations << " on process # " << me << std::endl;
        std::cout << "Number of particles killed = " << total_killed << " on process # " << me << std::endl;
	std::cout << "Death rate = " << double(total_killed)/N_0 << " on process # " << me << std::endl;
	std::cout << "Number of particles cloned = " << total_cloned << " on process # " << me << std::endl;
	std::cout << "Number of particles spawned = " << total_spawned << " on process # " << me << std::endl;
	std::cout << "Average spawning rate = " << double(total_spawned)/N_0 << " on process # " << me << std::endl;

        if (double(total_spawned)/N_0 > 1) {
            std::cout << "WARNING: the spawning rate is = " << double(total_spawned)/N_0 << " on process # " << me
                      << std::endl << "consider reducing the time step..." << std::endl;
        }
    }
}

void FCIQMC_mpi::save_walkers(fstream &out) {
    // This should probably be serialized (i.e. one process collects
    // all the walkers and saves them to a file)
    assert (out.is_open());
    auto n_bf = gb.get_basis_size();
    for (int i = 0; i < n_bf; i++) {
        out << m_walker_ensemble[i] << std::endl;;
    }
}

FCIQMC_mpi::~FCIQMC_mpi() {
    // No need to do anything :)
}
