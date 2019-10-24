#ifndef QESTIMATOR_H
#define QESTIMATOR_H

#include "qhamiltonian.h"
#include "qintegral.h"
#include <vector>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include "qorbitals.h"
#include "qparams_reader.h"

class Estimator {
    protected:
	Integral_factory &ig;
        Basis &bas_full;
	Hamiltonian h_full;
    public:
        Estimator(Integral_factory &ig_, Basis &bas_full_) : ig(ig_), bas_full(bas_full_), h_full(ig, bas_full) {} 
        virtual std::tuple<double, double> eval (std::vector<double> &wf) = 0;
        virtual std::tuple<double, double> eval (size_t idet) = 0;
};

class MixedBasisEstimator : public Estimator {
	private:
            ShellSet ss;
            Becke_grid g;
            TruncatedBasis *aux_bas_tr; 
            DetBasis *aux_bas_full;
            Aux_integrals *aux_int;
            std::vector<double> overlap;
            std::vector<double> trial_state;
            double trial_e;

            double calc_orb_overlap(size_t i, size_t j); // Evaluates orbital overlap

            double calc_overlap1(size_t i_aux, size_t j_full); // Evaluates determinant overlap for a number of special cases
            double calc_overlap2(size_t i_aux, size_t j_full);
            double calc_overlap(size_t i_aux, size_t j_full); // General function for Ne >= 3
            bool test_eval2(); // Implements a simple test of the eval function for 2e

        public:
            MixedBasisEstimator(Params_reader &q, Integral_factory &int_f, Basis &bas);
            ~MixedBasisEstimator();
            std::tuple<double, double> eval (std::vector<double> &wf);
            std::tuple<double, double> eval (size_t idet);
};

class ProjEstimator : public Estimator {
	private:
		TruncatedBasis *tr_basis; 
		std::vector<double> trial_state;
		double trial_e;

	public:
		ProjEstimator(std::map<string, int> &p, Integral_factory &int_f, Basis &bas) : Estimator(int_f, bas) {
                    auto d = h_full.build_diagonal();
                    size_t subspace_size = std::min(bas_full.get_basis_size(), size_t(p["fciqmc_projection_subspace"]));
                    if (subspace_size <= 0) subspace_size = 1;
                    tr_basis = new TruncatedBasis(p, ig.n1porb, subspace_size, d, bas_full);
                    trial_state.resize(tr_basis->get_basis_size());
                    Hamiltonian h_proj(ig, *tr_basis);
                    {
                        auto e = h_proj.diag(true); // Saving the ground state
			trial_e = e[0];
			auto v = h_proj.get_wfn();
			std::copy(v.begin(), v.end(), trial_state.begin());
                    }
                    printf("Energy of the trial function inside ProjEstimator is %13.6f\n", trial_e);
		    printf("Basis size inside ProjEstimator is %d\n", tr_basis->get_basis_size());
		}
		std::tuple<double, double>  eval(std::vector<double> &wf) {
		   	// Full w.f. version
			double num = 0.0, denom = 0.0;
			auto n_bf = tr_basis->get_basis_size();
			//auto n_bf = tr_basis.get_full_basis().get_basis_size();
			for (size_t i = 0; i < n_bf; i++) {
				auto &connected = tr_basis->get_neigh(i);
				//auto &connected = tr_basis.get_full_basis().get_neigh(i);
				auto id = tr_basis->get_id(i);
				//auto id = i;
				denom += trial_state[i] * wf[id];
				for (const auto &j : connected) {
					num += trial_state[i] * h_full.matrix(id, j) * wf[j];
				}
				/*
				for (size_t j = 0; j < n_bf; j++) {
					num += trial_state[i] * h_full.matrix(i, j) * wf[j];
				}
				*/
			}

			return std::make_tuple(num, denom);
		}
		std::tuple<double, double> eval(size_t idet) { 
			// Single determinant version
			// idet should be a valid determinant id
			// with respect to the full basis
			double num = 0.0, denom = 0.0;
			auto n_bf = tr_basis->get_basis_size();
			for (size_t i = 0; i < n_bf; i++) {
				// loop over small basis
				auto &connected = tr_basis->get_neigh(i);
				auto id = tr_basis->get_id(i);
				if (id == idet) denom += trial_state[i];
				if (std::binary_search(connected.begin(), connected.end(), idet)) {
					num += trial_state[i] * h_full.matrix(id, idet);
				}
			}

			return std::make_tuple(num, denom);
		}
                ~ProjEstimator() {
                    delete tr_basis;
                }
			
};


#endif
