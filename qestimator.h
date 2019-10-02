#ifndef QESTIMATOR_H
#define QESTIMATOR_H

#include "qhamiltonian.h"
#include "qintegral.h"
#include <vector>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <cstdio>

class ProjEstimator {
	private:
		TruncatedBasis &tr_basis;
		Integral_factory &ig;
		Hamiltonian h_full;

		std::vector<double> trial_state;
		double trial_e;

	public:
		ProjEstimator(Integral_factory &int_f, TruncatedBasis &tr_b) : tr_basis(tr_b), ig(int_f), h_full(ig, tr_b.get_full_basis()) {
			trial_state.resize(tr_basis.get_basis_size());
			//trial_state.resize(tr_b.get_full_basis().get_basis_size());
			//printf("Iside ProjEstimator!\n" );
			//printf("Truncated basis size is %d \n", tr_b.get_basis_size());
			Hamiltonian h_proj(ig, tr_basis);
			{
				auto e = h_proj.diag(true); // Saving the ground state
				//auto hd = h_full.build_diagonal();
				//std::cout << "Finished building diagonal" << std::endl;
				//auto e = h_full.diag(true); // Saving the ground state
				trial_e = e[0];
				auto v = h_proj.get_wfn();
				//auto v = h_full.get_wfn();
				std::copy(v.begin(), v.end(), trial_state.begin());
			}

			//printf("Energy of the trial function inside ProjEstimator is %13.6f\n", trial_e);
			//cout.flush();
			//double rel_e = h_full.check_wfn();
			//assert (abs(rel_e - trial_e) <= 1e-10);
			//printf("Basis size inside ProjEstimator is %d\n", tr_basis.get_basis_size());

		}
		std::tuple<double, double>  eval(std::vector<double> &wf) {
		   	// Full w.f. version
			double num = 0.0, denom = 0.0;
			auto n_bf = tr_basis.get_basis_size();
			//auto n_bf = tr_basis.get_full_basis().get_basis_size();
			for (size_t i = 0; i < n_bf; i++) {
				auto &connected = tr_basis.get_neigh(i);
				//auto &connected = tr_basis.get_full_basis().get_neigh(i);
				auto id = tr_basis.get_id(i);
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
			auto n_bf = tr_basis.get_basis_size();
			for (size_t i = 0; i < n_bf; i++) {
				// loop over small basis
				auto &connected = tr_basis.get_neigh(i);
				auto id = tr_basis.get_id(i);
				if (id == idet) denom += trial_state[i];
				if (std::binary_search(connected.begin(), connected.end(), idet)) {
					num += trial_state[i] * h_full.matrix(id, idet);
				}
			}

			return std::make_tuple(num, denom);
		}
			
};


#endif
