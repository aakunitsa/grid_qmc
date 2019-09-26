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
		ProjEstimator(Integral_factory &int_f, TruncatedBasis &tr_b) : tr_basis(tr_b), ig(int_f), h_full(ig, tr_basis.get_full_basis()) {
			trial_state.resize(tr_basis.get_basis_size());
			printf("Iside ProjEstimator!\n" );
			Hamiltonian h_proj(ig, tr_basis);
			{
				auto e = h_proj.diag(true); // Saving the ground state
				trial_e = e[0];
				auto v = h_proj.get_wfn();
				std::copy(v.begin(), v.end(), trial_state.begin());
			}

			printf("Energy of the trial function is %13.6f\n", trial_e);
			printf("Basis size inside ProjEstimator is %d\n", tr_basis.get_basis_size());

		}
		std::tuple<double, double>  eval(std::vector<double> &wf) {
		   	// Full w.f. version
			double num = 0.0, denom = 0.0;
			auto n_bf = tr_basis.get_basis_size();
			for (size_t i = 0; i < n_bf; i++) {
				auto &connected = tr_basis.get_neigh(i);
				auto id = tr_basis.get_id(i);
				denom += trial_state[i] * wf[id];
				for (const auto &j : connected) {
					num += trial_state[i] * h_full.matrix(id, j) * wf[j];
				}
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
