#include "qintegral.h"
#include <limits>
#include <list>
#include <array>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <utility>
#include <armadillo>
#include <gsl/gsl_math.h>
#include "poisson_fortran.h"
#include <omp.h>


extern "C" void read_orbitals_(double *coeff);


Integral_factory::Integral_factory(std::map<string, int> &p) : g(p), lp(p) {};

bool Integral_factory::check_orthogonality(double orth_thresh = 1e-6) {

	bool ortho = true;
	double max_dev = std::numeric_limits<double>::min();

	for (int i = 0; i < n1porb; i++) {
		for (int j = 0; j < n1porb; j++) {
			double d_prod = 0.0;
			for (int r = 0; r < g.nrad; r++) {
				for (int a = 0; a < g.nang; a++) {
					int g_idx = r * g.nang + a;
					d_prod += 4. * M_PI * g.gridw_r[r] * g.gridw_a[a] * paux_bf[ i * g.nang * g.nrad + g_idx] * paux_bf[ j * g.nang * g.nrad + g_idx ];
				}
			}

			if (i!= j) { 
				ortho = ortho && (abs(d_prod) <= orth_thresh);
				max_dev = std::max(max_dev, abs(d_prod));
			}
			if (i == j) { 
				ortho = ortho && (abs(d_prod - 1.0) <= orth_thresh);
				max_dev = std::max(max_dev, abs(d_prod - 1.0));
			}
		}
	}

	if (!ortho) printf("Maximum deviation from orthogonality %13.6f", max_dev);

	return ortho;

}

std::vector<double> Integral_factory::expand(std::vector<double> &vec) {
    // Returns the vector of expansion coefficients
    assert (vec.size() == g.nrad * g.nang);
    std::vector<double> exp_coef (n1porb, 0.0);
    for (int i = 0; i < n1porb; i++) {
        for (int r = 0; r < g.nrad;r++) {
            for (int a = 0; a < g.nang; a++) {
                int g_idx = r * g.nang + a;
                exp_coef[i] += 4. * M_PI * g.gridw_r[r] * g.gridw_a[a] * paux_bf[ i * g.nang * g.nrad + g_idx] * vec[g_idx];
            }
        }
    }

    return exp_coef;
}

void Integral_factory::fcidump() {

	fstream int_file;
	int_file.open("QFCIDUMP.POLY", std::ios::out);
	assert(int_file.is_open());

	size_t ngrid = g.nrad * g.nang;

	for (size_t i = 0; i < n1porb; i++) {
		for (size_t j = i; j < n1porb; j++) {

			double h = hc(i, j); // Will be symmetrized inside hc function! 

			int_file << std::scientific << std::setprecision(20) << std::setw(28) << h << "    " ;
			int_file << i + 1 << "    " << j + 1 << "    " << 0 << "    " << 0 << std::endl;
			
		}
	}

	// (Annoying) Two electron part 
	std::vector<int> pairs;

	// populate pairs same way as in Polymer
        /*
	for (int i = 0; i < n1porb; i++) {
		for (int j = i + 1; j < n1porb; j++) {
			pairs.push_back(i * n1porb + j);
		}
	}

	for (int j = 0; j < n1porb; j++) {
		pairs.push_back(j * n1porb + j);
	}
        */
	for (int i = 0; i < n1porb; i++) {
            for (int j = i; j < n1porb; j++) {
		pairs.push_back(i * n1porb + j);
            }
	}
	size_t numpairs = n1porb * (n1porb + 1) / 2;

	assert (pairs.size() == numpairs);

	for ( size_t p = 0; p < numpairs; p++) {
		int ij = pairs[p];
		int j = ij % n1porb;
		int i = (ij - j) / n1porb;
		for (size_t q = p; q < numpairs; q++) {
			int kl = pairs[q];
			int l = kl % n1porb;
			int k = (kl - l) / n1porb;

			double eri[8];
			eri[0] =  ce(i, j, k, l); // Density is built based on the first pair of indeces; same as POLYMER
			
			eri[1] = ce(j, i, k, l);
			eri[2] = ce(i, j, l, k);
			eri[3] = ce(j, i, l, k);
			eri[4] = ce(k, l, i, j);
			eri[5] = ce(l, k, i, j);
			eri[6] = ce(k, l, j, i);
			eri[7] = ce(l, k, j, i);

			double assym = 0.0;
			double assym_thresh = 1e-6;

			for (size_t i = 0 ; i < 8; i++) 
				for ( size_t j = 0; j < 8; j++) {
					assym = std::max(assym, std::abs(eri[i] - eri[j]));
				}

			if ( assym > assym_thresh ) {
				std::cout << " Detected significant symmetry breaking in ERI matrix " << std::endl;
				std::cout << " Assym = " << assym << std::endl;
			}
			

			int_file << std::scientific << std::setprecision(20) << std::setw(28) << eri[0] << "    " ; // See the note above
			int_file << i + 1 << "    " << j + 1 << "    " << k + 1 << "    " << l + 1 << std::endl;

		}
	}

	int_file.close();

}

bool Integral_factory::test_encode() {
    bool ok = true;
    size_t npair = n1porb * (n1porb + 1) / 2;
    size_t nquad = npair * (npair + 1) / 2;
    std::map<size_t, size_t> p;
    for (size_t i = 0; i < n1porb; i++) 
        for (size_t j = 0; j < n1porb; j++) 
            for (size_t k = 0; k < n1porb; k++)
                for (size_t l = 0; l < n1porb; l++) {
                    auto eri_idx = encode(i, j, k, l);
                    if (p.find(eri_idx) != p.end()) {
                        p[eri_idx] += 1;
                    } else {
                        p.insert(std::pair(eri_idx, 1));
                    }
                }

    ok = ok && (p.size() == nquad);

    return ok;
}

inline size_t Integral_factory::encode(size_t i, size_t j, size_t k, size_t l) {
    // Multiindex i, j, k, l represents the integral in Mulliken notation;
    // The input indeces will not be altered inside the function

    // Minor and major index definitions should be swapped

    size_t minor1 = (i > j ? j : i), major1 = (i > j ? i : j), 
           minor2 = (k > l ? l : k), major2 = (k > l ? k : l);

    size_t pair1 = (major1 * (major1 + 1) / 2) + minor1,
           pair2 = (major2 * (major2 + 1) / 2) + minor2;

    if (pair1 > pair2) std::swap(pair1, pair2);
    
    return pair1 + (pair2 * (pair2 + 1) / 2);

}


// Aux_integrals class


Aux_integrals::Aux_integrals(Params_reader &pr, ShellSet &aorb) : ss(aorb), Integral_factory(pr.params) {

	// For convenience
	auto p = pr.params;

	// Numerical thresholds
	//double lthresh = 1e-5; // Linear dependence threshold
	double lthresh = 1e-6; // Linear dependence threshold
	double orth_thresh = 1e-6; // Very mild orthogonalisation accuracy threshold
	//
	
	if (pr.orbital_file.size() != 0) ofile = pr.orbital_file;
	
	// Nuclear charge
	Znuc = double(p["Z"]);

	// Attempt to read orbitals if requested; 
	// Upon failure generate auxiliary basis set 
	// and populate paux_bf
	
	bool gen_orb = true;
	if (p["read_orb_file"] == 1) gen_orb = !read_orbs(); // Generate if reading was not successful

	if (gen_orb) {

		gen_aux_basis(lthresh, orth_thresh);
		n1porb = ss.size() * naux;
		paux_bf.resize(n1porb * g.nrad * g.nang);
		// Populate paux_bf using real spherical harmonics
		// This is a naive and inefficient approach; should probably
		// just represent all the Ylm-s on the Lebedev grid first and then 
		// combine those with aux_bf
	

		for (int ridx = 0; ridx < naux; ridx++) {
			for (int oidx = 0; oidx < ss.size(); oidx++) {
				auto &o = ss.aorb[oidx];
				for (int r = 0; r < g.nrad; r++) {
					for (int a = 0; a < g.nang; a++) {
						auto [th, p] = g.thetaphi_ang[a];
						double Y_real = rY(o.L, o.M, th, p);
						paux_bf[ (ridx * ss.size() + oidx) * g.nang * g.nrad + r * g.nang + a ] = aux_bf[ridx * g.nrad + r] * Y_real;
					}
				}
			}
		}

	}

	bool orth = check_orthogonality(orth_thresh);

	assert (orth);

	// Initialize Poisson solver before starting the calculation
	
	int iatom = int(Znuc), nrad = int(g.nrad) , nang = int(g.nang);
        max_cache_size = std::max(0, p["max_cache_size"]); // Integral cache
	initialize_poisson_(&nrad, &nang, &iatom);
        //bool encode_ok = test_encode();
        //assert(encode_ok);

}

bool Aux_integrals::read_orbs() {

	// Assumes that real valued orbitals are 
	// stored on disk ( the filename is supplied from cmd)
	// The fist three integers on the first line 
	// of the orbital file denonte:
	//
	// number of orbitals
	// maxngrid
	// number of atoms ( currently the code can only work with a single atom
	//
	// The orbitals are stored in one column
	

	// For now this will only properly read text files!!!!!!!!!!!
	// See the fortran subroutine for the binary reading mode definitions
	

	bool binary_mode = false;
	int dot = ofile.find_last_of(".");
	if (dot != string::npos) {
		auto suffix = ofile.substr(dot, ofile.find_first_not_of(" ") + 1);
		if (suffix.compare("bin") == 0 || suffix.compare("BIN") == 0) binary_mode = true;
	}  

	fstream orb_file;
	orb_file.open(ofile.c_str());

	if (!orb_file.is_open()) return false;

	int pmaxngrid, pnatom, porb;
	orb_file >> porb >> pmaxngrid >> pnatom;

	// porb has been set above
	// this needs to be refactored

	if ( size_t(pmaxngrid) != g.nang * g.nrad || pnatom != 1 ) {
		std::cout << " The grid used to generate the orbitals is inconsistent with the current grid! " << std::endl;
		std::cout << " The orbitals from ORBITALS.DAT file will not be used! " << std::endl;
	} else {
		std::cout << porb << " orbitals will be read " << std::endl;
		paux_bf.resize(porb * pmaxngrid);
		if (!binary_mode) {
			for (size_t i = 0; i < size_t (porb * pmaxngrid); i++ ) 
				orb_file >> paux_bf[i];
		} else {
			// use fortran subroutine
			std::cout << "Extracting orbitals from the binary file " << std::endl;
			read_orbitals_(paux_bf.data());
		} 
	}

	orb_file.close();

	return true;


	// Check orthogonality with respect to the grid used in the current code
	// Assemble the weights (radial and angular combined)
	// Orthogonality will be checked in the constructor 
    /*	
	std::vector<double> wei(g.nrad * g.nang, 0.0);

	for ( size_t r = 0; r < g.nrad; r++) 
		for ( size_t a = 0; a < g.nang; a++) 
			wei[r * g.nang + a] = 4. * M_PI * g.gridw_a[a] * g.gridw_r[r];

	arma::mat W = arma::diagmat(arma::vec(wei));
	arma::mat bas(paux_bf.data(), g.nrad * g.nang, porb, false);
	arma::mat id = arma::eye(porb, porb);

	arma::mat dp = arma::abs(bas.t() * W * bas - id);
	std::cout << "Maximum deviation from orthogonality is " << std::scientific << dp.max() << std::endl;
	*/

}

void Aux_integrals::gen_aux_basis(double lthresh, double orth_thresh) {

	std::vector<double> aux(g.nrad * g.nrad, 0.0);
	std::vector<double> bexp(g.nrad), bnorm(g.nrad);
	std::copy(g.r.begin(), g.r.end(), bexp.begin());
	double log2 = log(2.);
	std::transform(bexp.begin(), bexp.end(), bexp.begin(), [&](double r) {return log2 / r;});

	// Generate grid representations of the basis vectors
	
	for (size_t i = 0; i < g.nrad; i++) {
		for (size_t j = 0; j < g.nrad; j++) {
			aux[i * g.nrad + j] = exp(-bexp[i] * g.r[j]);
		}
	}

	// Calculate the overlap matrix of the basis functions
	// Will use armadillo classes instead of plain code
	
	arma::mat aux_bas(aux.data(), g.nrad, g.nrad, false); // create matrix without copying anything
	arma::mat W = arma::diagmat(arma::vec(g.gridw_r)); // Weight matrix

	arma::mat S = aux_bas.t() * W * aux_bas;
	arma::vec es;
	arma::mat ev;

	bool status = arma::eig_sym(es, ev, S);
	assert ( status );

	std::cout << " Maximum eigenvalue of the overlap matrix is " << std::scientific << es.max() << std::endl;
	std::cout << " Minimum eigenvalue of the overlap matrix is " << std::scientific << es.min() << std::endl;

	//ev.print("Eigenvectors of the original overlap matrix: ");
	//es.print("Eigenvalues : ");

	arma::uvec ps = arma::sort_index(es, "ascend");

	// Determine how many vectors should be discarded based on the
	// current lthresh
	
	size_t discard = 0;

	for (size_t k = 0; k < g.nrad; k++) {
		if ( abs(es[ps[k]]) >= lthresh && es[ps[k]] > 0.0 ) {
			discard = k;
			break;
		}
	}

	arma::vec es_sorted(g.nrad - discard);
	arma::mat ev_sorted(g.nrad, g.nrad - discard);

	for ( size_t i = discard; i < g.nrad; i++ ) {
		es_sorted[i - discard] = es[ps[i]];
		std::copy(ev.begin_col(ps[i]), ev.end_col(ps[i]), ev_sorted.begin_col(i - discard));
	}

	//ev_sorted.print(" Sorted eigenvectors : ");
	//es_sorted.print(" Sorted eigenvalues ( the ones below the threshold are excluded ) : ");

	arma::mat invsqrt = sqrt(arma::diagmat(1./es_sorted));
	//invsqrt.print(" Inv s^1/2 :");

	// Orthogonalize eigenvectors  (symmetric)
	
	arma::mat orth_aux_bas = aux_bas * ev_sorted * invsqrt;

	// Normally this should be enough; however, for larger grid sizes due to the relatively 
	// large eigenvalues of the overlap matrix - numerical errors tend to accumulate breaking
	// orthogonality; For this reason the basis needs to be Gram-Schmidt re-orthogonalized 
	// _after_ symmetric orthogonalisation; This is very important! Normally we would use QR 
	// decomposition but since the dot product definition involves radial grid weights - the
	// orthogonalization should be performed manually
	
	arma::mat q(orth_aux_bas.n_rows, orth_aux_bas.n_cols);

	q.col(0) = orth_aux_bas.col(0) / sqrt(arma::as_scalar(orth_aux_bas.col(0).t() * W * orth_aux_bas.col(0)));
	for (size_t ic = 1; ic < orth_aux_bas.n_cols; ic++) {
		arma::vec u = orth_aux_bas.col(ic);
		for (size_t prev = 0; prev < ic; prev++) 
			u -= arma::as_scalar(q.col(prev).t() * W * orth_aux_bas.col(ic)) * q.col(prev);
		q.col(ic) = u / sqrt(arma::as_scalar(u.t() * W * u));
	}

	orth_aux_bas = q; // Copy the new orthogonalized basis and work with it from now on

	// Check orthogonlity with respect to weight matrix
	
	arma::mat diff = arma::abs(orth_aux_bas.t() * W * orth_aux_bas - arma::mat(g.nrad - discard, g.nrad - discard, arma::fill::eye)) ;

	bool within_thresh = arma::all(arma::vectorise(diff) <= orth_thresh);

	if (!within_thresh) {
		std::cout << " Warning! Orthogonalisation error is large than the requested threshold! "  << std::endl;
		std::cout << diff.max() << std::endl;
	}

	// If the grid is small -- print new overlap matrix for visual 
	// inspection
	
	if ( g.nrad <= 25 ) {
		arma::mat overlap = orth_aux_bas.t() * W * orth_aux_bas;
		overlap.print(" Overlap matrix for the orthogonal basis: " );
	}

	// Flatten the array of the aux basis vectors to an std::vector;
	// This can be used/exposed without the need of invoking armadillo
	
	naux = orth_aux_bas.n_cols;
	aux_bf.resize(naux * g.nrad);

	// Storage convention: Fortran
	std::copy(orth_aux_bas.begin(), orth_aux_bas.end(), aux_bf.begin());

	// Check
	size_t col_idx = size_t(naux / 2);
	arma::vec test_col(&aux_bf[col_idx * g.nrad], g.nrad, false);
	double max_diff_col = arma::max(arma::abs(test_col - orth_aux_bas.col(col_idx)));
	assert ( max_diff_col < 1e-12);

}


double Aux_integrals::ce(size_t i_, size_t j_, size_t k_, size_t l_) {

    size_t i = i_, j = j_, k = k_, l = l_;
        if (i > j) std::swap(i, j);
        if (k > l) std::swap(k, l);
        if (i > k || (i == k && l < j)) {
                std::swap(i, k);
                std::swap(j, l);
        } // Ensures that the "left" pair is always less than the "right" pair; the major index of
                          // of a pair is the ___first___ index !!
	double m = 0.0; // Matrix element
	size_t ngrid = g.nrad * g.nang;
        auto eri_idx = encode(i, j, k, l);
        bool found = (cached_eri.find(eri_idx) != cached_eri.end());
        if (found) { 
        //if (false) { 
            // Will add this for debugging
            /*
            double thresh = 1e-2;
            double i_ijkl = eri_fortran_(&paux_bf[i * ngrid], &paux_bf[j * ngrid], &paux_bf[k * ngrid], &paux_bf[l * ngrid]); 
            double i_jikl = eri_fortran_(&paux_bf[j * ngrid], &paux_bf[i * ngrid], &paux_bf[k * ngrid], &paux_bf[l * ngrid]); 
            double i_jilk = eri_fortran_(&paux_bf[j * ngrid], &paux_bf[i * ngrid], &paux_bf[l * ngrid], &paux_bf[k * ngrid]); 
            double i_lkji = eri_fortran_(&paux_bf[l * ngrid], &paux_bf[k * ngrid], &paux_bf[j * ngrid], &paux_bf[i * ngrid]); 

            assert (abs(i_ijkl - i_jikl) < thresh && abs(i_jikl - i_jilk) < thresh && abs(i_jilk - i_lkji) < thresh); // This consistently fails...
            */

            //double i_klij = eri_fortran_(&paux_bf[k * ngrid], &paux_bf[l * ngrid], &paux_bf[i * ngrid], &paux_bf[j * ngrid]); 
            //assert (abs(cached_eri[eri_idx] - i_ijkl) < thresh || abs(cached_eri[eri_idx] - i_klij)); // So this would fail as well
            return cached_eri[eri_idx];
        }
/*
	std::vector<double> gridw(ngrid, 0.0);

	for(size_t r= 0; r < g.nrad; r++) 
		for (size_t a = 0; a < g.nang; a++) 
			gridw[r * g.nang + a] = 4. * M_PI * g.gridw_r[r] * g.gridw_a[a];

    arma::vec orb_i(&paux_bf[i * ngrid], ngrid, false);
    arma::vec orb_j(&paux_bf[j * ngrid], ngrid, false);
    arma::vec orb_k(&paux_bf[k * ngrid], ngrid, false);
    arma::vec orb_l(&paux_bf[l * ngrid], ngrid, false);

*/
	double i_ijkl = eri_fortran_(&paux_bf[i * ngrid], &paux_bf[j * ngrid], &paux_bf[k * ngrid], &paux_bf[l * ngrid]); 
#ifdef _OPENMP
        if (omp_in_parallel()) {
            #pragma omp critical
            if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx, i_ijkl));
        } else {
            if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx, i_ijkl));
        }
#else
        if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx, i_ijkl));
#endif

	// Will use Poisson solver to generate integrals
	// This is not practical for larger scale calculations; here I decided to try that for testing 
	// purposes
/*		
	// 1. Generate density for the first orbital pair, say (i, j)

	std::vector<double> den_ij(g.nrad * g.nang, 0.0);
	std::vector<double> pot_ij(g.nrad * g.nang, 0.0);

	for ( size_t g = 0; g < ngrid ; g++) 
		den_ij[g] = orb_i[g] * orb_j[g];
		

	// 2. Use Poisson solver to calculate corresponding potenials

	construct_potential_(den_ij.data(), pot_ij.data());

	// 3. Contract the potentials with the other orbital density;

	for ( size_t g = 0; g < ngrid ; g++) 
		m += orb_k[g] * orb_l[g] * pot_ij[g] * gridw[g];

	if (abs(m - i_ijkl) > 1e-8) {
		std::cout << "Fortran " << std::setprecision(20) << i_ijkl << std::endl;
		std::cout << "C++ " << std::setprecision(20) << m << std::endl;
		assert (abs(m - i_ijkl) <= 1e-8); // Make sure that C++ and Fortran produce similar results!
	}

	return m ;
*/
        return i_ijkl;

}

double Aux_integrals::hc(size_t i, size_t j) {

	double h = 0.0;
	// Assymetric kinetic part
	h += 0.5 * (ke(i, j) + ke(j, i));
	// Symmetric electron nucleus attraction
	h += vn(i, j);

	return h;
}

double Aux_integrals::ke(size_t i, size_t j) {

	// This function may need to be redesigned
	// to use std::complex instead of handling 
	// all the complex arithmetic manually

    size_t ngrid = g.nrad * g.nang;

    double m_re = 0.0, m_im = 0.0; // Matrix element

    arma::vec orb_i(&paux_bf[i * ngrid], ngrid, false);
    arma::vec orb_j(&paux_bf[j * ngrid], ngrid, false);

	std::vector<double> Ri_re(g.nrad * ss.size(), 0.0), 
		                Ri_im(g.nrad * ss.size(), 0.0),
                        Rj_re(g.nrad * ss.size(), 0.0),
		                Rj_im(g.nrad * ss.size(), 0.0);

	// Perform spherical harmonic expansion for both orbitals 
	
	for (size_t o = 0; o < ss.size(); o++) {
		auto &L = ss.aorb[o].L, &M = ss.aorb[o].M;
		for ( size_t r = 0; r < g.nrad ; r++ ) {
			double di_re = 0.0, di_im = 0.0, dj_re = 0.0, dj_im = 0.0;
			for ( size_t a = 0; a < g.nang ; a++ ) {

				auto [th, p] = g.thetaphi_ang[a];
				auto [re, im] = Y(L, M, th, p);

				// Note: we are multiplying by the complex conjugate; orbitals are real

				di_re += 4. * M_PI * g.gridw_a[a] * re * orb_i[r * g.nang + a];
				di_im -= 4. * M_PI * g.gridw_a[a] * im * orb_i[r * g.nang + a];
				dj_re += 4. * M_PI * g.gridw_a[a] * re * orb_j[r * g.nang + a];
				dj_im -= 4. * M_PI * g.gridw_a[a] * im * orb_j[r * g.nang + a];

			}

			Ri_re[o * g.nrad + r] = di_re;
			Ri_im[o * g.nrad + r] = di_im;
			Rj_re[o * g.nrad + r] = dj_re;
			Rj_im[o * g.nrad + r] = dj_im;
		}
	}

	// Calculate the kinetic energy integral
	
	std::vector<double> lapl_re(g.nrad, 0.0), lapl_im(g.nrad, 0.0);
	
	for ( size_t o = 0; o < ss.size(); o++) {
		std::fill(lapl_re.begin(), lapl_re.end(), 0.0);
		std::fill(lapl_im.begin(), lapl_im.end(), 0.0);
		int L = ss.aorb[o].L;
		lp.apply_fortran(&Rj_re[o * g.nrad], lapl_re.data()); 
		lp.apply_fortran(&Rj_im[o * g.nrad], lapl_im.data()); 
		for ( size_t r  = 0; r < g.nrad; r++) {

			double t_re = 0.0, t_im = 0.0;

			t_re = lapl_re[r] - (double(L * (L + 1))/ gsl_pow_2(g.r[r]) * Rj_re[o * g.nrad + r]);
			t_im = lapl_im[r] - (double(L * (L + 1))/ gsl_pow_2(g.r[r]) * Rj_im[o * g.nrad + r]);
			m_re += (g.gridw_r[r] * (Ri_re[o * g.nrad + r] * t_re + Ri_im[o * g.nrad + r] * t_im));
			m_im += (g.gridw_r[r] * (Ri_re[o * g.nrad + r] * t_im - Ri_im[o * g.nrad + r] * t_re));

		}
	}

	assert ( std::abs(m_im) < 1e-8);
	
	return -0.5 * m_re;

}

double Aux_integrals::vn(size_t i, size_t j) {

	double h = 0.0;
	int ngrid = g.nrad * g.nang;

	arma::vec orb_i(&paux_bf[i * ngrid], ngrid, false);
	arma::vec orb_j(&paux_bf[j * ngrid], ngrid, false);

	for ( size_t k = 0; k < g.nrad ; k++) 
		for ( size_t l = 0; l < g.nang; l++) 
			h += (-Znuc / g.r[k] * orb_j[k * g.nang + l] * orb_i[k * g.nang + l] * g.gridw_r[k] *  g.gridw_a[l] * 4. * M_PI);


	return h;

}


Aux_integrals::~Aux_integrals() {

	finalize_poisson_();

}

// Saved integrals class

Saved_integrals::Saved_integrals(Params_reader &pr) : Integral_factory(pr.params) {

	good = read_fcidump(pr.fcidump_file);
	if (!good) {
		std::cout << " WARNING: reading integrals from disk failed! " << std::endl;
	}

}


bool Saved_integrals::read_fcidump(string &int_file) {

	// Reads FCIDUMP from Polymer and creates arrays for Hcore and ERI
	
	int porb;

	fstream poly_dump;
	poly_dump.open(int_file.c_str());

	size_t nrecords =0 ;

	if (!poly_dump.is_open())  {
		std::cout << " FCIDUMP file was not found " << std::endl;
		return false;
	}

	std::cout << " Reading FCIDUMP " << std::endl;
	poly_dump >> porb;
	std::cout << " Number of basis function in the fcidump file is " << porb << std::endl;
	if (porb > 55109 ) std::cout << " WARNING: the number of orbitals may be too large for the hashing scheme!" << std::endl;

	n1porb = porb; // !!!!!!!!!!!!!!!!!!!!!!!!!!

	int i, j, k, l;
	double x;
	poly_dump >> x >> i >> j >> k >> l;

	do {
		nrecords++;
		//std::cout << "Original: " << x << '\t' << i << '\t' << j  << '\t' << k << '\t' << l << std::endl;

		// We care about two cases only:
		// 1. All indeces are non-zero => eri
		// 2. Two indeces are non-zero => core hamiltonian
		// Everything else will be ignored
		if ((i != 0) && (j != 0) && (k != 0) && (l != 0)) {
			if (i > j) std::swap(i, j);
			if (k > l) std::swap(k, l);
			if (i > k || (i == k && l < j)) {
				std::swap(i, k);
				std::swap(j, l);
			} // Ensures that the "left" pair is always less than the "right" pair; the major index of
                          // of a pair is the ___first___ index !!

			//std::cout << " New index : " << i << '\t' << j << '\t' << k << '\t' << l << std::endl;
			//
			// Index pairs follow lexical ordering as checked by the following assert statement

			assert (i <= j && k <= l && i <= k && j <= ( (i == k) ? l : porb));

			i -= 1; j -= 1; k -= 1; l -= 1;

		    std::cout << x << '\t' << i << '\t' << j  << '\t' << k << '\t' << l << std::endl;

			//int p1 = (porb - 1) * i + j - i * (i - 1) / 2,
			//	p2 = (porb - 1) * k + l - k * (k - 1) / 2;
			//int eri_counter = (p1 + 1) * p1 / 2  + p2;
			
			int eri_counter = gsl_pow_3(porb) * i +gsl_pow_2(porb) * j + porb * k + l;
			//std::cout << " ERI counter " << eri_counter << std::endl;

			assert(eri.find(eri_counter) == eri.end());

			eri.insert(std::pair(eri_counter, x)); // CXX 17 rocks!

		} else if ( (i != 0) && (j != 0) && (k == 0) && (l == 0)) {
			if (i > j) std::swap(i, j);
			i -= 1; 
			j -= 1; 
		    //std::cout << x << '\t' << i << '\t' << j << std::endl;
			assert (hcore.find(i * porb + j) == hcore.end());
			hcore.insert(std::pair(i * porb + j, x));
		}

		poly_dump >> x >> i >> j >> k >> l;

	} while (!poly_dump.eof()) ;

	std::cout << "Finished reading FCIDUMP! " << nrecords << " have been processed." <<  std::endl; 

	poly_dump.close();

	return true;

}

double Saved_integrals::hc(size_t i_, size_t j_) {
    
    size_t i = i_, j = j_;

	if (i > j) std::swap(i, j);
	return hcore[i*n1porb + j];

}


double Saved_integrals::ce(size_t i_, size_t j_, size_t k_, size_t l_) {

    size_t i = i_, j = j_, k = k_, l = l_;

	if (i > j) std::swap(i, j);
	if (k > l) std::swap(k, l);
	if (i > k || (i == k && l < j)) {
		std::swap(i, k);
		std::swap(j, l);
	} 

	assert (i <= j && k <= l && i <= k && j <= ( (i == k) ? l : n1porb));

	//int p1 = (porb - 1) * i + j - i * (i - 1) / 2,
	//	p2 = (porb - 1) * k + l - k * (k - 1) / 2;
	//int eri_index = (p1 + 1) * p1 / 2  + p2;
	
	int eri_index = gsl_pow_3(n1porb) * i +gsl_pow_2(n1porb) * j + n1porb * k + l;

	return eri[eri_index];

}

// Grid_integrals class


Grid_integrals::Grid_integrals(std::map<string, int> &p, ShellSet &aorb) : ss(aorb), r12(p), Integral_factory(p) {

	// Numerical thresholds
	double orth_thresh = 1e-6; // Very mild orthogonalisation accuracy threshold

	Znuc = double(p["Z"]);

	n1porb = ss.size() * g.nrad;
	paux_bf.resize(n1porb * g.nrad * g.nang);
	aux_bf.resize(g.nrad * g.nrad);
	std::fill(aux_bf.begin(), aux_bf.end(), 0.0);
	for (int i = 0; i < g.nrad; i++) aux_bf[i * g.nrad + i] = 1./ sqrt(g.gridw_r[i]);

	for (int ridx = 0; ridx < g.nrad; ridx++) {
		for (int oidx = 0; oidx < ss.size(); oidx++) {
			auto &o = ss.aorb[oidx];
			for (int r = 0; r < g.nrad; r++) {
				for (int a = 0; a < g.nang; a++) {
					auto [th, p] = g.thetaphi_ang[a];
					double Y_real = rY(o.L, o.M, th, p);
					paux_bf[ (ridx * ss.size() + oidx) * g.nang * g.nrad + r * g.nang + a ] = aux_bf[ridx * g.nrad + r] * Y_real;
				}
			}
		}
	}

	bool orth = check_orthogonality(orth_thresh);
	assert (orth);
        max_cache_size = std::max(0, p["max_cache_size"]); // Integral cache
}

double Grid_integrals::ce(size_t i, size_t j, size_t k, size_t l) {

	// Assumes that the orbital indeces are arranged in Mulliken's order, that is
	// i and j refer to the first electron whereas k and l refer to the second
	// Extract radial point and angular orbital index
	
        size_t i_ = i, j_ = j, k_ = k, l_ = l;
        if (i_ > j_) std::swap(i_, j_);
        if (k_ > l_) std::swap(k_, l_);
        if (i_ > k_ || (i_ == k_ && l_ < j_)) {
                std::swap(i_, k_);
                std::swap(j_, l_);
        } // Ensures that the "left" pair is always less than the "right" pair; the major index of
                          // of a pair is the ___first___ index !!
        auto eri_idx_cached = encode(i_, j_, k_, l_);
        bool found = (cached_eri.find(eri_idx_cached) != cached_eri.end());
        if (found) { 
            return cached_eri[eri_idx_cached];
        }

	auto [ir, iorb] = unpack_orb_index(i);
	auto [jr, jorb] = unpack_orb_index(j);
	auto [kr, korb] = unpack_orb_index(k);
	auto [lr, lorb] = unpack_orb_index(l);

	if ( kr != lr || ir != jr ) return 0;

	LM mask[4];
	mask[0] = ss.aorb[iorb];
	mask[1] = ss.aorb[jorb]; 
	mask[2] = ss.aorb[korb];
	mask[3] = ss.aorb[lorb];

	std::list<std::array<int, 4>> eri_idx {{0, 0, 0, 0}};
	for (size_t pass = 0;  pass < 4; pass++) {
		const auto current_size = eri_idx.size();
		for (size_t e = 0; e < current_size; e++) {
			const auto next_idx = eri_idx.front();
			if (mask[pass].M != 0) {
				std::array<int, 4> a1, a2;
				std::copy(next_idx.begin(), next_idx.end(), a1.begin());
				std::copy(next_idx.begin(), next_idx.end(), a2.begin());
				a1[pass] = 1; // m > 0 for true spherical harmonic
				a2[pass] = -1; // m < 0 for true spherical harmonic
				eri_idx.push_back(a1);
				eri_idx.push_back(a2);
			} else {
				eri_idx.push_back(next_idx);
			}
			eri_idx.pop_front(); // Very important
		}
	}

	std::complex<double> i_ijkl (0.0, 0.0);

	for (const auto &idx : eri_idx) {
		std::array<LM, 4> oset;
		std::complex< double > weight(1.0, 0.0);

		// Index is an array of 4 ints, so :

		for (size_t ipos = 0; ipos < 4; ipos++) {
			auto o = LM(mask[ipos]);
			o.M = o.M * idx[ipos];
			if (o.M != 0) {
				if (mask[ipos].M > 0 && o.M > 0) weight *= (gsl_pow_int(-1, mask[ipos].M)/ sqrt(2.));
				if (mask[ipos].M > 0 && o.M < 0) weight *= (1./sqrt(2.));
				if (mask[ipos].M < 0 && o.M > 0) weight *= std::complex< double >(0.0, -gsl_pow_int(-1, mask[ipos].M)/ sqrt(2.));
				if (mask[ipos].M < 0 && o.M < 0) weight *= std::complex< double >(0.0, 1./sqrt(2.));
			}
			oset[ipos] = o;
		}

		// ------- Make 1st and 3rd orbitals complex conjugate
		double phase  = gsl_pow_int(-1.0, oset[0].M + oset[2].M);
		oset[0].M *= -1;
		oset[2].M *= -1;
		// -------

		//i_ijkl += weight * r12.eval_simple(g.r[ir], g.r[kr], oset[0], oset[1], oset[2], oset[3]);
		i_ijkl += phase * weight * r12.eval_simple(g.r[ir], g.r[kr], oset[0], oset[1], oset[2], oset[3]);

	}
#ifdef _OPENMP
        if (omp_in_parallel()) {
            #pragma omp critical
            if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx_cached, std::real(i_ijkl)));
        } else {
            if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx_cached, std::real(i_ijkl)));
        }
#else 
            if (cached_eri.size() < max_cache_size) cached_eri.insert(std::make_pair(eri_idx_cached, std::real(i_ijkl)));
#endif


	//return r12.eval_simple(g.r[ir], g.r[kr], ss.aorb[iorb], ss.aorb[jorb], ss.aorb[korb], ss.aorb[lorb]);
	return std::real(i_ijkl);

}

double Grid_integrals::ce_ref(size_t i, size_t j, size_t k, size_t l) {

	// This function was created for benchmarking purposes; should produce the same result as ce

	// Independent way to evaluate the same "integral":
	// 1. Inegrate out the first coordinate by expanding the (ij) density and applying Laplace formula
	// 2. Pack the result of step 1 so it is represented on a grid
	// 3. Contract the (kl) density with the function obtained after step 2 
	
	// Make sure that the grid is large enough
	
	assert (g.L_max >= 3*ss.L_max);

	// Note: is it even correct to use all the spherical harmonics with L up to g.L_max?
	
	std::vector<double> den_ij(g.nrad * g.nang, 0.0), den_kl(g.nrad * g.nang, 0.0);

	for (size_t r = 0; r < g.nrad; r++) {
		for (size_t a = 0; a < g.nang; a++) {
			den_ij[r * g.nang + a] = paux_bf[i * g.nrad * g.nang +  r * g.nang + a] * paux_bf[j * g.nrad * g.nang +  r * g.nang + a];
			den_kl[r * g.nang + a] = paux_bf[k * g.nrad * g.nang +  r * g.nang + a] * paux_bf[l * g.nrad * g.nang +  r * g.nang + a];
		}
	}

	// Storage convention: grid representation of each radial function is contiguous in memory

	std::vector< std::complex< double > > rden_ij((g.L_max + 1) * (g.L_max + 1) * g.nrad), pot_ij((g.L_max + 1) * (g.L_max + 1) * g.nrad);

	for (int l = 0; l < g.L_max + 1; l++) {
		for (int m = -l; m < l + 1; m++) {
			int idx = l * l + l + m; // See the ShellSet definition for details
			for (size_t r = 0; r < g.nrad ;r++ ) {
				rden_ij[g.nrad * idx + r ] = std::complex< double > (0.0, 0.0);
				for (size_t a = 0; a < g.nang; a++) {
					auto [th, p] = g.thetaphi_ang[a];
					auto [re, im] = Y(l, m, th, p);
					// Note complex conjugation below
					rden_ij[g.nrad * idx + r] += std::complex< double >(re, -im) * den_ij[r * g.nang + a] * 4. * M_PI * g.gridw_a[a];
				}
			}
		}
	}

	// Apply Laplace expansion to calculate the potential 
	double rlarge, rsmall;
	
	for (int l =0 ; l < g.L_max + 1; l++) {
		for (int m= -l; m < l + 1; m++) {
			int idx = l * l + l + m; // See the ShellSet definition for details
			for (size_t r = 0; r < g.nrad ;r++ ) {
				pot_ij[idx * g.nrad + r] = std::complex< double >(0.0, 0.0);
				for (size_t r1 = 0; r1 < g.nrad; r1++) {
					rlarge = (r1 >= r ? g.r[r] : g.r[r1]);
					rsmall = (r1 >= r ? g.r[r1] : g.r[r]);
					pot_ij[idx * g.nrad + r] += 4. * M_PI / double(2*l + 1)  * gsl_pow_int(rsmall/rlarge, l) / rlarge * rden_ij[idx * g.nrad + r1] * g.gridw_r[r1];
				}
			}
		}
	}

	// Finally, contract with (kl) density; the latter will be constructed on the fly
	std::complex < double > i_ijkl (0.0, 0.0);
	for (int l =0 ; l < g.L_max + 1; l++) {
		for (int m= -l; m < l + 1; m++) {
			int idx = l * l + l + m; // See the ShellSet definition for details
			for (size_t r =0; r < g.nrad; r++) {
				for (size_t a = 0; a < g.nang; a++) {
					auto [th, p] = g.thetaphi_ang[a];
					auto [re, im] = Y(l, m, th, p);
					i_ijkl += 4. * M_PI * g.gridw_r[r] * g.gridw_a[a] * std::complex< double >(re, im) * pot_ij[idx * g.nrad + r] * den_kl[r * g.nang + a];
				}
			}
		}
	}


	return std::real(i_ijkl);

}

double Grid_integrals::hc(size_t i, size_t j) {

	double h = 0.0;
	// Assymetric kinetic part
	h += 0.5 * (ke(i, j) + ke(j, i));
	// Symmetric electron nucleus attraction
	h += vn(i, j);

	return h;
}

double Grid_integrals::ke(size_t i, size_t j) {

	// This function may need to be redesigned
	// to use std::complex instead of handling 
	// all the complex arithmetic manually

    size_t ngrid = g.nrad * g.nang;

    double m_re = 0.0, m_im = 0.0; // Matrix element

    arma::vec orb_i(&paux_bf[i * ngrid], ngrid, false);
    arma::vec orb_j(&paux_bf[j * ngrid], ngrid, false);

	std::vector<double> Ri_re(g.nrad * ss.size(), 0.0), 
		                Ri_im(g.nrad * ss.size(), 0.0),
                        Rj_re(g.nrad * ss.size(), 0.0),
		                Rj_im(g.nrad * ss.size(), 0.0);

	// Perform spherical harmonic expansion for both orbitals 
	
	for (size_t o = 0; o < ss.size(); o++) {
		auto &L = ss.aorb[o].L, &M = ss.aorb[o].M;
		for ( size_t r = 0; r < g.nrad ; r++ ) {
			double di_re = 0.0, di_im = 0.0, dj_re = 0.0, dj_im = 0.0;
			for ( size_t a = 0; a < g.nang ; a++ ) {

				auto [th, p] = g.thetaphi_ang[a];
				auto [re, im] = Y(L, M, th, p);

				// Note: we are multiplying by the complex conjugate; orbitals are real

				di_re += 4. * M_PI * g.gridw_a[a] * re * orb_i[r * g.nang + a];
				di_im -= 4. * M_PI * g.gridw_a[a] * im * orb_i[r * g.nang + a];
				dj_re += 4. * M_PI * g.gridw_a[a] * re * orb_j[r * g.nang + a];
				dj_im -= 4. * M_PI * g.gridw_a[a] * im * orb_j[r * g.nang + a];

			}

			Ri_re[o * g.nrad + r] = di_re;
			Ri_im[o * g.nrad + r] = di_im;
			Rj_re[o * g.nrad + r] = dj_re;
			Rj_im[o * g.nrad + r] = dj_im;
		}
	}

	// Calculate the kinetic energy integral
	
	std::vector<double> lapl_re(g.nrad, 0.0), lapl_im(g.nrad, 0.0);
	
	for ( size_t o = 0; o < ss.size(); o++) {
		std::fill(lapl_re.begin(), lapl_re.end(), 0.0);
		std::fill(lapl_im.begin(), lapl_im.end(), 0.0);
		int L = ss.aorb[o].L;
		lp.apply_fortran(&Rj_re[o * g.nrad], lapl_re.data()); 
		lp.apply_fortran(&Rj_im[o * g.nrad], lapl_im.data()); 
		for ( size_t r  = 0; r < g.nrad; r++) {

			double t_re = 0.0, t_im = 0.0;

			t_re = lapl_re[r] - (double(L * (L + 1))/ gsl_pow_2(g.r[r]) * Rj_re[o * g.nrad + r]);
			t_im = lapl_im[r] - (double(L * (L + 1))/ gsl_pow_2(g.r[r]) * Rj_im[o * g.nrad + r]);
			m_re += (g.gridw_r[r] * (Ri_re[o * g.nrad + r] * t_re + Ri_im[o * g.nrad + r] * t_im));
			m_im += (g.gridw_r[r] * (Ri_re[o * g.nrad + r] * t_im - Ri_im[o * g.nrad + r] * t_re));

		}
	}

	assert ( std::abs(m_im) < 1e-8);
	
	return -0.5 * m_re;

}

double Grid_integrals::vn(size_t i, size_t j) {

	double h = 0.0;
	int ngrid = g.nrad * g.nang;

	arma::vec orb_i(&paux_bf[i * ngrid], ngrid, false);
	arma::vec orb_j(&paux_bf[j * ngrid], ngrid, false);

	for ( size_t k = 0; k < g.nrad ; k++) 
		for ( size_t l = 0; l < g.nang; l++) 
			h += (-Znuc / g.r[k] * orb_j[k * g.nang + l] * orb_i[k * g.nang + l] * g.gridw_r[k] *  g.gridw_a[l] * 4. * M_PI);


	return h;

}
