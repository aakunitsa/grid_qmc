#ifndef QINT_H
#define QINT_H

#include <map>
#include <string>
#include "qgrid.h"
#include "qorbitals.h"
#include "qparams_reader.h"

class Integral_factory {

	public:

		Integral_factory(std::map<string, int> &p, ShellSet &orb);

		// Member functions
		void fcidump();

		// Virtual member functions
		virtual double hc(size_t i, size_t j) {return 0.0;}; // Evaluates one particle kinetic energy integral 
		virtual double ce(size_t i, size_t j, size_t k, size_t l) {return 0.0;}; // Evaluates coulomb integral in chemists notation; i.e. (ij|kl)

		// Member variables
	    Becke_grid g;
		int n1porb;
		std::vector<double> paux_bf;

	protected:

		// Member functions
		bool check_orthogonality(double orth_thresh);

		// Member variables
		ShellSet &ss;
	    Laplacian lp;

};


class Aux_integrals : public Integral_factory {
	// This class can either read orbital coefficients
	// from disk and calculate the integrals using Poisson 
	// solver or read the inegral file directly

	public:
		// Member functions 
		Aux_integrals(Params_reader &pr, ShellSet &orb);
		~Aux_integrals(); // Clean up after Poisson solver
	    void gen_aux_basis(double lthresh, double othresh); // Putting this here temporarily

		double hc(size_t i, size_t j); 
		double ce(size_t i, size_t j, size_t k, size_t l);

		double ke(size_t i, size_t j); 
		double vn(size_t i, size_t j); 


	private: 

		// Member variables
		int naux;
		std::vector<double> aux_bf;
		double Znuc;
		string ofile; // Orbital file name 

		// Member functions 
		// WARNING: read_orbs needs to be refactored and should not be used for now
		bool read_orbs(); // attempts to read orbitals; returns false on failure

};

class Grid_integrals : public Integral_factory {

	// Don't forget to issue a warning if the L_max for the
	// grid is not large enough

	public:

		// Member functions 
		Grid_integrals(std::map<string, int> &p, ShellSet &orb);

		double hc(size_t i, size_t j); 
		double ce(size_t i, size_t j, size_t k, size_t l);

		double ce_ref(size_t i, size_t j, size_t k, size_t l);
		double ke(size_t i, size_t j); 
		double vn(size_t i, size_t j); 

	private:

		// Member variables
		Coulomb r12;
		double Znuc;
		std::vector<double> aux_bf;

		// Member functions
		inline std::tuple<size_t, size_t> unpack_orb_index(size_t i) {
			size_t iorb = i % ss.size(), ir = (i - iorb) / ss.size();
			return std::make_tuple(ir, iorb);
		}

};


#endif  
