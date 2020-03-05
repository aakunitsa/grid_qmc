#ifndef QINT_H
#define QINT_H

#include <map>
#include <unordered_map>
#include <string>
#include "qgrid.h"
#include "qorbitals.h"
#include "qparams_reader.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

class Integral_factory {

	public:

		Integral_factory(std::map<string, int> &p); // Make the second parameter optional?

		// Member functions
		virtual void fcidump();

		// Virtual member functions
		virtual double hc(size_t i, size_t j) {return 0.0;}; // Evaluates one particle kinetic energy integral 
		virtual double ce(size_t i, size_t j, size_t k, size_t l) {return 0.0;}; // Evaluates coulomb integral in chemists notation; i.e. (ij|kl)

                std::vector<double> expand(std::vector<double> &vec); // Expands a vector wrt to the orbital basis defined in paux_bf

		// Member variables
                Becke_grid g;
		int n1porb;
		std::vector<double> paux_bf;

	protected:

		// Member functions
		bool check_orthogonality(double orth_thresh);

		// Member variables
                Laplacian lp;

                // Memoize integrals
                size_t max_cache_size;
                std::unordered_map<size_t, double> cached_eri;
                // Function to convert ERI index into a number
                size_t encode(size_t i, size_t j, size_t k, size_t l);
                bool test_encode();
                // Reorder basis functions accrording to their core energies (i.e. 2e part is ignored)
                // This will only be used with fcidump for now...
                std::vector<size_t> paux_bf_map; 
                void gen_bf_map(bool sort);
                int me; // to control printout when running with MPI
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
		ShellSet &ss;
		std::vector<double> aux_bf;
		double Znuc;
		string ofile; // Orbital file name 

		// Member functions 
		// WARNING: read_orbs needs to be refactored and should not be used for now
		bool read_orbs(); // attempts to read orbitals; returns false on failure


};

class Saved_integrals : public Integral_factory {

	// This class will just read FCIDUMP from disc
	// and provide integrals upon request;
	
	public:

		// Member functions
		Saved_integrals(Params_reader &pr);

		double hc(size_t i, size_t j); 
		double ce(size_t i, size_t j, size_t k, size_t l);

		bool good;

	private:

		// Member variables
		std::map<int, double> eri, hcore;  

		// Member functions
		bool read_fcidump(string &int_file);

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
		ShellSet &ss;

		// Member functions
		inline std::tuple<size_t, size_t> unpack_orb_index(size_t i) {
			size_t iorb = i % ss.size(), ir = (i - iorb) / ss.size();
			return std::make_tuple(ir, iorb);
		}

};


#endif  
