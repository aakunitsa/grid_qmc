#ifndef QHAM_H
#define QHAM_H

#include "qorbitals.h"
#include "qgrid.h"
#include <map>
#include <string>


class Hamiltonian {

    public:

        Hamiltonian(std::map<string, int> &p, ShellSet &orb);
		~Hamiltonian() {
			for (auto p : basis) delete [] p;
			basis.resize(0);
		}

        //void build_basis();
        //vector<double> diag();

        int get_basis_size() { return nel * basis.size(); }

        vector<double *> basis;


    private:

        ShellSet &ss;
		Becke_grid g;
		Laplacian lp;
		Coulomb r12;
		size_t nel, nalpha, nbeta;


};



#endif
