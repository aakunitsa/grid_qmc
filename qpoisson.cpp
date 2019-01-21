#include <gsl/gsl_math.h>
#include "qpoisson.h"
#include <armadillo>
#include <fstream>
#include <cstdlib>

// Note : this lapack subroutine destroys matrix A !!!
extern "C" void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void second_deriv2_(double *R, int *IA, int *NRAD);
extern "C" void second_deriv2_debug_(double *R, int *IA, int *NRAD, double *DD1, double *DD2);

Poisson_solver::Poisson_solver(std::map<string, int> &p) : g(p), ss(g.L_max) {


    // Auxiliary arrays
    
    // Set up a stencil for a given grid
    poi_lhs.resize(g.nrad * g.nrad); // Stencil matrix

    poi_rhs.resize(g.nrad);

    rrho_re.resize(ss.size() * g.nrad);
    rrho_im.resize(ss.size() * g.nrad);
    d1.resize(g.nrad);
    d2.resize(g.nrad);
	pc.resize(g.nrad); // temporary variable to generate stencil matrix
	ipiv.resize(g.nrad); // IPIV array; see LAPACK manual


}

void Poisson_solver::density(const std::vector<double> &rho_re, const std::vector<double> &rho_im) {

	assert ( rho_re.size() == g.nrad * g.nang);
	assert ( rho_im.size() == g.nrad * g.nang);

    // Calculate fictitious charges first
    q_re = 0.0, q_im = 0.0;

    for (size_t i = 0; i < g.nrad; i++) {
        for (size_t j = 0; j < g.nang; j++) {
            q_re += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_re[i * g.nang + j];
            q_im += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_im[i * g.nang + j];
        }
    }

	std::fill (rrho_re.begin(), rrho_re.end(), 0.0);
	std::fill (rrho_im.begin(), rrho_im.end(), 0.0);

    for (size_t i = 0; i < ss.size(); i++) {
        const auto &o = ss.aorb[i];
		int L_ = o.L, M_ = o.M;
        for (size_t ir = 0; ir < g.nrad; ir++) {
            double d_re = 0.0, d_im = 0.0;
            for (size_t ia = 0; ia < g.nang; ia++) {
                auto [th, p] = g.thetaphi_ang[ia];
				assert (g.thetaphi_ang[ia][0] == th && g.thetaphi_ang[ia][1] == p);
				assert(o.L == L_);
				assert(o.M == M_);
                auto [y_re, y_im] = Y(o.L, o.M, th, p);
				//assert ((g.gridw_a[ia] > 0.) && (rho_re[ir * g.nang + ia] >= 0.));
				assert (g.gridw_a[ia] > 0.);
				if (o.L == 0) assert (y_re > 0);
                d_re += 4. * M_PI * g.gridw_a[ia] * (y_re * rho_re[ir * g.nang + ia] + y_im * rho_im[ir * g.nang + ia]);
                d_im += 4. * M_PI * g.gridw_a[ia] * (y_re * rho_im[ir * g.nang + ia] - y_im * rho_re[ir * g.nang + ia]);
            }
            rrho_re[i * g.nrad + ir] = d_re;
            rrho_im[i * g.nrad + ir] = d_im;
        }

		/*std::cout << " Printing the radial densities for orbital " << i << std::endl;

		for (size_t jr = 0; jr < g.nrad; jr++) {
			std::cout << rrho_re[i * g.nrad + jr] << " ";
		}

		std::cout << std::endl;
		std::cout <<  "------ end of real part ------" << std::endl;

		for (size_t jr = 0; jr < g.nrad; jr++) {
			std::cout << rrho_im[i * g.nrad + jr] << " ";
		}

		std::cout << std::endl;
		std::cout <<  "------ end of imaginary part ------" << std::endl; */
    }

	//printf("q_re = %18.10f q_im = %18.10f \n", q_re, q_im); // Added for debugging purposes; looks fine
}

void Poisson_solver::potential(std::vector<double> &el_pot_re, std::vector<double> &el_pot_im) {

	assert (el_pot_re.size() == g.nrad * g.nang);
	assert (el_pot_im.size() == g.nrad * g.nang);

	std::fill(el_pot_re.begin(), el_pot_re.end(), 0.0); // Very important to make sure that there is nothing sitting in the 
	std::fill(el_pot_im.begin(), el_pot_im.end(), 0.0); // input vectors in the first place...

    for (size_t i = 0; i < ss.size(); i++) {
        const auto &o = ss.aorb[i];

        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);

        for(size_t j = 0; j < g.nrad; j++) {
			std::fill(pc.begin(), pc.end(), 0.0); pc[j] = 1.0;
            //second_deriv(pc.data());
			int iat = 2, nrad = g.nrad;
            second_deriv2_(pc.data(), &iat, &nrad);
			// After the previous call array p contains the j-th column of the stencil matrix
			// Transfer the contents of p to poi_rhs; this can be done directly since LAPACK 
			// assumes the column major ordering of the matrix: 
			std::copy (pc.begin(), pc.end(), std::next(poi_lhs.begin(), j * g.nrad));
            poi_lhs[j * g.nrad + j] -= double(o.L * (o.L + 1)) / gsl_pow_2(g.r[j]);
            poi_rhs[j] = -4. * M_PI * g.r[j] * rrho_re[ i * g.nrad + j ];
        }


        {
			int info, nrhs = 1, N = int(g.nrad);
			dgesv_(&N, &nrhs, poi_lhs.data(), &N, ipiv.data(), poi_rhs.data(), &N, &info);
            assert(info == 0); // If it is not zero - something went wrong

			const auto &x = poi_rhs;


            // Go over the combined radial+angular grid and include the appropriate contributions
            // to the electrostatic potential
            for (size_t ir = 0; ir < g.nrad; ir++) {
                double c = gsl_pow_int(g.r[ir], -1) * x[ir];
                for (size_t ia = 0; ia < g.nang; ia++) {
                    auto [th, p] = g.thetaphi_ang[ia];
                    auto [y_re, y_im] = Y(o.L, o.M, th, p);
                    el_pot_re[ir * g.nang + ia] += c * y_re;
                    el_pot_im[ir * g.nang + ia] += c * y_im;
                }
            }
        }

        // reset all the relevant matrices (ie the lhs and the rhs of the Poisson equation)
		
        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);

        for(size_t j = 0; j < g.nrad; j++) {
			std::fill(pc.begin(), pc.end(), 0.0); pc[j] = 1.0;
            //second_deriv(pc.data());
			int iat = 2, nrad = g.nrad;
            second_deriv2_(pc.data(), &iat, &nrad);
			// After the previous call array p contains the j-th column of the stencil matrix
			// Transfer the contents of p to poi_rhs; this can be done directly since LAPACK 
			// assumes the column major ordering of the matrix: 
			std::copy (pc.begin(), pc.end(), std::next(poi_lhs.begin(), j * g.nrad));
            poi_lhs[j * g.nrad + j] -= double(o.L * (o.L + 1)) / gsl_pow_2(g.r[j]);
            poi_rhs[j] = -4. * M_PI * g.r[j] * rrho_im[ i * g.nrad + j ];
        }

        {

			int info, nrhs = 1, N = int(g.nrad);
			dgesv_(&N, &nrhs, poi_lhs.data(), &N, ipiv.data(), poi_rhs.data(), &N, &info);
            assert(info == 0); // If it is not zero - something went wrong

			const auto &x = poi_rhs;

            // Go over the combined radial+angular grid and include the appropriate contributions
            // to the electrostatic potential
            for (size_t ir = 0; ir < g.nrad; ir++) {
                double c = gsl_pow_int(g.r[ir], -1) * x[ir];
                for (size_t ia = 0; ia < g.nang; ia++) {
                    auto [th, p] = g.thetaphi_ang[ia];
                    auto [y_re, y_im] = Y(o.L, o.M, th, p);
                    el_pot_re[ir * g.nang + ia] -= c * y_im;
                    el_pot_im[ir * g.nang + ia] += c * y_re;
                }
            }
        }
    }
}

void Poisson_solver::second_deriv_debug(double *f, double *dd1, double *dd2) {

	// Note: this subroutine is different from the one used in the Laplacian 
	// class as it implements the calculation of the radial second derivative; not the actual 
	// Laplcian
	// Finite-difference formulas were taken from Polymer and are due to So Hirata

    size_t &N_G = g.nrad;
    double h = M_PI/(N_G + 1), h2 = gsl_pow_2(h);
     


    dd1[0] = (-1764.*f[0]+4320.*f[1]-5400.*f[2]+4800.*f[3]-2700.*f[4]+864.*f[5]-120.*f[6])/(720.*h);
    //R1D(1)=(-1764.0D0*R(1)+4320.0D0*R(2)-5400.0D0*R(3)+4800.0D0*R(4) -2700.0D0*R(5)+864.0D0*R(6)-120.0D0*R(7))/(720.0D0*H)

    dd1[1] = (-120.*f[0]-924.*f[1]+1800.*f[2]-1200.*f[3]+600.*f[4]-180.*f[5]+24.*f[6])/(720.*h);
    //R1D(2)=(-120.0D0*R(1)-924.0D0*R(2)+1800.0D0*R(3)-1200.0D0*R(4) +600.0D0*R(5)-180.0D0*R(6)+24.0D0*R(7))/(720.0D0*H)

    dd1[2] = (24.*f[0]-288.*f[1]-420.*f[2]+960.*f[3]-360.*f[4]+96.*f[5]-12.*f[6])/(720.*h);
    //R1D(3)=(24.0D0*R(1)-288.0D0*R(2)-420.0D0*R(3)+960.0D0*R(4) -360.0D0*R(5)+96.0D0*R(6)-12.0D0*R(7))/(720.0D0*H)

    dd1[3] = (-12.*f[0]+108.*f[1]-540.*f[2]+540.*f[4]-108.*f[5]+12.*f[6]) / (720.*h);
    //R1D(4)=(-12.0D0*R(1)+108.0D0*R(2)-540.0D0*R(3) +540.0D0*R(5)-108.0D0*R(6)+12.0D0*R(7))/(720.0D0*H)


    dd1[N_G - 4] = (-12.*f[N_G-7]+108.*f[N_G-6]-540.*f[N_G-5]+540.*f[N_G-3]-108.*f[N_G-2]+12.*f[N_G-1]) / (720.*h);
    //R1D(RG-3)=(-12.0D0*R(RG-6)+108.0D0*R(RG-5)-540.0D0*R(RG-4) +540.0D0*R(RG-2)-108.0D0*R(RG-1)+12.0D0*R(RG))/(720.0D0*H)

    dd1[N_G - 3] = (-12.*f[N_G-6] +108.*f[N_G-5] -540.*f[N_G-4]   +540.*f[N_G-2]  -108.*f[N_G-1]) / (720.*h);
    //R1D(RG-2)=(-12.0D0*R(RG-5)+108.0D0*R(RG-4)-540.0D0*R(RG-3) +540.0D0*R(RG-1)-108.0D0*R(RG))/(720.0D0*H)

    dd1[N_G - 2] = (12.*f[N_G-6] -96.*f[N_G-5] +360.*f[N_G-4]-960.*f[N_G-3] +420.*f[N_G-2] +288.*f[N_G-1]) / (720.*h);
    //R1D(RG-1)=(+12.0D0*R(RG-5)-96.0D0*R(RG-4)+360.0D0*R(RG-3)-960.0D0*R(RG-2) +420.0D0*R(RG-1)+288.0D0*R(RG))/(720.0D0*H)

    dd1[N_G - 1] = (-24.*f[N_G-6]+180.*f[N_G-5]-600.*f[N_G-4] +1200.*f[N_G-3]-1800.*f[N_G-2]+924.*f[N_G-1]) / (720.*h);
    //R1D(RG)=(-24.0D0*R(RG-5)+180.0D0*R(RG-4)-600.0D0*R(RG-3)+1200.0D0*R(RG-2) -1800.0D0*R(RG-1)+924.0D0*R(RG))/(720.0D0*H)
	

    
    dd2[0] = (1624.*f[0]-6264.*f[1]+10530.*f[2]-10160.*f[3]+5940.*f[4]-1944.*f[5]+274.*f[6]) / (360.*h2);
    //R2D(1)=(1624.0D0*R(1)-6264.0D0*R(2)+10530.0D0*R(3)-10160.0D0*R(4)+5940.0D0*R(5)-1944.0D0*R(6)+274.0D0*R(7))/(360.0D0*H**2)
    dd2[1] = (274.*f[0]-294.*f[1]-510.*f[2]+940.*f[3]-570.*f[4]+186.*f[5]-26.*f[6]) / (360.*h2);
    //R2D(2)=(274.0D0*R(1)-294.0D0*R(2)-510.0D0*R(3)+940.0D0*R(4)-570.0D0*R(5)+186.0D0*R(6)-26.0D0*R(7))/(360.0D0*H**2)
    dd2[2] = (-26.*f[0]+456.*f[1]-840.*f[2]+400.*f[3]+30.*f[4]-24.*f[5]+4.*f[6]) / (360.*h2);
    //R2D(3)=(-26.0D0*R(1)+456.0D0*R(2)-840.0D0*R(3)+400.0D0*R(4)+30.0D0*R(5)-24.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)
    dd2[3] = (4.*f[0]-54.*f[1]+540.*f[2]-980.*f[3]+540.*f[4]-54.*f[5]+4.*f[6]) / (360.*h2);
    //R2D(4)=(4.0D0*R(1)-54.0D0*R(2)+540.0D0*R(3)-980.0D0*R(4)+540.0D0*R(5)-54.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)

    
    dd2[N_G - 4] = (4.*f[N_G-7]-54.*f[N_G-6]+540.*f[N_G-5]-980.*f[N_G-4]+540.*f[N_G-3]-54.*f[N_G-2]+4.*f[N_G-1]) / (360.*h2);
    //R2D(RG-3)=(4.0D0*R(RG-6)-54.0D0*R(RG-5)+540.0D0*R(RG-4)-980.0D0*R(RG-3)+540.0D0*R(RG-2)-54.0D0*R(RG-1)+4.0D0*R(RG))/(360.0D0*H**2)
    dd2[N_G - 3] = ( 4.*f[N_G-6]-54.*f[N_G-5]+540.*f[N_G-4]-980.*f[N_G-3]+540.*f[N_G-2]-54.*f[N_G-1]) / (360.*h2);
    //R2D(RG-2)=(+4.0D0*R(RG-5)-54.0D0*R(RG-4)  +540.0D0*R(RG-3)-980.0D0*R(RG-2)+540.0D0*R(RG-1)-54.0D0*R(RG))/(360.0D0*H**2)
    dd2[N_G - 2] = (+4.*f[N_G-6]-24.*f[N_G-5]  +30.*f[N_G-4]+400.*f[N_G-3]-840.*f[N_G-2]+456.*f[N_G-1]) / (360.*h2);
    //R2D(RG-1)=(+4.0D0*R(RG-5)-24.0D0*R(RG-4)   +30.0D0*R(RG-3)+400.0D0*R(RG-2)-840.0D0*R(RG-1)+456.0D0*R(RG))/(360.0D0*H**2)
    dd2[N_G - 1] = (-26.*f[N_G-6]+186.*f[N_G-5] -570.*f[N_G-4] +940.*f[N_G-3]  -510.*f[N_G-2]-294.*f[N_G-1]) / (360.*h2);
    //R2D(RG)=(-26.0D0*R(RG-5)+186.0D0*R(RG-4)  -570.0D0*R(RG-3)+940.0D0*R(RG-2)-510.0D0*R(RG-1)-294.0D0*R(RG))/(360.0D0*H**2)

    for (size_t i = 4; i < N_G - 4; i++) {
        dd1[i] = (+144.*f[i-4]-1536.*f[i-3]+8064.*f[i-2]-32256.*f[i-1]+32256.*f[i+1]-8064.*f[i+2]+1536.*f[i+3]-144.*f[i+4])/(40320.*h);
        //R1D(I)=(+144.0D0*R(I-4)-1536.0D0*R(I-3)+8064.0D0*R(I-2)-32256.0D0*R(I-1) +32256.0D0*R(I+1)-8064.0D0*R(I+2)+1536.0D0*R(I+3)-144.0D0*R(I+4))/(40320.0D0*H)
        dd2[i] = (-36.*f[i-4]+512.*f[i-3]-4032.*f[i-2]+32256.*f[i-1]-57400.*f[i]+32256.*f[i+1]-4032.*f[i+2]+512.*f[i+3]-36.*f[i+4])/(20160.*h2);
        //R2D(I)=(-36.0D0*R(I-4) +512.0D0*R(I-3)-4032.0D0*R(I-2)+32256.0D0*R(I-1)-57400.0D0*R(I)+32256.0D0*R(I+1)-4032.0D0*R(I+2)+512.0D0*R(I+3)-36.0D0*R(I+4))/(20160.0D0*H**2)
    }


    for (size_t i = 0; i < N_G; i++) {
        double z = 3.14159265358979 / (N_G + 1) * (i + 1);
        //double fac2 = gsl_pow_2(cos(z) - 1.0), fac4 = gsl_pow_4(cos(z) - 1.0);
        //double fac3 = gsl_pow_3(1 - cos(z));
        double fac2 = (cos(z) - 1.0)*(cos(z) - 1.0), fac4 = fac2 * fac2;
        double fac3 = (1.0 - cos(z)) * fac2;
		//f[i] = dd1[i] * (fac3/(2.0*g.r_at*g.r_at*sin(z)) - fac4*cos(z)/(4.0*g.r_at*g.r_at*sin(z)*sin(z)*sin(z)));
        //f[i] = fac4*(cos(z) + 2)/(4*gsl_pow_2(g.r_at)*gsl_pow_3(sin(z))) * d1[i];
        //f[i] += fac4 / (4*gsl_pow_2(g.r_at)*gsl_pow_2(sin(z))) * dd2[i];
        //f[i] += fac4 / (4*(g.r_at*sin(z))*(g.r_at*sin(z))) * dd2[i];
        //f[i] = fac4 / (4*(g.r_at*sin(z))*(g.r_at*sin(z))) * dd2[i];
        //R(I)=+R1D(I)*((1.0D0-X)**3/(2.0D0*R0**2*DSIN(OMEGA)) - (1.0D0-X)**4*DCOS(OMEGA)/(4.0D0*R0**2*DSIN(OMEGA)**3)) +R2D(I)*((1.0D0-X)**2/(2.0D0*R0*DSIN(OMEGA)))**2
		f[i] = dd1[i] * ((1.0 - cos(z)) * (1.0 - cos(z)) * (1.0 - cos(z))/(2.0*g.r_at*g.r_at*sin(z)) 
				-(1.0 - cos(z))*(1.0 - cos(z))*(1.0 - cos(z))*(1.0 - cos(z))*cos(z)/(4.0*g.r_at*g.r_at*sin(z)*sin(z)*sin(z))) 
			    +(1.0 - cos(z))*(1.0 - cos(z))*(1.0 - cos(z))*(1.0 - cos(z)) / (4*(g.r_at*sin(z))*(g.r_at*sin(z))) * dd2[i];
    }

}

void Poisson_solver::second_deriv(double *f) {

	// Note: this subroutine is different from the one used in the Laplacian 
	// class as it implements the calculation of the radial second derivative; not the actual 
	// Laplcian
	// Finite-difference formulas were taken from Polymer and are due to So Hirata

    size_t &N_G = g.nrad;
    double h = M_PI/(N_G + 1), h2 = gsl_pow_2(h);
     


    d1[0] = (-1764.*f[0]+4320.*f[1]-5400.*f[2]+4800.*f[3]-2700.*f[4]+864.*f[5]-120.*f[6])/(720.*h);
    //R1D(1)=(-1764.0D0*R(1)+4320.0D0*R(2)-5400.0D0*R(3)+4800.0D0*R(4) -2700.0D0*R(5)+864.0D0*R(6)-120.0D0*R(7))/(720.0D0*H)

    d1[1] = (-120.*f[0]-924.*f[1]+1800.*f[2]-1200.*f[3]+600.*f[4]-180.*f[5]+24.*f[6])/(720.*h);
    //R1D(2)=(-120.0D0*R(1)-924.0D0*R(2)+1800.0D0*R(3)-1200.0D0*R(4) +600.0D0*R(5)-180.0D0*R(6)+24.0D0*R(7))/(720.0D0*H)

    d1[2] = (24.*f[0]-288.*f[1]-420.*f[2]+960.*f[3]-360.*f[4]+96.*f[5]-12.*f[6])/(720.*h);
    //R1D(3)=(24.0D0*R(1)-288.0D0*R(2)-420.0D0*R(3)+960.0D0*R(4) -360.0D0*R(5)+96.0D0*R(6)-12.0D0*R(7))/(720.0D0*H)

    d1[3] = (-12.*f[0]+108.*f[1]-540.*f[2]+540.*f[4]-108.*f[5]+12.*f[6]) / (720.*h);
    //R1D(4)=(-12.0D0*R(1)+108.0D0*R(2)-540.0D0*R(3) +540.0D0*R(5)-108.0D0*R(6)+12.0D0*R(7))/(720.0D0*H)


    d1[N_G - 4] = (-12.*f[N_G-7]+108.*f[N_G-6]-540.*f[N_G-5]+540.*f[N_G-3]-108.*f[N_G-2]+12.*f[N_G-1]) / (720.*h);
    //R1D(RG-3)=(-12.0D0*R(RG-6)+108.0D0*R(RG-5)-540.0D0*R(RG-4) +540.0D0*R(RG-2)-108.0D0*R(RG-1)+12.0D0*R(RG))/(720.0D0*H)

    d1[N_G - 3] = (-12.*f[N_G-6] +108.*f[N_G-5] -540.*f[N_G-4]   +540.*f[N_G-2]  -108.*f[N_G-1]) / (720.*h);
    //R1D(RG-2)=(-12.0D0*R(RG-5)+108.0D0*R(RG-4)-540.0D0*R(RG-3) +540.0D0*R(RG-1)-108.0D0*R(RG))/(720.0D0*H)

    d1[N_G - 2] = (12.*f[N_G-6] -96.*f[N_G-5] +360.*f[N_G-4]-960.*f[N_G-3] +420.*f[N_G-2] +288.*f[N_G-1]) / (720.*h);
    //R1D(RG-1)=(+12.0D0*R(RG-5)-96.0D0*R(RG-4)+360.0D0*R(RG-3)-960.0D0*R(RG-2) +420.0D0*R(RG-1)+288.0D0*R(RG))/(720.0D0*H)

    d1[N_G - 1] = (-24.*f[N_G-6]+180.*f[N_G-5]-600.*f[N_G-4] +1200.*f[N_G-3]-1800.*f[N_G-2]+924.*f[N_G-1]) / (720.*h);
    //R1D(RG)=(-24.0D0*R(RG-5)+180.0D0*R(RG-4)-600.0D0*R(RG-3)+1200.0D0*R(RG-2) -1800.0D0*R(RG-1)+924.0D0*R(RG))/(720.0D0*H)
	

    
    d2[0] = (1624.*f[0]-6264.*f[1]+10530.*f[2]-10160.*f[3]+5940.*f[4]-1944.*f[5]+274.*f[6]) / (360.*h2);
    //R2D(1)=(1624.0D0*R(1)-6264.0D0*R(2)+10530.0D0*R(3)-10160.0D0*R(4)+5940.0D0*R(5)-1944.0D0*R(6)+274.0D0*R(7))/(360.0D0*H**2)
    d2[1] = (274.*f[0]-294.*f[1]-510.*f[2]+940.*f[3]-570.*f[4]+186.*f[5]-26.*f[6]) / (360.*h2);
    //R2D(2)=(274.0D0*R(1)-294.0D0*R(2)-510.0D0*R(3)+940.0D0*R(4)-570.0D0*R(5)+186.0D0*R(6)-26.0D0*R(7))/(360.0D0*H**2)
    d2[2] = (-26.*f[0]+456.*f[1]-840.*f[2]+400.*f[3]+30.*f[4]-24.*f[5]+4.*f[6]) / (360.*h2);
    //R2D(3)=(-26.0D0*R(1)+456.0D0*R(2)-840.0D0*R(3)+400.0D0*R(4)+30.0D0*R(5)-24.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)
    d2[3] = (4.*f[0]-54.*f[1]+540.*f[2]-980.*f[3]+540.*f[4]-54.*f[5]+4.*f[6]) / (360.*h2);
    //R2D(4)=(4.0D0*R(1)-54.0D0*R(2)+540.0D0*R(3)-980.0D0*R(4)+540.0D0*R(5)-54.0D0*R(6)+4.0D0*R(7))/(360.0D0*H**2)

    
    d2[N_G - 4] = (4.*f[N_G-7]-54.*f[N_G-6]+540.*f[N_G-5]-980.*f[N_G-4]+540.*f[N_G-3]-54.*f[N_G-2]+4.*f[N_G-1]) / (360.*h2);
    //R2D(RG-3)=(4.0D0*R(RG-6)-54.0D0*R(RG-5)+540.0D0*R(RG-4)-980.0D0*R(RG-3)+540.0D0*R(RG-2)-54.0D0*R(RG-1)+4.0D0*R(RG))/(360.0D0*H**2)
    d2[N_G - 3] = ( 4.*f[N_G-6]-54.*f[N_G-5]+540.*f[N_G-4]-980.*f[N_G-3]+540.*f[N_G-2]-54.*f[N_G-1]) / (360.*h2);
    //R2D(RG-2)=(+4.0D0*R(RG-5)-54.0D0*R(RG-4)  +540.0D0*R(RG-3)-980.0D0*R(RG-2)+540.0D0*R(RG-1)-54.0D0*R(RG))/(360.0D0*H**2)
    d2[N_G - 2] = (+4.*f[N_G-6]-24.*f[N_G-5]  +30.*f[N_G-4]+400.*f[N_G-3]-840.*f[N_G-2]+456.*f[N_G-1]) / (360.*h2);
    //R2D(RG-1)=(+4.0D0*R(RG-5)-24.0D0*R(RG-4)   +30.0D0*R(RG-3)+400.0D0*R(RG-2)-840.0D0*R(RG-1)+456.0D0*R(RG))/(360.0D0*H**2)
    d2[N_G - 1] = (-26.*f[N_G-6]+186.*f[N_G-5] -570.*f[N_G-4] +940.*f[N_G-3]  -510.*f[N_G-2]-294.*f[N_G-1]) / (360.*h2);
    //R2D(RG)=(-26.0D0*R(RG-5)+186.0D0*R(RG-4)  -570.0D0*R(RG-3)+940.0D0*R(RG-2)-510.0D0*R(RG-1)-294.0D0*R(RG))/(360.0D0*H**2)

    for (size_t i = 4; i < N_G - 4; i++) {
        d1[i] = (+144.*f[i-4]-1536.*f[i-3]+8064.*f[i-2]-32256.*f[i-1]+32256.*f[i+1]-8064.*f[i+2]+1536.*f[i+3]-144.*f[i+4])/(40320.*h);
        //R1D(I)=(+144.0D0*R(I-4)-1536.0D0*R(I-3)+8064.0D0*R(I-2)-32256.0D0*R(I-1) +32256.0D0*R(I+1)-8064.0D0*R(I+2)+1536.0D0*R(I+3)-144.0D0*R(I+4))/(40320.0D0*H)
        d2[i] = (-36.*f[i-4]+512.*f[i-3]-4032.*f[i-2]+32256.*f[i-1]-57400.*f[i]+32256.*f[i+1]-4032.*f[i+2]+512.*f[i+3]-36.*f[i+4])/(20160.*h2);
        //R2D(I)=(-36.0D0*R(I-4) +512.0D0*R(I-3)-4032.0D0*R(I-2)+32256.0D0*R(I-1)-57400.0D0*R(I)+32256.0D0*R(I+1)-4032.0D0*R(I+2)+512.0D0*R(I+3)-36.0D0*R(I+4))/(20160.0D0*H**2)
    }

	/*

	std::cout << " Printing d arrays before calculating derivatives " << std::endl;
	std::cout << " D1: " << std::endl;

	for (size_t i = 0; i < N_G; i++) {
		std::cout << d1[i] << " " ;
	}

	std::cout << std::endl;

	std::cout << " D2: " << std::endl;

	for (size_t i = 0; i < N_G; i++) {
		std::cout << d2[i] << " " ;
	}

	std::cout << std::endl;
	*/

    for (size_t i = 0; i < N_G; i++) {
        double z = M_PI / (N_G + 1) * (i + 1);
        double fac2 = gsl_pow_2(cos(z) - 1), fac4 = gsl_pow_4(cos(z) - 1);
        double fac3 = gsl_pow_3(1 - cos(z));
		f[i] = d1[i] * (fac3/(2.0*gsl_pow_2(g.r_at)*sin(z)) - fac4*cos(z)/(4.0*gsl_pow_2(g.r_at)*gsl_pow_3(sin(z))));
        //f[i] = fac4*(cos(z) + 2)/(4*gsl_pow_2(g.r_at)*gsl_pow_3(sin(z))) * d1[i];
        f[i] += fac4 / (4*gsl_pow_2(g.r_at)*gsl_pow_2(sin(z))) * d2[i];
        //R(I)=+R1D(I)*((1.0D0-X)**3/(2.0D0*R0**2*DSIN(OMEGA)) - (1.0D0-X)**4*DCOS(OMEGA)/(4.0D0*R0**2*DSIN(OMEGA)**3)) +R2D(I)*((1.0D0-X)**2/(2.0D0*R0*DSIN(OMEGA)))**2
    }

}

void Poisson_solver::test_poisson() {

	std::cout << "Poisson solver will be tested by evaluating a set of ERI-s between the orbitals of the following form: " 
		      << " psi_lm = e**(-r**2 / 2.) * Ylm, where Ylm is a spherical harmonic with L = l and M = m ; L is less or equal to L_max allowed by a chosen grid " << std::endl;

	int L_max_t = 1;
	assert ( L_max_t <= g.L_max );
	ShellSet st(L_max_t);

	std::vector<double> prho_re (g.nang * g.nrad), prho_im (g.nang * g.nrad);
	std::vector<double> pot_re (g.nang * g.nrad), pot_im (g.nang * g.nrad);

	// Integrals will be calculated following chemists notation
	// (i1i2|i3i4)
	for (size_t i1 = 0; i1 < st.size(); i1++) {
		for (size_t i2 = 0; i2 < st.size(); i2++) {
			for (size_t i3 = 0; i3 < st.size(); i3++) {
				for (size_t i4 = 0; i4 < st.size(); i4++) {
					// treat i1 i2 as the first orbital pair
					// represent density cc (psi_i1) * psi_i2 on the grid
					LM &o1 = st.aorb[i1], &o2 = st.aorb[i2];
					// Just in case
					std::fill(prho_re.begin(), prho_re.end(), 0.0);
					std::fill(prho_im.begin(), prho_im.end(), 0.0);
					//std::cout << "Preparing tmp density...";
					for (size_t ir = 0; ir < g.nrad; ir++) {
						double rho_rad = exp(-g.r[ir] * g.r[ir]); 
						for ( size_t ia = 0; ia < g.nang; ia++) {
							auto [th, p] = g.thetaphi_ang[ia];
							auto [y_re1, y_im1] = Y(o1.L, o1.M, th, p);
							auto [y_re2, y_im2] = Y(o2.L, o2.M, th, p);
							// Real part of the product of two spherical harmonics
							// y1_re * y2_re + y1_im * y2_im
							// Imaginary:
							// y1_re * y2_im - y1_im * y2_re
							prho_re[ir * g.nang + ia] = rho_rad * (y_re1 * y_re2 + y_im1 * y_im2);
							prho_im[ir * g.nang + ia] = rho_rad * (y_re1 * y_im2 - y_im1 * y_re2);
						}
					}

					//std::cout << " Done! " << std::endl;
					//std::cout << "Submitting the density to the solver... ";

					density(prho_re, prho_im);

					//std::cout << " Done! " << std::endl;
					//std::cout << "Starting the solver... ";
					//std::cout.flush();

					potential(pot_re, pot_im);

					//std::cout << " potential has been evaluated!" << std::endl;

					double eri_re = 0.0, eri_im = 0.0;

					LM &o3 = st.aorb[i3], &o4 = st.aorb[i4];

					for (size_t ir = 0; ir < g.nrad; ir++) {
						double rho_rad = exp(-g.r[ir] * g.r[ir]); 
						for ( size_t ia = 0; ia < g.nang; ia++) {
							auto [th, p] = g.thetaphi_ang[ia];
							auto [y_re1, y_im1] = Y(o3.L, o3.M, th, p);
							auto [y_re2, y_im2] = Y(o4.L, o4.M, th, p);
							// Real part of the product of two spherical harmonics
							// y1_re * y2_re + y1_im * y2_im
							// Imaginary:
							// y1_re * y2_im - y1_im * y2_re
							eri_re += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_re[ir * g.nrad + ia] - 
									              (y_re1 * y_im2 - y_im1 * y_re2) * pot_im[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];

							eri_im += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_im[ir * g.nrad + ia] +
									              (y_re1 * y_im2 - y_im1 * y_re2) * pot_re[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];
						}
					}

					printf("Re (%d%d|%d%d) = %18.10f \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_re);
					printf("Im (%d%d|%d%d) = %18.10f \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_im);

				}
			}
		}
	}

	std::cout << "ERI calculation has been completed successfully! " << std::endl;

}

void Poisson_solver::test_stencil() {
	// The function performs numerical second derivatives of 
	// functions using two different approaches - direct finite difference (as implemented in Laplacian class)
	// and stencil matrix based calculation - and compares the accuracy of the two
	//
	
	std::cout << " Testing the correctness of the stencil matrix calculation " << std::endl;

    double tol = 1e-10;	
	// Prepare stencil matrix; I will call it A
	
    std::vector<double> st(g.nrad * g.nrad); 
    std::fill(st.begin(), st.end(), 0.0);
    for (size_t i = 0; i < g.nrad; i++) {
        st[i * g.nrad + i] = 1.0;
		int iat = 2, nrad = g.nrad;
        //second_deriv(st.data() + i * g.nrad);
        second_deriv2_(st.data() + i * g.nrad, &iat, &nrad);
    }

	arma::mat A(st.data(), g.nrad, g.nrad, false);
	//std::cout << " The stencil matrix (calculated by poisson solver) is printed below for reference purposes: " << std::endl; // Can generate a lot of data if the matrix is large..
	//A.print();
    
	// Perform tests for hydrogen atom states first
	std::vector<double> psi1(g.nrad, 0.0), psi2(g.nrad, 0.0), psi3(g.nrad, 0.0);
	std::transform(psi1.begin(), psi1.end(), psi1.begin(), psi_1s_H);
	std::transform(psi2.begin(), psi2.end(), psi2.begin(), psi_2s_H);
	std::transform(psi3.begin(), psi3.end(), psi3.begin(), psi_3s_H);

	arma::vec p1(psi1.data(), g.nrad, false), p2(psi2.data(), g.nrad, false), p3(psi3.data(), g.nrad, false);

	std::vector<double> lpsi1(g.nrad, 0.0), lpsi2(g.nrad, 0.0), lpsi3(g.nrad, 0.0);

	std::copy(psi1.begin(), psi1.end(), lpsi1.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3.begin());

	// Calculating second derivatives using second_deriv method of Poisson_solver class
	second_deriv(lpsi1.data());
	second_deriv(lpsi2.data());
	second_deriv(lpsi3.data());

	arma::vec x1(lpsi1.data(), g.nrad, false);
	arma::vec x2(lpsi2.data(), g.nrad, false);
	arma::vec x3(lpsi3.data(), g.nrad, false);

	// Calculating second derivatives with the stencil matrix; That procedure amounts to the simple matrix-vector 
	// multiplication
	
	arma::vec x1_s = A * p1, x2_s = A * p2, x3_s =  A * p3;
	double err1 = arma::max(arma::abs(x1_s - x1)), 
		   err2 = arma::max(arma::abs(x2_s - x2)),
		   err3 = arma::max(arma::abs(x3_s - x3));

	printf(" Maximum error for 1S function is %18.10f \n", err1);
	printf(" Maximum error for 2S function is %18.10f \n", err2);
	printf(" Maximum error for 3S function is %18.10f \n", err3);

	std::cout << " Comparing second_deriv2 subroutine from Polymer with second_deriv method of the Poisson_solver class " << std::endl;

	std::vector<double> lpsi1_f(g.nrad), lpsi2_f(g.nrad), lpsi3_f(g.nrad);
	std::copy(psi1.begin(), psi1.end(), lpsi1_f.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2_f.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3_f.begin());

	int ia = 2;
	int nrad = int(g.nrad);

	second_deriv2_(lpsi1_f.data(), &ia, &nrad);
	second_deriv2_(lpsi2_f.data(), &ia, &nrad);
	second_deriv2_(lpsi3_f.data(), &ia, &nrad);

	arma::vec x1_f(lpsi1_f.data(), g.nrad, false);
	arma::vec x2_f(lpsi2_f.data(), g.nrad, false);
	arma::vec x3_f(lpsi3_f.data(), g.nrad, false);

	double err1_f = arma::max(arma::abs(x1_f - x1)), 
		   err2_f = arma::max(arma::abs(x2_f - x2)),
		   err3_f = arma::max(arma::abs(x3_f - x3));

	printf(" Maximum error for 1S function is %18.10f (Fortran) \n", err1_f);
	printf(" Maximum error for 2S function is %18.10f (Fortran) \n", err2_f);
	printf(" Maximum error for 3S function is %18.10f (Fortran) \n", err3_f);

	std::cout << " Stencil matrix generated by polymer will be read from the text file stencil.dat " << std::endl;

	fstream input;
	input.open("stencil.dat");
	std::vector<double> ref_st(g.nrad * g.nrad);
	if ( input.is_open() ) {
		for (size_t i = 0; i < g.nrad; i++ ) {
			for (size_t j = 0; j < g.nrad; j++)  {
				input >> ref_st[i + j * g.nrad ]; // Transpose while reading 
			}
		}
	}

	arma::vec ref_A_vec(ref_st.data(), g.nrad * g.nrad, false),
		          A_vec(st.data(), g.nrad * g.nrad, false);

	
	std::cout << " Comparing the two stencil matrices " << std::endl;

	double err = arma::max(arma::abs(ref_A_vec - A_vec)); 

	printf(" The stencil matrix has maximum error of %18.10f \n", err);
	

}

void Poisson_solver::test_against_poly() {

	// This function reads reference data produced by polymer code 
	// and compares it the results generated by the Poisson_solver class
	// It is assumed that the calculations were performed on the equvalent 
	// grids (up to rigid rotation)

	
	std::vector<double> ref_rho_re(g.nrad * g.nang, 0.0), ref_rho_im(g.nrad * g.nang, 0.0); // Reference data is assumed to be real
	std::vector<double> ref_pot_re(g.nrad * g.nang, 0.0), ref_pot_im(g.nrad * g.nang, 0.0); 
	std::vector<double> pot_re(g.nrad * g.nang), pot_im(g.nrad * g.nang);
	
	fstream input, output;
	input.open("RHO.DAT");
	output.open("RHO_IN_CPP.DAT", ios::out);
	if (input.is_open()) {
		for (size_t i = 0; i < g.nrad; i++) {
			for (size_t j = 0; j < g.nang; j++) {
				input >> ref_rho_re[i * g.nang + j];
				assert(ref_rho_re[i * g.nang + j] >= 0.);
				output << ref_rho_re[i * g.nang + j] << std::endl;
			}
		}
	}
	input.close();
	output.close();

	input.open("POTENTIAL.DAT");
	output.open("POTENTIAL_IN_CPP.DAT", ios::out);
	if (input.is_open()) {
		for (size_t i = 0; i < g.nrad; i++) {
			for (size_t j = 0; j < g.nang; j++) {
				input >> ref_pot_re[i * g.nang + j];
				output << ref_pot_re[i * g.nang + j] << std::endl;
			}
		}
	}
	input.close();
	output.close();

	density(ref_rho_re, ref_rho_im);
	potential(pot_re, pot_im);

	arma::vec ref_apot_re(ref_pot_re.data(), g.nrad * g.nang, false),
		      ref_apot_im(ref_pot_im.data(), g.nrad * g.nang, false),
			  apot_re(pot_re.data(), g.nrad * g.nang, false),
			  apot_im(pot_im.data(), g.nrad  * g.nang, false);

	// Calculate maximum errors
	
	double err_re = arma::max(arma::abs(ref_apot_re - apot_re)),
		   err_im = arma::max(arma::abs(ref_apot_im - apot_im));

	printf(" Maximum difference between the real parts of the potential is %18.10f (%18.10f for imaginary) \n", err_re, err_im);

}

std::tuple<double, double> Poisson_solver::calc_eri(const LM &o1, const LM &o2, const LM &o3, const LM &o4) {

	std::vector<double> prho_re (g.nrad * g.nang), prho_im(g.nrad * g.nang);
	std::vector<double> pot_re (g.nrad * g.nang), pot_im(g.nrad * g.nang);

	std::fill(prho_re.begin(), prho_re.end(), 0.0);
	std::fill(prho_im.begin(), prho_im.end(), 0.0);
	for (size_t ir = 0; ir < g.nrad; ir++) {
		double rho_rad = exp(-g.r[ir] * g.r[ir]); 
		for ( size_t ia = 0; ia < g.nang; ia++) {
			auto [th, p] = g.thetaphi_ang[ia];
			auto [y_re1, y_im1] = Y(o1.L, o1.M, th, p);
			auto [y_re2, y_im2] = Y(o2.L, o2.M, th, p);
			prho_re[ir * g.nang + ia] = rho_rad * (y_re1 * y_re2 + y_im1 * y_im2);
			prho_im[ir * g.nang + ia] = rho_rad * (y_re1 * y_im2 - y_im1 * y_re2);
		}
	}

	density(prho_re, prho_im);
	potential(pot_re, pot_im);

	assert (pot_re.size() == g.nang * g.nrad);
	assert (pot_im.size() == g.nang * g.nrad);
	assert (prho_re.size() == g.nrad * g.nang);
	assert (prho_im.size() == g.nrad * g.nang);

	double eri_re = 0.0, eri_im = 0.0;
	
	for (size_t ir = 0; ir < g.nrad; ir++) {
		double rho_rad = exp(-g.r[ir] * g.r[ir]); 
		for (size_t ia = 0; ia < g.nang; ia++) {
			auto [th, p] = g.thetaphi_ang[ia];
			auto [y_re1, y_im1] = Y(o3.L, o3.M, th, p);
			auto [y_re2, y_im2] = Y(o4.L, o4.M, th, p);
			eri_re += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_re[ir * g.nrad + ia] - 
								(y_re1 * y_im2 - y_im1 * y_re2) * pot_im[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];

			eri_im += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_im[ir * g.nrad + ia] +
								(y_re1 * y_im2 - y_im1 * y_re2) * pot_re[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];
		}
	}
	
	return make_tuple(eri_re, eri_im);
}

void Poisson_solver::test_second_deriv() {
	// This method compares the intermediate quantities involved in the second derivative calculation 
	// (as implemented in the second_deriv method of the Poisson_solver class) with thier conterparts 
	// from the Polymer code; debug versions of the subroutines will be used;


	// Perform tests for hydrogen atom states first
	std::vector<double> psi1(g.nrad, 0.0), psi2(g.nrad, 0.0), psi3(g.nrad, 0.0);
	std::vector<double> tmp1(g.nrad, 0.0), tmp2(g.nrad, 0.0), tmp1_(g.nrad, 0.0), tmp2_(g.nrad, 0.0);
	std::vector<double> lpsi1_f(g.nrad), lpsi2_f(g.nrad), lpsi3_f(g.nrad);
	std::vector<double> lpsi1_c(g.nrad), lpsi2_c(g.nrad), lpsi3_c(g.nrad);

	std::copy(g.r.begin(), g.r.end(), psi1.begin());
	std::copy(g.r.begin(), g.r.end(), psi2.begin());
	std::copy(g.r.begin(), g.r.end(), psi3.begin());

	std::transform(psi1.begin(), psi1.end(), psi1.begin(), psi_1s_H);
	std::transform(psi2.begin(), psi2.end(), psi2.begin(), psi_2s_H);
	std::transform(psi3.begin(), psi3.end(), psi3.begin(), psi_3s_H);

	std::copy(psi1.begin(), psi1.end(), lpsi1_f.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2_f.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3_f.begin());

	std::copy(psi1.begin(), psi1.end(), lpsi1_c.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2_c.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3_c.begin());

	// Apply second derivative (debug version) from Poisson_solver
	
	int nrad = g.nrad, iatom = 2;
	
	second_deriv_debug(lpsi1_c.data(), tmp1.data(), tmp2.data());
	second_deriv2_debug_(lpsi1_f.data(), &iatom, &nrad, tmp1_.data(), tmp2_.data());

	{
		double max_err = 0.0, max_err1 = 0.0, max_err2 = 0.0;
		for (size_t i = 0; i < g.nrad; i++) {
			max_err = std::max(max_err, std::abs(lpsi1_f[i] - lpsi1_c[i]));
			max_err1 = std::max(max_err1, std::abs(tmp1[i] - tmp1_[i]));
			max_err2 = std::max(max_err2, std::abs(tmp2[i] - tmp2_[i]));
		}

		// Print summary 

		printf(" max derivative error is %18.10f \n", max_err);
		printf(" d1 max error is %18.10f \n", max_err1);
		printf(" d2 max error is %18.10f \n", max_err2);
	}

	std::fill (tmp1.begin(), tmp1.end(), 0.0);
	std::fill (tmp1_.begin(), tmp1_.end(), 0.0);
	std::fill (tmp2.begin(), tmp2.end(), 0.0);
	std::fill (tmp2_.begin(), tmp2_.end(), 0.0);

	second_deriv_debug(lpsi2_c.data(), tmp1.data(), tmp2.data());
	second_deriv2_debug_(lpsi2_f.data(), &iatom, &nrad, tmp1_.data(), tmp2_.data());

	{
		double max_err = 0.0, max_err1 = 0.0, max_err2 = 0.0;
		for (size_t i = 0; i < g.nrad; i++) {
			max_err = std::max(max_err, std::abs(lpsi1_f[i] - lpsi1_c[i]));
			max_err1 = std::max(max_err1, std::abs(tmp1[i] - tmp1_[i]));
			max_err2 = std::max(max_err2, std::abs(tmp2[i] - tmp2_[i]));
		}

		// Print summary 

		printf(" max derivative error is %18.10f \n", max_err);
		printf(" d1 max error is %18.10f \n", max_err1);
		printf(" d2 max error is %18.10f \n", max_err2);
	}

	std::fill (tmp1.begin(), tmp1.end(), 0.0);
	std::fill (tmp1_.begin(), tmp1_.end(), 0.0);
	std::fill (tmp2.begin(), tmp2.end(), 0.0);
	std::fill (tmp2_.begin(), tmp2_.end(), 0.0);

	second_deriv_debug(lpsi3_c.data(), tmp1.data(), tmp2.data());
	second_deriv2_debug_(lpsi3_f.data(), &iatom, &nrad, tmp1_.data(), tmp2_.data());

	{
		// The errors will be saved to file for the third test
		
		fstream err_file1, err_file2, err_file3;
		err_file1.open("ERR1.DAT", ios::out);
		err_file2.open("ERR2.DAT", ios::out);
		err_file3.open("ERR3.DAT", ios::out);

		assert( err_file1.is_open() && err_file2.is_open() && err_file3.is_open() );

		double max_err = 0.0, max_err1 = 0.0, max_err2 = 0.0;
		for (size_t i = 0; i < g.nrad; i++) {
			max_err = std::max(max_err, std::abs(lpsi3_f[i] - lpsi3_c[i]));
			max_err1 = std::max(max_err1, std::abs(tmp1[i] - tmp1_[i]));
			max_err2 = std::max(max_err2, std::abs(tmp2[i] - tmp2_[i]));
			err_file1 << i << " " << std::abs(tmp1[i] - tmp1_[i]) << std::endl;
			err_file2 << i << " " << std::abs(tmp2[i] - tmp2_[i]) << std::endl;
			err_file3 << i << " " << std::abs(lpsi3_f[i] - lpsi3_c[i]) << std::endl;
		}

		err_file1.close();
		err_file2.close();
		err_file3.close();

		// Print summary 

		printf(" max derivative error is %18.10f \n", max_err);
		printf(" d1 max error is %18.10f \n", max_err1);
		printf(" d2 max error is %18.10f \n", max_err2);
	}

}

void Poisson_solver::test_second_deriv2() {

	std::vector<double> psi1(g.nrad, 0.0), psi2(g.nrad, 0.0), psi3(g.nrad, 0.0);
	std::vector<double> tmp1(g.nrad, 0.0), tmp2(g.nrad, 0.0), tmp1_(g.nrad, 0.0), tmp2_(g.nrad, 0.0);

	std::copy(g.r.begin(), g.r.end(), psi1.begin());
	std::copy(g.r.begin(), g.r.end(), psi2.begin());
	std::copy(g.r.begin(), g.r.end(), psi3.begin());

	std::transform(psi1.begin(), psi1.end(), psi1.begin(), psi_1s_H);
	std::transform(psi2.begin(), psi2.end(), psi2.begin(), psi_2s_H);
	std::transform(psi3.begin(), psi3.end(), psi3.begin(), psi_3s_H);

	std::vector<double> lpsi1_c(g.nrad), lpsi2_c(g.nrad), lpsi3_c(g.nrad);

	std::copy(psi1.begin(), psi1.end(), lpsi1_c.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2_c.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3_c.begin());

	second_deriv(psi1.data());
	second_deriv(psi2.data());
	second_deriv(psi3.data());
	second_deriv_debug(lpsi1_c.data(), tmp1.data(), tmp2.data());
	second_deriv_debug(lpsi2_c.data(), tmp1.data(), tmp2.data());
	second_deriv_debug(lpsi3_c.data(), tmp1.data(), tmp2.data());

	// Comparing two subroutines from the Poisson_solver class

	{
		double max_err = 0.0, max_err1 = 0.0, max_err2 = 0.0;
		for (size_t i = 0; i < g.nrad; i++) {
			max_err = std::max(max_err, std::abs(lpsi1_c[i] - psi1[i]));
			max_err1 = std::max(max_err, std::abs(lpsi2_c[i] - psi2[i]));
			max_err2 = std::max(max_err, std::abs(lpsi3_c[i] - psi3[i]));
		}

		// Print summary 
		std::cout << std::scientific << max_err << " " << max_err1 << " " << max_err2 << std::endl;

	}

	// Comparing two fortran subroutines 

}

void Poisson_solver::test_second_deriv3() {

	std::vector<double> psi1(g.nrad, 0.0), psi2(g.nrad, 0.0), psi3(g.nrad, 0.0);
	std::vector<double> tmp1(g.nrad, 0.0), tmp2(g.nrad, 0.0), tmp1_(g.nrad, 0.0), tmp2_(g.nrad, 0.0);

	std::copy(g.r.begin(), g.r.end(), psi1.begin());
	std::copy(g.r.begin(), g.r.end(), psi2.begin());
	std::copy(g.r.begin(), g.r.end(), psi3.begin());

	std::transform(psi1.begin(), psi1.end(), psi1.begin(), psi_1s_H);
	std::transform(psi2.begin(), psi2.end(), psi2.begin(), psi_2s_H);
	std::transform(psi3.begin(), psi3.end(), psi3.begin(), psi_3s_H);

	std::vector<double> lpsi1_f(g.nrad), lpsi2_f(g.nrad), lpsi3_f(g.nrad);

	std::copy(psi1.begin(), psi1.end(), lpsi1_f.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2_f.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3_f.begin());

	int iatom = 2, nrad = g.nrad;

	second_deriv2_(psi1.data(), &iatom, &nrad);
	second_deriv2_(psi2.data(), &iatom, &nrad);
	second_deriv2_(psi3.data(), &iatom, &nrad);
	second_deriv2_debug_(lpsi1_f.data(), &iatom, &nrad, tmp1.data(), tmp2.data());
	second_deriv2_debug_(lpsi2_f.data(), &iatom, &nrad, tmp1.data(), tmp2.data());
	second_deriv2_debug_(lpsi3_f.data(), &iatom, &nrad, tmp1.data(), tmp2.data());

	// Comparing two subroutines from the Poisson_solver class

	{
		double max_err = 0.0, max_err1 = 0.0, max_err2 = 0.0;
		for (size_t i = 0; i < g.nrad; i++) {
			max_err = std::max(max_err, std::abs(lpsi1_f[i] - psi1[i]));
			max_err1 = std::max(max_err, std::abs(lpsi2_f[i] - psi2[i]));
			max_err2 = std::max(max_err, std::abs(lpsi3_f[i] - psi3[i]));
		}

		// Print summary 
		std::cout << std::scientific << max_err << " " << max_err1 << " " << max_err2 << std::endl;

	}

}
