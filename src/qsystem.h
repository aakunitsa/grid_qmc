#ifndef QSYSTEM_H
#define QSYSTEM_H

#include <tuple>
#include <vector>
#include <assert.h>
#include "qparams_reader.h"
//#include <cblas.h>
#include <gsl/gsl_cblas.h>

using namespace std;

class QSystem {
    public:
        QSystem(Params_reader &p) : electrons(p.params["electrons"]) {} 
        QSystem(const QSystem &q) { electrons = q.electrons; }
        virtual double Psi_T(double *x) { return 0.0; } // Trial function
        virtual tuple<double, double, double> E_local(vector<double> &x) { return make_tuple(0.0, 0.0, 0.0); }

    public:
        size_t electrons;
        double alpha; // Variational parameter


};

template<int Ch, int Mult>
class Atom : public QSystem {
    public:
        Atom(Params_reader &p) : QSystem(p) {
            assert (Ch == electrons); // For now
            nalpha = static_cast<size_t>(Ch + Mult - 1)/2;
            nbeta = static_cast<size_t>(Ch - Mult + 1)/2;
            alpha = 0.0;
        }
        /*
        // Requires a call of the base class constructor but how do I ensure consistency?
        Atom(const Atom &a) {
            electrons = a.electrons;
            nalpha = a.nalpha;
            nbeta = a.nbeta;
            alpha = a.alpha;
        }
        */
    public:
        size_t nalpha, nbeta;
};


template<>
class Atom<2, 1>: public QSystem {
    public:
        Atom(Params_reader &p) : QSystem(p) {
            assert (electrons == 2);
            nalpha = 1;
            nbeta = 1;
            alpha = 0.4;
        }
        /*
        Atom(const Atom &a) {
            electrons = a.electrons;
            nalpha = a.nalpha;
            nbeta = a.nbeta;
        }
        */
        double Psi_T(double *xx) {

            const size_t DIM_1e = 3;
            double *r1 = xx, *r2 = xx + DIM_1e;
            double r12[DIM_1e];
            double r1_abs, r2_abs, r12_abs;
            cblas_dcopy(DIM_1e, r1, 1, r12, 1);
            cblas_daxpy(DIM_1e, -1.0, r2, 1, r12, 1);
            r12_abs = cblas_dnrm2(DIM_1e, r12, 1);
            r1_abs = cblas_dnrm2(DIM_1e, r1, 1);
            r2_abs = cblas_dnrm2(DIM_1e, r2, 1);

            return exp(-2.*(r1_abs+r2_abs)) * exp(r12_abs*0.5/(1.+alpha*r12_abs));

        }

        tuple<double, double, double> E_local(vector<double> &xx) {

            double E = 0.0, E_dphidalpha = 0.0, dphidalpha = 0.0;
            auto xx_ptr = xx.data();
            constexpr size_t DIM_1e = 3;
            size_t Nw = xx.size() / electrons / DIM_1e;

            for (size_t i = 0; i < Nw; i++) {

                double *r1 = xx_ptr + i*3*electrons, *r2 = xx_ptr + i*3*electrons + DIM_1e;
                double r12[DIM_1e];
                double r1_abs, r2_abs, r12_abs;
                //std::transforms(r1.begin(), r1.end(), r1.bein(), r12.begin(), std::minus<>);
                cblas_dcopy(DIM_1e, r1, 1, r12, 1);
                cblas_daxpy(DIM_1e, -1.0, r2, 1, r12, 1);
                r12_abs = cblas_dnrm2(DIM_1e, r12, 1);
                r1_abs = cblas_dnrm2(DIM_1e, r1, 1);
                r2_abs = cblas_dnrm2(DIM_1e, r2, 1);

                //printf("r1_abs = %18.10f r2_abs = %18.10f r12_abs = %18.10f\n", r1_abs, r2_abs, r12_abs);
                //printf("%18.10f %18.10f %18.10f\n", r1_abs, r2_abs, r12_abs);

                double fac = 1. / (1. + alpha * r12_abs) ;
                double fac2 =  cblas_ddot(DIM_1e, r1, 1, r12, 1)/r1_abs -  cblas_ddot(DIM_1e, r2, 1, r12, 1)/r2_abs;
                double E_int = -4. + 1./r12_abs * (fac2*fac*fac - fac*fac*fac + 1.) - 0.25 * fac * fac * fac * fac, 
                    dphidalpha_int = -0.5*r12_abs*r12_abs * fac * fac;
                E += E_int;
                dphidalpha += dphidalpha_int;
                E_dphidalpha += E_int * dphidalpha_int;

            }

            return make_tuple(E / Nw, E_dphidalpha / Nw,  dphidalpha / Nw );
        }

    public:
        size_t nalpha, nbeta;
};

template<>
class Atom<2, 3> : public QSystem {
    public:
        Atom(Params_reader &p) : QSystem(p) {
            assert (electrons == 2);
            nalpha = 2;
            nbeta = 0;
            alpha = 0.2;
        }
        /*
        Atom(const Atom &a) {
            electrons = a.electrons;
            nalpha = a.nalpha;
            nbeta = a.nbeta;
        }
        */
        double Psi_T(double *xx) {

            const size_t DIM_1e = 3;
            double *r1 = xx, *r2 = xx + DIM_1e;
            double r12[DIM_1e];
            double r1_abs, r2_abs, r12_abs;
            cblas_dcopy(DIM_1e, r1, 1, r12, 1);
            cblas_daxpy(DIM_1e, -1.0, r2, 1, r12, 1);
            r12_abs = cblas_dnrm2(DIM_1e, r12, 1);
            r1_abs = cblas_dnrm2(DIM_1e, r1, 1);
            r2_abs = cblas_dnrm2(DIM_1e, r2, 1);

            //return (exp(-2.*(r1_abs+ 0.5 * r2_abs)) - exp(-2.*(r2_abs+ 0.5 * r1_abs))) * exp(r12_abs*0.25/(1.+alpha*r12_abs));
            return (exp(-2.*(r1_abs+ 0.5 * r2_abs)) * (2.0 - 2.0 * r2_abs) - 
                    exp(-2.*(r2_abs+ 0.5 * r1_abs)) *  (2.0 - 2.0 * r1_abs)) * exp(r12_abs*0.25/(1.+alpha*r12_abs));

        }

        tuple<double, double, double> E_local(vector<double> &xx) {

            double E = 0.0, E_dphidalpha = 0.0, dphidalpha = 0.0;
            auto xx_ptr = xx.data();
            constexpr size_t DIM_1e = 3;
            size_t Nw = xx.size() / electrons / DIM_1e;

            for (size_t i = 0; i < Nw; i++) {

                double *r1 = xx_ptr + i*3*electrons, *r2 = xx_ptr + i*3*electrons + DIM_1e;
                double r12[DIM_1e];
                double r1_abs, r2_abs, r12_abs;
                //std::transforms(r1.begin(), r1.end(), r1.bein(), r12.begin(), std::minus<>);
                cblas_dcopy(DIM_1e, r1, 1, r12, 1);
                cblas_daxpy(DIM_1e, -1.0, r2, 1, r12, 1);
                r12_abs = cblas_dnrm2(DIM_1e, r12, 1);
                r1_abs = cblas_dnrm2(DIM_1e, r1, 1);
                r2_abs = cblas_dnrm2(DIM_1e, r2, 1);

                //printf("r1_abs = %18.10f r2_abs = %18.10f r12_abs = %18.10f\n", r1_abs, r2_abs, r12_abs);
                //printf("%18.10f %18.10f %18.10f\n", r1_abs, r2_abs, r12_abs);

                double invD[4], LD[4], ND[12]; // 4 (*3) - for two electrons
                //D[0] = exp(-2. * r1_abs), D[1] = exp(-2. * r2_abs);
                //D[2] = exp(-1. * r1_abs), D[3] = exp(-1. * r2_abs);
                double d = exp(-2.*(r1_abs+ 0.5 * r2_abs)) * (2.0 - 2.0 * r2_abs) - 
                    exp(-2.*(r2_abs+ 0.5 * r1_abs)) *  (2.0 - 2.0 * r1_abs);
                invD[0] = exp(-1. * r2_abs) *  (2.0 - 2.0 * r2_abs)/ d; invD[1] = -1. * exp(-2. * r2_abs) / d;
                invD[2] = -1. * exp(-1. * r1_abs) *  (2.0 - 2.0 * r1_abs)/ d; invD[3] = exp(-2. * r1_abs) / d;
                //double d  = exp(-2. * r1_abs - r2_abs) - exp(-2. * r2_abs - r1_abs);
                //invD[0] = exp(-1. * r2_abs) / d; invD[1] = -1. * exp(-2. * r2_abs) / d;
                //invD[2] = -1. * exp(-1. * r1_abs) / d; invD[3] = exp(-2. * r1_abs) / d;

                //printf("invD:\n");
                //printf("%13.6f %13.6f\n%13.6f %13.6f\n", invD[0], invD[1], invD[2], invD[3]);
                assert (d != 0.0 );
                //assert (invD[1] < 0.0);
                //assert (invD[2] < 0.0);

                double D[4], unit[4], buff[4], buff1[6];
                //D[0] = exp(-2. * r1_abs); D[1] = exp(-2. * r2_abs);
                //D[2] = exp(-1. * r1_abs); D[3] = exp(-1. * r2_abs);
                D[0] = exp(-2. * r1_abs); D[1] = exp(-2. * r2_abs);
                D[2] = exp(-1. * r1_abs) * (2.0 - 2.0 * r1_abs); D[3] = exp(-1. * r2_abs) * (2.0 - 2.0 * r2_abs);
                /*
                printf("%13.6f %13.6f\n%13.6f %13.6f\n", D[0] * invD[0] + D[1] * invD[2], 
                                                         D[0] * invD[1] + D[1] * invD[3], 
                                                         D[2] * invD[0] + D[3] * invD[2], 
                                                         D[2] * invD[1] + D[3] * invD[3]);
                */
                // Just a quick test
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0, D, 2, invD, 2, 0.0, unit, 2); 
                assert (d != 0.0);
                //assert (invD[1] < 0 && invD[2] < 0);
               /* 
                printf("Checking if the inverse matrix is correct\n");
                printf("invD:\n");
                printf("%13.6f %13.6f\n%13.6f %13.6f\n", invD[0], invD[1], invD[2], invD[3]);
                printf("D:\n");
                printf("%13.6f %13.6f\n%13.6f %13.6f\n", D[0], D[1], D[2], D[3]);
                printf("D * invD:\n");
                printf("%13.6f %13.6f\n%13.6f %13.6f\n", unit[0], unit[1], unit[2], unit[3]);
                printf("%13.6f %13.6f\n%13.6f %13.6f\n", D[0] * invD[0] + D[1] * invD[2], 
                                                         D[0] * invD[1] + D[1] * invD[3], 
                                                         D[2] * invD[0] + D[3] * invD[2], 
                                                         D[2] * invD[1] + D[3] * invD[3]);
                */
                
                //LD[0] = (4. - 4./r1_abs) * exp(-2. * r1_abs), LD[1] = (4. - 4./r2_abs) * exp(-2. * r2_abs);
                //LD[2] = (1. - 2./r1_abs) * exp(-1. * r1_abs), LD[3] = (1. - 2./r2_abs) * exp(-1. * r2_abs);
                LD[0] = (4. - 4./r1_abs) * exp(-2. * r1_abs), LD[1] = (4. - 4./r2_abs) * exp(-2. * r2_abs);
                LD[2] = ( 10. - 8./r1_abs - 2. * r1_abs) * exp(-1. * r1_abs), LD[3] = (10. - 8./r2_abs - 2. * r2_abs) * exp(-1. * r2_abs);

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0,  invD, 2, LD, 2, 0.0, buff, 2); 

                // ND should be stored in a column major order and processed in a loop via ddot? Or should it be matrix-matrix multiply?
/*
                ND[0] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[0], ND[1] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[0];
                ND[2] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[0], ND[3] =  - 1./r2_abs * exp(-1. * r2_abs) * r2[0];
                ND[4] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[1], ND[5] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[1];
                ND[6] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[1], ND[7] =  - 1./r2_abs * exp(-1. * r2_abs) * r2[1];
                ND[8] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[2], ND[9] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[2];
                ND[10] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[2], ND[11] =  - 1./r2_abs * exp(-1. * r2_abs) * r2[2];
*/
                /*
                ND[0] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[0], ND[1] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[0],
                ND[2] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[1], ND[3] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[1],
                ND[4] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[2], ND[5] =  - 1./r1_abs * exp(-1. * r1_abs) * r1[2],
                ND[6] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[0], ND[7] =  - 1./r2_abs * exp(-1. * r2_abs) * r2[0],
                ND[8] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[1], ND[9] =  - 1./r2_abs * exp(-1. * r2_abs) * r2[1],
                ND[10] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[2], ND[11] =- 1./r2_abs * exp(-1. * r2_abs) * r2[2];
                */
                ND[0] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[0], ND[1] =  (-4.0 + 2.0 * r1_abs) * 1./r1_abs * exp(-1. * r1_abs) * r1[0],
                ND[2] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[1], ND[3] =  (-4.0 + 2.0 * r1_abs) * 1./r1_abs * exp(-1. * r1_abs) * r1[1],
                ND[4] =  - 2./r1_abs * exp(-2. * r1_abs) * r1[2], ND[5] =  (-4.0 + 2.0 * r1_abs) * 1./r1_abs * exp(-1. * r1_abs) * r1[2],
                ND[6] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[0], ND[7] =   (-4.0 + 2.0 * r2_abs) *1./r2_abs * exp(-1. * r2_abs) * r2[0],
                ND[8] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[1], ND[9] =  (-4.0 + 2.0 * r2_abs) * 1./r2_abs * exp(-1. * r2_abs) * r2[1],
                ND[10] =  - 2./r2_abs * exp(-2. * r2_abs) * r2[2], ND[11] = (-4.0 + 2.0 * r2_abs) *1./r2_abs * exp(-1. * r2_abs) * r2[2];

                buff1[0] = cblas_ddot(2, invD, 1, ND, 1);
                buff1[1] = cblas_ddot(2, invD, 1, ND + 2, 1);
                buff1[2] = cblas_ddot(2, invD, 1, ND + 4, 1);
                buff1[3] = cblas_ddot(2, invD + 2, 1, ND + 6, 1);
                buff1[4] = cblas_ddot(2, invD + 2, 1, ND + 8, 1);
                buff1[5] = cblas_ddot(2, invD + 2, 1, ND + 10, 1);

                // Taking care of Jastrow factor here
                double fac = 1. / (1. + alpha * r12_abs);
                double LJ = 1./16. * fac * fac * fac * fac + 0.5 / r12_abs * fac * fac * fac; 
                double NJ[] {0.25 * r12[0] / r12_abs * fac * fac, 0.25 * r12[1] / r12_abs * fac * fac, 0.25 * r12[2] / r12_abs * fac * fac,   
                                   -0.25 * r12[0] / r12_abs * fac * fac, -0.25 * r12[1] / r12_abs * fac * fac, -0.25 * r12[2] / r12_abs * fac * fac};

                // E_L = -LJ - 0.5 * \sum_{k = 1}^{N_e} L_k D / D - \sum_{k = 1}^{N_e} N_k D \dot N_k J / (J * D)
                // General formulas can be found in W. Schattke & R. Diez Muino, p.77
                double E_int = -LJ - cblas_ddot(6, buff1, 1, NJ, 1) - 0.5 * ( buff[0] + buff[3] ) -2.0 * (1./ r1_abs + 1./r2_abs) + 1/r12_abs, 
                       dphidalpha_int = -0.25*r12_abs*r12_abs * fac * fac;

                E += E_int;
                dphidalpha += dphidalpha_int;
                E_dphidalpha += E_int * dphidalpha_int;

            }

            return make_tuple(E / Nw, E_dphidalpha / Nw,  dphidalpha / Nw );
        }

    public:
        size_t nalpha, nbeta;
};

#endif
