/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : ExactRPSolver.cpp                          */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 23/01/2022                                                */
/*                                                                    */
/*   Programmer :                                         */
/*                                                                    */
/*   Description : Riemann Problem exact solution routines specific   */
/*                 for ideal and stiffened gases.                     */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "ExactRPSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

/**
 * @brief Constructor of the ExactRPSolver class
 */
ExactRPSolver::ExactRPSolver() {
}

/**
 * @brief Destructor of the ExactRPSolver class
 */
ExactRPSolver::~ExactRPSolver() {
}

/**
 * @brief Set variables
 *
 * @param   gamma   (Constant) adiabatic index (ratio of specific heats)
 * @param   nvar    Number of conservative variables
 */
void ExactRPSolver::SetVariables(const A4Double& wl, const A4Double& wr, double gammaL, double gammaR, double pinfL, double pinfR) {

    GammaL = gammaL;
    GammaR = gammaR;
    PinfL = pinfL;
    PinfR = pinfR;

    equationofStatesL.SetVariables(GammaL, PinfL);
    equationofStatesR.SetVariables(GammaR, PinfR);

    G1 = (GammaL-1.0)/(2.0*GammaL);
    G2 = (GammaL+1.0)/(2.0*GammaL);
    G3 = 2.0*GammaL/(GammaL-1.0);
    G4 = 2.0/(GammaL-1.0);
    G5 = 2.0/(GammaL+1.0);
    G6 = (GammaL-1.0)/(GammaL+1.0);
    G7 = (GammaL-1.0)/2.0;

    G9 = (GammaR-1.0)/(2.0*GammaR);
    G10 = (GammaR+1.0)/(2.0*GammaR);
    G11 = 2.0*GammaR/(GammaR-1.0);
    G12 = 2.0/(GammaR-1.0);
    G13 = 2.0/(GammaR+1.0);
    G14 = (GammaR-1.0)/(GammaR+1.0);
    G15 = (GammaR-1.0)/2.0;

    rhol = wl[0];
    rhor = wr[0];
    ul = wl[1];
    ur = wr[1];
    pl = wl[3];
    pr = wr[3];
    al = equationofStatesL.computeSoundSpeedFromEoS(rhol, pl);
    ar = equationofStatesR.computeSoundSpeedFromEoS(rhor, pr);

}

double ExactRPSolver::ComputePStar() {

    // Newton-Raphson method
    int iter = 0;
    int iterMax = 100;

    double Err = 1.0;
    double Tol = 1.0E-6;

    double Pold = GuessPressure();

    double Piter;

    while (Err > Tol && iter < iterMax) {

        Piter = Pold - (f_NR(Pold)/df_NR(Pold));

        Err = std::max(2.0 * std::abs((Piter-Pold)/(Piter+Pold)), std::abs(Piter - Pold));

        Pold = std::max(Tol, Piter);

        iter += 1;
    }

    // std::cout << "No. of iterations: " << iter << std::endl;

    return Pold;
}

double ExactRPSolver::ComputeUStar(const double& pstar) {

    return 0.5*(ul+ur)+0.5*(fr_NR(pstar)-fl_NR(pstar));

}

void ExactRPSolver::FindPattern(const double& pstar, int& Pattern) const {
    if (pstar < pl && pstar < pr) {
        // The pattern is rarefaction-contact-rarefaction
        Pattern = 1;
    } else if (pstar < pl && pstar >= pr) {
        // The pattern is rarefaction-contact-shock
        Pattern = 2;
    } else if (pstar >= pl && pstar < pr) {
        // The pattern is shock-contact-rarefaction
        Pattern = 3;
    } else if (pstar >= pl && pstar >= pr) {
        // The pattern is shock-contact-shock
        Pattern = 4;
    }
}

void ExactRPSolver::ComputeExactSolution(const double& time, const VDouble& x, const double& xd, VectA4& uExact) {

    // Check that vacuum is not generated
    double DeltaU = ur - ul;

    bool Vacuum1 = G4*(al+ar) <= DeltaU;
    bool Vacuum2 = G12*(al+ar) <= DeltaU;
    if (Vacuum1 || Vacuum2) {
        std::cerr << "The initial data is such that vacuum is generated. Program stopped.";
        std::_Exit(EXIT_FAILURE);
    }

    double Pstar = ComputePStar();
    double ustar = ComputeUStar(Pstar);

    // Wave pattern
    int Pattern;
    FindPattern(Pstar, Pattern);

    // std::cout<< "Pattern = " << Pattern << " pstar = " << Pstar << " ustar = " << ustar << std::endl;

    // Sample the solution
    for (int i = 0; i < x.size(); i++) {

        double csi = (x[i] - xd)/time;

        A4Double wExact;
        wExact = Sample(Pstar, ustar, Pattern, csi);

        uExact(i, 0) = wExact;
    }
}

double ExactRPSolver::GuessPressure() const {

    double PGuess;

    // Guess pressure based on PVRS RP Solver
    int Quser = 2;
    double cup, ppv, pmin, pmax, qmax;

    cup = 0.25*(rhol+rhor)*(al+ar);
    ppv = std::max(0.0,0.5*(pl+pr)+0.5*(ul-ur)*cup);

    pmin = std::min(pl,pr);
    pmax = std::max(pl,pr);

    qmax = pmax/pmin;

    bool Cond1 = (qmax <= Quser) && ((pmin <= ppv) && (pmax >= ppv));

    if (Cond1) {
        /// Select PVRS RP Solver
        PGuess = ppv;
        // std::cout << "PVRS selected." << std::endl;
    } else {
        if (ppv < pmin) {
            /// Select two-rarefaction RP Solver
//            double pq,um,ptl,ptr;
//
//            pq = pow(pl/pr,G1);
//            um = (pq*ul/al+ur/ar+G4*(pq-1.0))/(pq/al+1.0/ar);
//
//            ptl = 1.0+G7*(ul-um)/al;
//            ptr = 1.0+G7*(um-ur)/ar;
//
//            PGuess = 0.5*(pl*pow(ptl,G3)+pr*pow(ptr,G3));
            PGuess = 0.5 * (pl + pr);

            // std::cout << "TRRS selected." << std::endl;

        } else {
            /// Select two-shock RP solver with PVRS as estimate
            double gel,ger;

            gel = sqrt((G5/rhol)/(G6*pl+ppv));
            ger = sqrt((G13/rhor)/(G14*pr+ppv));

            PGuess = (gel*pl+ger*pr-(ur-ul))/(gel+ger);

            // std::cout << "TSRS selected." << std::endl;
        }
    }

    return PGuess;
}

double ExactRPSolver::fl_NR(const double& p) const {

    double fl;

    double Al = G5/rhol;
    /// Stiffened gas
    double Bl = G6 * pl + (2.0 * GammaL * PinfL)/(GammaL + 1.0);

    // Newton-Raphson function
    fl = (p > pl) ? (p-pl)*sqrt(Al/(p+Bl)) : G4*al*(pow((p + PinfL)/(pl + PinfL),G1)-1.0);

    return fl;
}

double ExactRPSolver::fr_NR(const double& p) const {

    double fr;

    double Ar = G13/rhor;
    /// Stiffened gas
    double Br = G14 * pr + (2.0 * GammaR * PinfR)/(GammaR + 1.0);

    // Newton-Raphson function
    fr = (p > pr) ? (p-pr)*sqrt(Ar/(p+Br)) : G12*ar*(pow((p + PinfR)/(pr + PinfR),G9)-1.0);

    return fr;
}

double ExactRPSolver::f_NR(const double& p) const {

    double f;

    double fl = fl_NR(p);
    double fr = fr_NR(p);

    // Newton-Raphson function
    f = fl + fr + (ur - ul);

    return f;
}

double ExactRPSolver::df_NR(const double& p) const {

    double df;

    double Al = G5/rhol;
    double Ar = G13/rhor;

    /// Stiffened gas
    double Bl = G6 * pl + (2.0 * GammaL * PinfL)/(GammaL + 1.0);
    double Br = G14 * pr + (2.0 * GammaR * PinfR)/(GammaR + 1.0);

    // Newton-Raphson function
    double dfL = (p > pl) ? sqrt(Al/(Bl+p))*(1.0-0.5*(p-pl)/(Bl+p))
                          : 1.0/(rhol*al)*pow(((p + PinfL)/(pl + PinfL)),-G2);
    double dfR = (p > pr) ? sqrt(Ar/(Br+p))*(1.0-0.5*(p-pr)/(Br+p))
                          : 1.0/(rhor*ar)*pow(((p + PinfR)/(pr + PinfR)),-G10);

    df = dfL + dfR;

    return df;
}

void ExactRPSolver::ComputeRhoStar(const double& pstar, const double& ustar, const int& Pattern,
                                   double& rholstar, double& alstar, double& Shl, double& Stl, double& Sl, double& rhorstar, double& arstar, double& Shr, double& Str, double& Sr) const {

    /// Left of contact discontinuity
    if (Pattern == 1 || Pattern == 2) {  //  rcr or rcs
        rholstar = rhol*pow((pstar + PinfL)/(pl + PinfL),1.0/GammaL);
        alstar   = al*pow((pstar + PinfL)/(pl + PinfL),G1);
        Shl      = ul-al;
        Stl      = ustar-alstar;
    } else {  //  scr or scs
        rholstar = rhol*((pstar + PinfL)/(pl + PinfL)+G6)/(1.0+G6*(pstar + PinfL)/(pl + PinfL));
        Sl       = ul-al*sqrt(G1+G2*(pstar + PinfL)/(pl + PinfL));
    }

    /// Right of contact discontinuity
    if (Pattern == 1 || Pattern == 3) {  //  rcr or scr
        rhorstar = rhor*pow((pstar + PinfR)/(pr + PinfR),1.0/GammaR);
        arstar   = ar*pow((pstar + PinfR)/(pr + PinfR),G9);
        Shr      = ur+ar;
        Str      = ustar+arstar;
    } else {  //  rcs or scs
        rhorstar = rhor*((pstar + PinfR)/(pr + PinfR)+G14)/(1.0+G14*(pstar + PinfR)/(pr + PinfR));
        Sr       = ur+ar*sqrt(G9+G10*(pstar + PinfR)/(pr + PinfR));
    }
}

A4Double ExactRPSolver::Sample(const double& pstar, const double& ustar, const int& Pattern, const double& csi) const {

    A4Double Wexact;

    // Sampling the solution
    double rholstar, alstar, Shl, Stl, Sl;
    double rhorstar, arstar, Shr, Str, Sr, C;

    ComputeRhoStar(pstar, ustar, Pattern, rholstar, alstar, Shl, Stl, Sl, rhorstar, arstar, Shr, Str, Sr);

    switch (Pattern) {
        case 1:
        {
            if (csi < Shl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[3] = pl;
            }
            else if (csi >= Shl && csi < Stl)
            {
                C = G5 * (al + G7 * (ul - csi));
                Wexact[0] = rhol*pow(C/al, G4);
                Wexact[1] = G5*(al+G7*ul+csi);
                Wexact[3] = (pl + PinfL)*pow(C/al, G3) - PinfL;
            }
            else if (csi >= Stl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= ustar && csi < Str)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= Str && csi < Shr)
            {
                C = G13 * (ar - G15 * (ur - csi));
                Wexact[0] = rhor*pow(C/ar, G12);
                Wexact[1] = G13*(-ar+G15*ur+csi);
                Wexact[3] = (pr + PinfR)*pow(C/ar, G11) - PinfR;
            }
            else if (csi >= Shr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[3] = pr;
            }

            break;
        }

        case 2:
        {
            if (csi < Shl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[3] = pl;
            }
            else if (csi >= Shl && csi < Stl)
            {
                C = G5 * (al + G7 * (ul - csi));
                Wexact[0] = rhol*pow(C/al, G4);
                Wexact[1] = G5*(al+G7*ul+csi);
                Wexact[3] = (pl + PinfL)*pow(C/al, G3) - PinfL;
            }
            else if (csi >= Stl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= ustar && csi < Sr)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= Sr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[3] = pr;
            }

            break;
        }

        case 3:
        {
            if (csi < Sl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[3] = pl;
            }
            else if (csi >= Sl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= ustar && csi < Str)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= Str && csi < Shr)
            {
                C = G13 * (ar - G15 * (ur - csi));
                Wexact[0] = rhor*pow(C/ar, G12);
                Wexact[1] = G13*(-ar+G15*ur+csi);
                Wexact[2] = (pr + PinfR)*pow(C/ar, G11) - PinfR;
            }
            else if (csi >= Shr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[3] = pr;
            }

            break;
        }

        case 4:
        {
            if (csi < Sl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[3] = pl;
            }
            else if (csi >= Sl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= ustar && csi < Sr)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[3] = pstar;
            }
            else if (csi >= Sr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[3] = pr;
            }

            break;
        }

        default:
        {
            cout << "default \n" << endl;
            break;
        }
    }

    return Wexact;
}