/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : NumericalFlux.cpp                          */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 22/01/2022                                                */
/*                                                                    */
/*   Programmer :                                      */
/*                                                                    */
/*   Description : Various centred and Riemann schemes to             */
/*                 compute the numerical fluxes at i+1/2 and i-1/2    */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "NumericalFlux.h"

using namespace std;

/**
 * @brief Constructor of the NumericalFlux class
 */
NumericalFlux::NumericalFlux() {
}

/**
 * @brief Destructor of the NumericalFlux class
 */
NumericalFlux::~NumericalFlux() {
}

/**
 * @brief Set variables
 *
 * @param   nvar     Number of conservative variables
 * @param   gamma    (Constant) adiabatic index (ratio of specific heats)
 */
void NumericalFlux::SetVariables(const A4Double& Ql, const A4Double& Qr, unsigned int nvar, double gamma, double pinf) {

    nVar = nvar;
    Gamma = gamma;
    Pinf = pinf;

    convertVars.SetVariables(nVar, Gamma, Pinf);
    equationofStates.SetVariables(Gamma, Pinf);

    wl = convertVars.conservativeToprimitive((A4Double) Ql);
    wr = convertVars.conservativeToprimitive((A4Double) Qr);

    rhol = wl[0];
    rhor = wr[0];
    uxl = wl[1];
    uxr = wr[1];
    uyl = wl[2];
    uyr = wr[2];
    pl = wl[3];
    pr = wr[3];

    /** Ideal Gas EoS **/
    al = equationofStates.computeSoundSpeedFromEoS((double) rhol, (double) pl);
    ar = equationofStates.computeSoundSpeedFromEoS((double) rhor, (double) pr);
}

/**
 * @brief Estimate wave speeds S_L and S_R in x and y direction
 *
 *        Toro's book pg 327 section 10.5
 *
 * @param   Ql      Left vector of conservative variables
 * @param   Qr      Right vector of conservative variables
 * @param   Sl      Left wave speed
 * @param   Sr      Right wave speed
 */
void NumericalFlux::WaveEstimates(double& Sleft, double& Sright, int dim) const {

    /// More restrictive wave speed definitions
    // Sleft = -max(sqrt(uxl*uxl + uyl*uyl) + al, sqrt(uxr*uxr + uyr*uyr) + ar);
    // Sright = max(sqrt(uxl*uxl + uyl*uyl) + al, sqrt(uxr*uxr + uyr*uyr) + ar);

    if (dim == 1) {
        Sright = std::max(fabs(uxl) + al, fabs(uxr) + ar); // uxl - al;
        Sleft = -Sright; // uxr + ar;
    } else if (dim == 2) {
        Sright = std::max(fabs(uyl) + al, fabs(uyr) + ar); // uyl - al;
        Sleft = -Sright; // Sright = uyr + ar;
    }

    // Add here different Wave speed estimates ...
    // delete[] wl;
}

/**
 * @brief HLL Approx Riemann Solver
 *
 *        Toro's book pg 320 section 10.3
 *
 * @param   Ql        Left vector of conservative variables
 * @param   Qr        Right vector of conservative variables
 * @param   Fhll      HLL intercell flux for the approximate Godunov method
 */
void NumericalFlux::HLLFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fhll, int dim) {

    WaveEstimates(Sl, Sr, dim); // calculate Sl and Sr

    A4Double Fl, Fr;

    Fl = convertVars.EulerFlux((A4Double) Ql, dim);
    Fr = convertVars.EulerFlux((A4Double) Qr, dim);

    /** Eq (10.21) **/
    if (Sl >= 0) {
        Fhll = Fl;
        return;
    }

    if (Sr <= 0) {
        Fhll = Fr;
        return;
    }

    for (int n = 0; n < nVar; n++) {
        Fhll[n] = (Sr * Fl[n] - Sl * Fr[n] + Sl * Sr * (Qr[n] - Ql[n]))/(Sr - Sl);
    }
}

/**
 * @brief HLLC Approx Riemann Solver
 *
 *        Toro's book pg 322 section 10.4
 *
 * @param   Ql         Left vector of conservative variables
 * @param   Qr         Right vector of conservative variables
 * @param   Fhllc      HLLC intercell flux for the approximate Godunov method
 */
void NumericalFlux::HLLCFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fhllc, int dim) {

    double Sstar, pLR;

    WaveEstimates(Sl, Sr, dim); // calculate Sl and Sr

    A4Double Fl, Fr, Dstar, Qstarl, Qstarr, Fstarl, Fstarr;

    Fl = convertVars.EulerFlux((A4Double) Ql, dim);
    Fr = convertVars.EulerFlux((A4Double) Qr, dim);

    /** Eq (10.26) **/
    if (Sl >= 0) {
        Fhllc = Fl;
        return;
    }

    if (Sr <= 0) {
        Fhllc = Fr;
        return;
    }

    if (dim == 1) {

        Sstar = (pr - pl + rhol * uxl * (Sl - uxl) - rhor * uxr * (Sr - uxr)) / (rhol * (Sl - uxl) - rhor * (Sr - uxr)); // eq (10.37)

        /** Variant 2 **/
        pLR = 0.5 * (pl + pr + rhol * (Sl - uxl) * (Sstar - uxl) + rhor * (Sr - uxr) * (Sstar - uxr)); // eq (10.42)
        Dstar = {0, 1, 0, Sstar}; // 1D, eq (10.40)

    } else if (dim == 2) {

        Sstar = (pr - pl + rhol * uyl * (Sl - uyl) - rhor * uyr * (Sr - uyr)) / (rhol * (Sl - uyl) - rhor * (Sr - uyr)); // eq (10.37)

        /** Variant 2 **/
        pLR = 0.5 * (pl + pr + rhol * (Sl - uyl) * (Sstar - uyl) + rhor * (Sr - uyr) * (Sstar - uyr)); // eq (10.42)
        Dstar = {0, 0, 1, Sstar}; // 1D, eq (10.40)

    }

    /** Eq (10.43) **/
    for (int n = 0; n < nVar; n++) {
        Qstarl[n] = (Sl * Ql[n] - Fl[n] + pLR * Dstar[n]) / (Sl - Sstar);
        Qstarr[n] = (Sr * Qr[n] - Fr[n] + pLR * Dstar[n]) / (Sr - Sstar);
    }
    /** Eq (10.44) **/
    for (int n = 0; n < nVar; n++) {
        Fstarl[n] = (Sstar * (Sl * Ql[n] - Fl[n]) + Sl * pLR * Dstar[n]) / (Sl - Sstar);
        Fstarr[n] = (Sstar * (Sr * Qr[n] - Fr[n]) + Sr * pLR * Dstar[n]) / (Sr - Sstar);
    }

    if (Sl <=0 && Sstar >= 0) {
        Fhllc = Fstarl;
        return;
    }

    if (Sstar <=0 && Sr >= 0) {
        Fhllc = Fstarr;
        return;
    }
}

/**
 * @brief Force flux
 *
 *        Toro's book pg 600 section 18.2
 *
 * @param   Ql          Left vector of conservative variables
 * @param   Qr          Right vector of conservative variables
 * @param   Fforce      FORCE intercell flux
 * @param   dx          Spatial discretisation length in x-direction
 * @param   dt          Time step
 */
void NumericalFlux::FORCEFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fforce, const double& dx, const double& dt, int dim) {

    A4Double Fl, Fr, Flf, Qtemp, Flw;

    Fl = convertVars.EulerFlux((A4Double) Ql, dim);
    Fr = convertVars.EulerFlux((A4Double) Qr, dim);

    for (int n = 0; n < nVar; n++) {
        Flf[n] = 0.5 * (Fl[n] + Fr[n]) - 0.5 * (dx/dt) * (Qr[n] - Ql[n]);   // Lax–Friedrichs flux, eq (18.6)
        Qtemp[n] = 0.5 * (Ql[n] + Qr[n]) - 0.5 * (dt/dx) * (Fr[n] - Fl[n]); // Two–step Lax–Wendroff flux, eq (18.7)
    }

    Flw = convertVars.EulerFlux(Qtemp, dim); // Two-step Lax–Wendroff = Richtmyer flux, eq (18.7)

    for (int n = 0; n < nVar; n++) {
        Fforce[n] = 0.5 * (Flf[n] + Flw[n]); // FORCE flux, eq (18.8)
    }
}