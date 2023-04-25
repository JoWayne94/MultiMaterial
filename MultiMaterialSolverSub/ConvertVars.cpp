/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : ConvertVars.cpp                            */
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
/*   Description : Convert primitive variables to conservative        */
/*                 variables and vice versa. Convert conservative     */
/*                 variables to Euler flux.                           */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "ConvertVars.h"

using namespace std;

/**
 * @brief Constructor of the ConvertVars class
 */
ConvertVars::ConvertVars() {
}

/**
 * @brief Destructor of the ConvertVars class
 */
ConvertVars::~ConvertVars() {
}

/**
 * @brief Set variables
 *
 * @param   nvar     Number of conservative variables
 * @param   gamma    (Constant) adiabatic index (ratio of specific heats)
 */
void ConvertVars::SetVariables(unsigned int nvar, double gamma, double pinf) {

    nVar = nvar;
    Gamma = gamma;
    Pinf = pinf;

    equationofStates.SetVariables(Gamma, Pinf);

}

/**
 * @brief  Conservative to primitive routine
 *
 * @return  Primitive variables vector
 */
A4Double ConvertVars::conservativeToprimitive(const A4Double& U) {

    double rho = U[0];
    double u_x = U[1]/U[0];
    double u_y = U[2]/U[0];
    double E   = U[3];

    double kin = (double) 0.5*rho*(u_x*u_x + u_y*u_y);
    double e   = (double) (E - kin)/rho;
    double p   = equationofStates.computePressureFromEoS((double) rho, (double) e);

    w[0] = rho;
    w[1] = u_x;
    w[2] = u_y;
    w[3] = p;

    return w;
}

/**
 * @brief  Primitive to conservative routine
 *
 * @return  Conservative variables vector
 */
A4Double ConvertVars::primitiveToconservative(const A4Double& W) {

    double rho = W[0];
    double u_x = W[1];
    double u_y = W[2];
    double p   = W[3];

    double kin = (double) 0.5*rho*(u_x*u_x + u_y*u_y);
    double e   = equationofStates.computeInternalEnergyFromEoS((double) rho, (double) p);
    double E   = (double) rho*e+kin;

    Q[0] = rho;
    Q[1] = rho*u_x;
    Q[2] = rho*u_y;
    Q[3] = E;

    return Q;
}

/**
 * @brief  Conservative variables to flux vector
 *
 * @return  MHD flux vector
 */
A4Double ConvertVars::EulerFlux(const A4Double& U, int dim) {

    Wtemp = conservativeToprimitive(U);
    double rho = Wtemp[0];
    double u_x = Wtemp[1];
    double u_y = Wtemp[2];
    double p   = Wtemp[3];
    double E   = U[3];

    if (dim == 1) {

        F[0] = rho*u_x;
        F[1] = rho*u_x*u_x+p;
        F[2] = rho*u_x*u_y;
        F[3] = (E+p)*u_x;

    } else if (dim == 2) {

        F[0] = rho*u_y;
        F[1] = rho*u_y*u_x;
        F[2] = rho*u_y*u_y+p;
        F[3] = (E+p)*u_y;

    }

    return F;
}