/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : Euler.cpp                                  */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 10/02/2022                                                */
/*                                                                    */
/*   Programmer :                                  */
/*                                                                    */
/*   Description : Driver for the Euler solver                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "Euler.h"
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

/**
 * @brief Constructor for the Euler class
 */
Euler::Euler() {
}

/**
 * @brief Destructor for the Euler class
 *
 * 	      Clean up memory
 */
Euler::~Euler() {
    delete u1;
    delete uBc1;
    delete u2;
    delete uBc2;
    delete wExact;
    delete xCells;
    delete yCells;
    delete phi;
    delete phiBc;
}

/**
 * @brief Set domain size
 *
 *	      Initialize Lx and Ly, xLeftDomain, xRightDomain, yBottomDomain, and yTopDomain variables
 *
 * @param   xlength        Length of the domain in the x-direction
 * @param   ylength        Length of the domain in the y-direction
 * @param   xLeft          Left domain boundary in the x-direction
 * @param   xRight         Right domain boundary in the x-direction
 * @param   yBottom        Bottom domain boundary in the y-direction
 * @param   yTop           Top domain boundary in the y-direction
 */
void Euler::SetDomainSize(double xlength, double ylength, double xLeft, double xRight, double yBottom, double yTop) {
    Lx = xlength;
    Ly = ylength;
    xLeftDomain = xLeft;
    xRightDomain = xRight;
    yBottomDomain = yBottom;
    yTopDomain = yTop;
}

/**
 * @brief Set number of cells
 *
 *  	  Initialize Nx and Ny arrays, and other variables
 *
 * @param   nx             Number of cells in the x-direction
 * @param   ny             Number of cells in the y-direction
 * @param   nxGhost        Number of fictitious cells in the x-direction
 * @param   nyGhost        Number of fictitious cells in the y-direction
 */
void Euler::SetNumCells(unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost) {
    Nx = nx;
    Ny = ny;
    N = Nx * Ny;
    NxGhost = nxGhost;
    NyGhost = nyGhost;
}

/**
 * @brief Set Courant–Friedrichs–Lewy number
 *
 * @param   cfl     Courant–Friedrichs–Lewy number
 */
void Euler::SetCourantNumber(double cfl) {
    Cfl = cfl;
}

/**
 * @brief Set mesh size
 *
 * @param   deltax        Mesh size in the x-direction
 * @param   deltay        Mesh size in the y-direction
 */
void Euler::SetMeshSize(double deltax, double deltay) {
    dx = deltax;
    dy = deltay;
}

/**
 * @brief Set number of conservative variables
 *
 * @param   nvar     Number of conservative variables
 */
void Euler::SetnVar(unsigned int nvar) {
    nVar = nvar;
}

/**
 * @brief Set test number
 *
 * @param   testNum   Choice of pre-defined test cases
 */
void Euler::SetTestNum(unsigned int testNum) {
    testNumber = testNum;
}

/**
 * @brief Set output file name
 *
 * @param   Name   Data .dat file name
 */
void Euler::SetFileName(std::string Name) {
    name = std::move(Name);
}

/**
 * @brief Set BC types
 *
 * @param   xbc   BC type in x-direction
 * @param   ybc   BC type in y-direction
 */
//void Euler::SetBCType(BCtype xbc, BCtype ybc) {
//    xBC = xbc;
//    yBC = ybc;
//}

/**
 * @brief Initialise conservative variable vectors
 *
 *	      Initialise u, uBc, uExact, xCells, and yCells vectors to satisfy initial conditions
 */
void Euler::Initialise() {

    u1 = new VectA4;
    u1->SetSize(Nx, Ny);
    uBc1 = new VectA4;
    uBc1->SetSize(Nx + 2 * NxGhost, Ny + 2 * NyGhost);
    u2 = new VectA4;
    u2->SetSize(Nx, Ny);
    uBc2 = new VectA4;
    uBc2->SetSize(Nx + 2 * NxGhost, Ny + 2 * NyGhost);
    wExact = new VectA4;
    wExact->SetSize(Nx, Ny);

    xCells = new VDouble(N); // x coordinates
    yCells = new VDouble(N); // y coordinates
    phi = new VDouble(N);
    phiBc = new VDouble((Nx + 2 * NxGhost) * (Ny + 2 * NyGhost));

}

/**
 * @brief Construction of computational domain
 */
void Euler::ConstructFVDomain(VDouble& xCellCentres, VDouble& yCellCentres) const {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            xCellCentres[i + j*Nx] = xLeftDomain + dx * (i + 0.5);
            yCellCentres[i + j*Nx] = yBottomDomain + dy * (j + 0.5);
        }
    }
}

/**
 * @brief Set initial conditions
 *
 *        Toro's tests 1, 2, 3, 4, 5, 6, 7 (Page 334), and other tests
 */
void Euler::SetInitialConditions1(const VDouble& xCellCentres, const VDouble& yCellCentres, VectA4& uIC, double& finalt, double& xDiscontinuity, double& yDiscontinuity,
                                  double& gammaL, double& gammaR, double& pinfL, double& pinfR) {

    A4Double wll, wl, wr, Qll, Ql, Qr;

    switch(testNumber) {
        case(1):
        {
            finalt = 0.25;
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.4;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 1.0};
            wr = {0.125, 0.0, 0.0, 0.1};

            break;
        }

        case(2):
        {
            finalt = 0.15;
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.4;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, -2.0, 0.0, 0.4};
            wr = {1.0, 2.0, 0.0, 0.4};

            break;
        }

        case(3):
        {
            finalt = 0.012;
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.4;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 1000.0};
            wr = {1.0, 0.0, 0.0, 0.01};

            break;
        }
        /// Fedkiw Test B
        case(4):
        {
            finalt = 0.0012;
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.67;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wll = {1.3333, 0.3535 * sqrt(1e5), 0.0, 1.5e5};
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            RP_Shock = wll;

            break;
        }
        /// Wang Test B
        case (5):
        {
            finalt = 0.000099; // 0.001; // 0.0014; 1.01e-4
            // 0.00078268328 + 0.00160407972 + 0.00044724759
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.67;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            // wll = {1.3333, 0.3535 * sqrt(1e5), 0.0, 1.5e5};
            wll = {5.92593, 6220.51, 0.0, 4.665e7}; // for Mach 10 shock
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }
        case (6):
        {
            finalt = 0.0014;
            yDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.67;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wll = {1.3333, 0.0, 0.3535 * sqrt(1e5), 1.5e5};
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }
        case (7):
        {
            finalt = 0.0014;
            xDiscontinuity = 0.5;
            gammaL = 1.4;
            gammaR = 1.67;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Neumann;

            wll = {1.3333, 0.3535 * sqrt(1e5), 0.0, 1.5e5};
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }
        /// Helium Bubble Collapse
        case(8):
        {
            finalt = 427e-6;
            xDiscontinuity = 225;
            gammaL = 1.4;
            gammaR = 1.67;
            pinfL = 0.0;
            pinfR = 0.0;

            xBC = Neumann;
            yBC = Reflective;

            wll = {1.0, 0.0, 0.0, 1.0};
            wl = {1.3764, -0.394, 0.0, 1.5698};
            wr = {0.138, 0.0, 0.0, 1.0};

            break;
        }

        default:
        {
            std::cerr << "Test case number invalid.";
            std::exit(EXIT_FAILURE);
        }
    }

    // Save Left and Right states of the RP
    RP_LeftState  = wl;

//    double C = (wl[3])/pow(wl[0], gammaL);
//    wr[0] = exp(log((wr[3])/C)/gammaL);

    convertVarsL.SetVariables((unsigned int) nVar, (double) gammaL, (double) pinfL);
    convertVarsR.SetVariables((unsigned int) nVar, (double) gammaR, (double) pinfR);
    ghostFluidMethod.SetVariables(nVar, Nx, Ny, NxGhost, NyGhost, gammaL, gammaR, pinfL, pinfR, dx, dy, xLeftDomain, yBottomDomain, xBC, yBC);

    Ql = convertVarsL.primitiveToconservative((A4Double) wl);
    Qr = convertVarsL.primitiveToconservative((A4Double) wr);

    /// Assign initial condition values
    if (testNumber == 4) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx] > 0.05 && xCellCentres[i + j*Nx] <= xDiscontinuity) { uIC(i, j) = Ql; }
                else { uIC(i, j) = Qr; }

            }
        }

    } else if (testNumber == 5) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx] > 0.4 && xCellCentres[i + j*Nx] <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

            }
        }

    } else if (testNumber == 6) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( yCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (yCellCentres[i + j*Nx] > 0.4 && yCellCentres[i + j*Nx] <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

            }
        }

    } else if (testNumber == 7) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx]/cos(35*M_PI/180) <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx]/cos(35*M_PI/180) > 0.4 && xCellCentres[i + j*Nx]/cos(35*M_PI/180) <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

            }
        }

    } else if (testNumber == 8) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= xDiscontinuity ) { uIC(i, j) = Qll; }
                else if (-25 + sqrt((xCellCentres[i + j*Nx] - 175)*(xCellCentres[i + j*Nx] - 175) + yCellCentres[i + j*Nx]*yCellCentres[i + j*Nx]) <= 0) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

            }
        }

    } else {

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                uIC(i, j) = Ql;

            }
        }
    }
}

void Euler::SetInitialConditions2(const VDouble& xCellCentres, const VDouble& yCellCentres, VectA4& uIC, VDouble& phiIC) {

    A4Double wll, wl, wr, Qll, Ql, Qr;

    switch(testNumber) {
        case(1):
        {
            wl = {1.0, 0.0, 0.0, 1.0};
            wr = {0.125, 0.0, 0.0, 0.1};

            break;
        }

        case(2):
        {
            wl = {1.0, -2.0, 0.0, 0.4};
            wr = {1.0, 2.0, 0.0, 0.4};

            break;
        }

        case(3):
        {
            wl = {1.0, 0.0, 0.0, 1000.0};
            wr = {1.0, 0.0, 0.0, 0.01};

            break;
        }
        /// Fedkiw Test B & Wang Test B
        case 4 ... 5:
        {
            // wll = {1.3333, 0.3535 * sqrt(1e5), 0.0, 1.5e5};
            wll = {5.92593, 6220.51, 0.0, 4.665e7}; // for Mach 10 shock
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }
        case(6):
        {
            wll = {1.3333, 0.0, 0.3535 * sqrt(1e5), 1.5e5};
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }
        case(7):
        {
            wll = {1.3333, 0.3535 * sqrt(1e5), 0.0, 1.5e5};
            wl = {1.0, 0.0, 0.0, 1.0e5};
            wr = {0.1379, 0.0, 0.0, 1.0e5};

            break;
        }

        case(8):
        {
            wll = {1.0, 0.0, 0.0, 1.0};
            wl = {1.3764, -0.394, 0.0, 1.5698};
            wr = {0.138, 0.0, 0.0, 1.0};

            break;
        }

        default:
        {
            std::cerr << "Test case number invalid.";
            std::exit(EXIT_FAILURE);
        }
    }

    RP_RightState = wr;
//    double C = (wr[3])/pow(wr[0], GammaR);
//    wl[0] = exp(log((wl[3])/C)/GammaR);

    /// Assign initial condition values
    Ql = convertVarsR.primitiveToconservative((A4Double) wl);
    Qr = convertVarsR.primitiveToconservative((A4Double) wr);

    if (testNumber == 4) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx] > 0.05 && xCellCentres[i + j*Nx] <= RP_xDisc) { uIC(i, j) = Ql; }
                else { uIC(i, j) = Qr; }

                phiIC[i + j*Nx] = xCellCentres[i + j*Nx] - RP_xDisc;

            }
        }

    } else if (testNumber == 5) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx] > 0.4 && xCellCentres[i + j*Nx] <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

                phiIC[i + j*Nx] = 0.1 - std::abs(xCellCentres[i + j*Nx] - RP_xDisc);

            }
        }

    } else if (testNumber == 6) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( yCellCentres[i + j*Nx] <= 0.05 ) { uIC(i, j) = Qll; }
                else if (yCellCentres[i + j*Nx] > 0.4 && yCellCentres[i + j*Nx] <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

                phiIC[i + j*Nx] = 0.1 - std::abs(yCellCentres[i + j*Nx] - RP_yDisc);

            }
        }

    } else if (testNumber == 7) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx]/cos(35*M_PI/180) <= 0.05 ) { uIC(i, j) = Qll; }
                else if (xCellCentres[i + j*Nx]/cos(35*M_PI/180) > 0.4 && xCellCentres[i + j*Nx]/cos(35*M_PI/180) <= 0.6) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

                phiIC[i + j*Nx] = (0.1 - std::abs(xCellCentres[i + j*Nx] - RP_xDisc))/cos(35*M_PI/180);

            }
        }

    } else if (testNumber == 8) {

        Qll = convertVarsL.primitiveToconservative((A4Double) wll);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] <= RP_xDisc ) { uIC(i, j) = Qll; }
                else if (-25 + sqrt((xCellCentres[i + j*Nx] - 175)*(xCellCentres[i + j*Nx] - 175) + yCellCentres[i + j*Nx]*yCellCentres[i + j*Nx]) <= 0) { uIC(i, j) = Qr; }
                else { uIC(i, j) = Ql; }

                phiIC[i + j*Nx] = 25 - sqrt((xCellCentres[i + j*Nx] - 175)*(xCellCentres[i + j*Nx] - 175) + yCellCentres[i + j*Nx]*yCellCentres[i + j*Nx]);

            }
        }

    } else {

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                uIC(i, j) = Qr;

                phiIC[i + j*Nx] = xCellCentres[i + j*Nx] - RP_xDisc;

            }
        }
    }
}

/**
 * @brief Set boundary conditions
 */
void Euler::SetBoundaryConditions(const VectA4& uIC, VectA4& ubc) const {

    /// Copy the valid state
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            ubc(i + NxGhost, j + NyGhost) = uIC(i, j);
        }
    }

    switch(xBC)
    {
        case(Neumann):
        {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < NxGhost; i++) {
                    ubc(i, j + NyGhost) = uIC(0, j);
                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(Nx - 1, j);
                }
            }
            break;
        }
        case(Periodic):
        {
            /// Periodic BC in x-direction
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < NxGhost; i++) {
                    ubc(i, j + NyGhost) = uIC(Nx - NxGhost + i, j);
                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(NxGhost + i, j);
                }
            }
            break;
        }
        default:
        {
            std::cerr << "Boundary conditions not recognized.";
            std::exit(EXIT_FAILURE);
        }
    }

    switch(yBC)
    {
        case(Neumann):
        {
            /// Neumann boundary conditions
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {
                    ubc(i + NxGhost, j) = uIC(i, 0);
                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, Ny - 1);
                }
            }
            break;
        }
        case(Periodic):
        {
            /// Implement here other boundary conditions ... (reflective for example)
            /// Periodic BC in y-direction
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {
                    ubc(i + NxGhost, j) = uIC(i, Ny - NyGhost + j);
                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, NyGhost + j);
                }
            }
            break;
        }
        case(Reflective):
        {
            /// Reflective BC in y-direction
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {

                    for (int k = 0; k < nVar; k++) {
                        /// Reflecting u_y
                        if (k == 2) {
                            ubc(i + NxGhost, j, k) = -uIC(i, NyGhost - 1 - j, k);
                            ubc(i + NxGhost, Ny + NyGhost + j, k) = -uIC(i, Ny - 1 - j, k);
                        } else {
                            ubc(i + NxGhost, j, k) = uIC(i, NyGhost - 1 - j, k);
                            ubc(i + NxGhost, Ny + NyGhost + j, k) = uIC(i, Ny - 1 - j, k);
                        }
                    }

                }
            }
            break;
        }
        default:
        {
            std::cerr << "Boundary conditions not recognized.";
            std::exit(EXIT_FAILURE);
        }
    }
}

/**
 * @brief Solution to the user-defined problem to compute conservative variable values at time t + dt
 *        Update all the steps together using given EoS to solve for the conservative variables iteratively
 * 	      Create a new RiemannSolver instance and solve the Riemann problem
 */
void Euler::Update() {

    /// Construct computational domain
    cout << "Constructing spatial domain... " << endl;
    ConstructFVDomain((*xCells), (*yCells));

    /// Set initial conditions
    cout << "Setting initial conditions... " << endl;
    SetInitialConditions1((*xCells), (*yCells), (*u1), T, RP_xDisc, RP_yDisc, GammaL, GammaR, PinfL, PinfR);
    SetInitialConditions2((*xCells), (*yCells), (*u2), (*phi));  // can I get away without setting phiBc here (*phiBc)

    /// 1D Implementation to get exact solution
    if (testNumber == 4) exactRPSolver.SetVariables(RP_Shock, RP_RightState, GammaL, GammaR, PinfL, PinfR);
    else exactRPSolver.SetVariables(RP_LeftState, RP_RightState, GammaL, GammaR, PinfL, PinfR);

    /// Initialise time loop
    int iter = 0;
    double time = 0;

    /// Initialise timing output dat files
    std::remove("Data/Timing.dat"); // delete file
    std::ofstream outTime("Data/Timing.dat");

    /// FluxUpdate routines
    fluxUpdate1.SetVariables((unsigned int) Nx, (unsigned int) Ny, (double) Cfl, (unsigned int) NxGhost, (unsigned int) NyGhost,
                             (unsigned int) nVar, (double) T, GammaL, PinfL);
    fluxUpdate1.Initialise();

    fluxUpdate2.SetVariables((unsigned int) Nx, (unsigned int) Ny, (double) Cfl, (unsigned int) NxGhost, (unsigned int) NyGhost,
                             (unsigned int) nVar, (double) T, GammaR, PinfR);
    fluxUpdate2.Initialise();

    /// Timing parameters
    clock_t start, end;
    double elapsed_time = 0.0;
    PairDouble timerxUpdate, timeryUpdate; // dimensional split timings

    cout << "Starting time loop... " << endl;
    do {
        /// Interface indices
        // cout << "Finding interface location... " << endl;
        interfaceList = ghostFluidMethod.FindInterfaceLocation((*phi));

        /// Apply Boundary Conditions
        // cout << "Setting BCs... " << endl;
        SetBoundaryConditions((*u1), (*uBc1));
        SetBoundaryConditions((*u2), (*uBc2));

        // cout << "Setting ghost state... " << endl;
        ghostFluidMethod.SetGhostFluid((*u1), (*u2), (*uBc1), (*uBc2), (*phi), (*phiBc), (*xCells), (*yCells), interfaceList);

        /// Compute time step
        // cout << "Computing time-step... " << endl;
        start = clock();
        fluxUpdate1.ComputeDt((*uBc1), dx, dy, time, dt1);
        fluxUpdate2.ComputeDt((*uBc2), dx, dy, time, dt2);
        if (dt1 < 1e-10) dt = dt2;
        else if (dt2 < 1e-10) dt = dt1;
        else dt = std::min<double>(dt1, dt2);
        end = clock();
        elapsed_time = (end - start)/(double)CLOCKS_PER_SEC ;

        // cout << "Updating level set... " << endl;
        ghostFluidMethod.UpdateLevelSet((*uBc1), (*uBc2), dt, (*phi), (*phiBc));
        ghostFluidMethod.Reinitialisation(interfaceList, (*phi), (*phiBc)); // every 5 - 10 time steps

        /// Dimensional split x-direction
        // Compute numerical fluxes and update solution
        // cout << "Update solution with numerical fluxes... " << endl;
        fluxUpdate1.UpdatewithFluxes((*uBc1), (*u1), dx, dt, 1, timerxUpdate);
        fluxUpdate2.UpdatewithFluxes((*uBc2), (*u2), dx, dt, 1, timerxUpdate);

        /// Dimensional split y-direction
        // cout << "Apply BCs for y split... " << endl;
        SetBoundaryConditions((*u1), (*uBc1));
        SetBoundaryConditions((*u2), (*uBc2));

        // Compute numerical fluxes and update solution
        // cout << "Update solution with numerical fluxes for y split... " << endl;
        fluxUpdate1.UpdatewithFluxes((*uBc1), (*u1), dy, dt, 2, timeryUpdate);
        fluxUpdate2.UpdatewithFluxes((*uBc2), (*u2), dy, dt, 2, timeryUpdate);

        iter += 1;
        time += dt;
        std::cout << "Iteration # " << iter << " Time = " << time << std::endl;

        // Output the data
        name = "Mach10T" + std::to_string(time) + ".dat";
        std::ofstream output("Data/Mach10T" + std::to_string(time) + ".dat");
        GenerateData((*xCells), (*yCells), (*u1), (*u2), (*phi), (*wExact));

        outTime << iter << " " << elapsed_time << " " << timerxUpdate.first << " " << timerxUpdate.second << " " << timeryUpdate.first << " " << timeryUpdate.second << std::endl;

    } while (time < T);

    if (testNumber == 4) {
        double Sspeed = (RP_Shock[0]*RP_Shock[1] - RP_LeftState[0]*RP_LeftState[1])/(RP_Shock[0] - RP_LeftState[0]);
        double deductTime = 0.45/Sspeed;
        time -= deductTime;
    }

    /// Compute Exact Solution
    if (testNumber == 1 || testNumber == 2 || testNumber == 3 || testNumber == 4) {
        exactRPSolver.ComputeExactSolution(time, (*xCells), RP_xDisc, (*wExact));
    }

    /// Write to output file
    cout << "Writing data file... " << endl;
    GenerateData((*xCells), (*yCells), (*u1), (*u2), (*phi), (*wExact));
}

/**
 *  @brief  Output x & y coordinates and conservative variable values to dat file
 */
void Euler::GenerateData(const VDouble& xCellCentres, const VDouble& yCellCentres, const VectA4& uOut1, const VectA4& uOut2, const VDouble& phiOut, const VectA4& wExactOut) {

    std::ofstream outFile;

    outFile.open(name);

    outFile << "# This file is the output of Euler.cpp" << std::endl;

    outFile << "# " << std::endl;

    VectA4 prim1, prim2;
    prim1.SetSize(Nx, Ny); // primitive for output
    prim2.SetSize(Nx, Ny);
    double rhoExact, uExact, pExact, eExact;

    equationofStatesL.SetVariables(GammaL, PinfL);
    equationofStatesR.SetVariables(GammaR, PinfR);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            prim1(i, j) = convertVarsL.conservativeToprimitive(uOut1(i, j));
            prim2(i, j) = convertVarsR.conservativeToprimitive(uOut2(i, j));

            if (testNumber == 1 || testNumber == 2 || testNumber == 3 || testNumber == 4) {
                rhoExact = wExactOut(i, j, 0);
                uExact = wExactOut(i, j, 1);
                pExact = wExactOut(i, j, 3);
                if (i < interfaceList[0].first) eExact = equationofStatesL.computeInternalEnergyFromEoS(rhoExact, pExact);
                else eExact = equationofStatesR.computeInternalEnergyFromEoS(rhoExact, pExact);
            } else {
                rhoExact = 0.0;
                uExact = 0.0;
                pExact = 0.0;
                eExact = 0.0;
            }

            outFile << xCellCentres[i + j*Nx] << " " << yCellCentres[i + j*Nx]
            << " " << prim1(i, j, 0) << " " << prim1(i, j, 1) << " " << prim1(i, j, 2) << " " << prim1(i, j, 3) << " " << equationofStatesL.computeInternalEnergyFromEoS(prim1(i, j, 0), prim1(i, j, 3))
            << " " << prim2(i, j, 0) << " " << prim2(i, j, 1) << " " << prim2(i, j, 2) << " " << prim2(i, j, 3) << " " << equationofStatesR.computeInternalEnergyFromEoS(prim2(i, j, 0), prim2(i, j, 3))
            << " " << phiOut[i + j*Nx] << " " << rhoExact << " " << uExact << " " << pExact << " " << eExact << "\n";
        }
        // outFile << "\n";
    }

    outFile.close();
}