/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : GhostFluidMethod.cpp                       */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 23/01/2022                                                */
/*                                                                    */
/*   Programmer :                                        */
/*                                                                    */
/*   Description : Ghost fluid method routines                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "GhostFluidMethod.h"
#include <cmath>

using namespace std;

/**
 * @brief Constructor of the GhostFluidMethod class
 */
GhostFluidMethod::GhostFluidMethod() {
}

/**
 * @brief Destructor of the GhostFluidMethod class
 */
GhostFluidMethod::~GhostFluidMethod() {
}

/**
 * @brief Set variables
 *
 */
void GhostFluidMethod::SetVariables(unsigned int nvar, unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost,
                                    double gammaL, double gammaR, double pinfL, double pinfR, double deltax, double deltay, double xLeft, double yBot, BCtype xbc, BCtype ybc) {
    Nx = nx;
    Ny = ny;
    NxGhost = nxGhost;
    NyGhost = nyGhost;
    GammaL = gammaL;
    GammaR = gammaR;
    PinfL = pinfL;
    PinfR = pinfR;
    nVar = nvar;
    dx = deltax;
    dy = deltay;
    xLeftDomain = xLeft;
    yBottomDomain = yBot;
    xBC = xbc;
    yBC = ybc;

    convertVars1.SetVariables(nVar, GammaL, PinfL);
    convertVars2.SetVariables(nVar, GammaR, PinfR);
}

/// Find interface cells and return vector of pairs i, j index (not including ghost cells)
VectPair GhostFluidMethod::FindInterfaceLocation(const VDouble& phi) const {

    VectPair interfaceInt;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (phi[i + j*Nx] * phi[i+1 + j*Nx] <= 0.0)  {
                interfaceInt.emplace_back(i, j);
                interfaceInt.emplace_back(i+1, j);
            }
            if (phi[i + j*Nx] * phi[i + (j+1)*Nx] <= 0.0) {
                interfaceInt.emplace_back(i, j);
                interfaceInt.emplace_back(i, j+1);
            }
        }
    }
    // cout << "interfaceInt length " << interfaceInt.size() << endl;
    // cout << interfaceInt[0].first << " " << interfaceInt[0].second << " " << interfaceInt[1].first << " " << interfaceInt[1].second << endl;

    return interfaceInt;
}

/// Set ghost fluid state
void GhostFluidMethod::SetGhostFluid(VectA4& u1, VectA4& u2, VectA4& ubc1, VectA4& ubc2, const VDouble& phi, VDouble& phibc,
                                     const VDouble& xCellCentres, const VDouble& yCellCentres, const VectPair& interfaceInt) {
    /// Copy the valid state
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            phibc[i + NxGhost + (j + NyGhost)*(Nx + 2 * NxGhost)] = phi[i + j*Nx];
        }
    }

    /// Zero Neumann BC for level set in x and y-directions
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < NxGhost; i++) {
            phibc[i + (j + NyGhost)*(Nx + 2 * NxGhost)] = phi[j*Nx];
            phibc[Nx + NxGhost + i + (j + NyGhost)*(Nx + 2 * NxGhost)] = phi[Nx - 1 + j*Nx];
        }
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < NyGhost; j++) {
            phibc[i + NxGhost + j*(Nx + 2 * NxGhost)] = phi[i];
            phibc[i + NxGhost + (Ny + NyGhost + j)*(Nx + 2 * NxGhost)] = phi[i + (Ny - 1)*Nx];
        }
    }

    VectA4 w1, w2, wbc1, wbc2;
    std::array<double, 2> normal{0, 0};
    w1.SetSize(Nx, Ny);
    w2.SetSize(Nx, Ny);
    wbc1.SetSize(Nx + 2 * NxGhost, Ny + 2 * NyGhost);
    wbc2.SetSize(Nx + 2 * NxGhost, Ny + 2 * NyGhost);

    // double ds = 1.0 * std::min<double>(dx, dy);
    // std::vector<std::array<double, 2> > newVel((Nx + 2 * NxGhost) * (Ny + 2 * NyGhost));

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            w1(i, j) = convertVars1.conservativeToprimitive(u1(i, j));
            w2(i, j) = convertVars2.conservativeToprimitive(u2(i, j));
        }
    }

    for (int i = 0; i < Nx + 2 * NxGhost; i++) {
        for (int j = 0; j < Ny + 2 * NyGhost; j++) {
            wbc1(i, j) = convertVars1.conservativeToprimitive(ubc1(i, j));
            wbc2(i, j) = convertVars2.conservativeToprimitive(ubc2(i, j));
        }
    }

    // IterativeMethod(2, phibc, w1, w2, wbc1, wbc2);  // in case sampling ghost cells in the mixed RP

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            if (isInterfaceCell) {

                double phival = phibc[i + j*(Nx + 2 * NxGhost)];
                normal = GetNormalVector(phibc, i, j); // Centred scheme

                A4Double interfaceW1, interfaceW2;
                GetInterpolatedMixedRPStates(xCellCentres, yCellCentres, phival, normal, wbc1, wbc2, i, j, interfaceW1, interfaceW2);

                std::array<double, 2> vel1{interfaceW1[1], interfaceW1[2]}, vel2{interfaceW2[1], interfaceW2[2]};
                double vn1 = normal[0]*vel1[0] + normal[1]*vel1[1];
                double vn2 = normal[0]*vel2[0] + normal[1]*vel2[1];
                std::array<double, 2> vel1tan{vel1[0] - vn1*normal[0], vel1[1] - vn1*normal[1]}, vel2tan{vel2[0] - vn2*normal[0], vel2[1] - vn2*normal[1]};

                A4Double wL, wR, wlstar, wrstar;

                wL[0] = interfaceW1[0];
                wL[1] = vn1;
                wL[2] = 0.0;
                wL[3] = std::max(1e-6, interfaceW1[3]);

                wR[0] = interfaceW2[0];
                wR[1] = vn2;
                wR[2] = 0.0;
                wR[3] = std::max(1e-6, interfaceW2[3]);

                exactRPSolver.SetVariables(wL, wR, GammaL, GammaR, PinfL, PinfR);
                ComputeIntermediateStates(wlstar, wrstar);

                std::array<double, 2> u1star{wlstar[1] * normal[0] + vel1tan[0], wlstar[1] * normal[1] + vel1tan[1]}, u2star{wrstar[1] * normal[0] + vel2tan[0], wrstar[1] * normal[1] + vel2tan[1]};

                if (phival > 0.0) {  // Set left star state to Material 1 ghost cell, right of interface
                    wbc1(i, j, 0) = wlstar[0];
                    wbc1(i, j, 1) = u1star[0];
                    wbc1(i, j, 2) = u1star[1];
                    wbc1(i, j, 3) = wlstar[3];
                } else {  // Set right star state to Material 2 ghost cell, left of interface
                    wbc2(i, j, 0) = wrstar[0];
                    wbc2(i, j, 1) = u2star[0];
                    wbc2(i, j, 2) = u2star[1];
                    wbc2(i, j, 3) = wrstar[3];
                }

                // newVel[i + j*(Nx + 2 * NxGhost)] = {wlstar[1] * normal[0], wlstar[1] * normal[1]};
            }
        }
    }

    // IterativeMethod(6, phibc, w1, w2, wbc1, wbc2);
    PopulatingGhostCells(interfaceInt, phibc, wbc1, wbc2);

    for (int i = 0; i < Nx + 2 * NxGhost; i++) {
        for (int j = 0; j < Ny + 2 * NyGhost; j++) {
            ubc1(i, j) = convertVars1.primitiveToconservative(wbc1(i, j));
            ubc2(i, j) = convertVars2.primitiveToconservative(wbc2(i, j));
        }
    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {
            w1(i-NxGhost, j-NyGhost) = wbc1(i, j);
            w2(i-NxGhost, j-NyGhost) = wbc2(i, j);
            u1(i-NxGhost, j-NyGhost) = ubc1(i, j);
            u2(i-NxGhost, j-NyGhost) = ubc2(i, j);
        }
    }

//    std::ofstream output;
//
//    output.open("Extrapolate.dat");
//
//    output << "# This file is the output of Euler.cpp" << std::endl;
//
//    output << "# " << std::endl;
//
//    for (int i = 0; i < Nx; i++) {
//        for (int j = 0; j < Ny; j++) {
//
//            output << xCellCentres[i + j*Nx] << " " << yCellCentres[i + j*Nx] << " " << w1(i, j, 0) << " " << w1(i, j, 1) << " " << w1(i, j, 2)
//                    << " " << w1(i, j, 3) << " " << w2(i, j, 0) << " " << w2(i, j, 1) << " " << w2(i, j, 2) << " " << w2(i, j, 3)
//                    << " " << phi[i + j*Nx] << "\n";
//        }
//    }

//    output.close();

//    for (int i = 0; i < Nx + 2 * NxGhost; i++) {
//        for (int j = 0; j < Ny + 2 * NyGhost; j++) {
//            if (phibc[i + j*(Nx + 2 * NxGhost)] <= 0.0) newVel[i + j*(Nx + 2 * NxGhost)] = {wbc1(i, j, 1), wbc1(i, j, 2)};
//            else newVel[i + j*(Nx + 2 * NxGhost)] = {wbc2(i, j, 1), wbc2(i, j, 2)};
//        }
//    }
}

/// Bilinear interpolation to set up one-dimensional Riemann problem
void GhostFluidMethod::GetInterpolatedMixedRPStates(const VDouble& xCellCentres, const VDouble& yCellCentres, double phival, const std::array<double, 2>& normal,
                                                    const VectA4& wbc1, const VectA4& wbc2, int i, int j, A4Double& LeftState, A4Double& RightState) const {

    std::array<double, 2> interfacePos{0, 0}, sampleL{0, 0}, sampleR{0, 0};
    std::vector<std::array<double, 4> > wBc1, wBc2;  // new vectors to fit template
    double ds = 1.5 * std::max(dx, dy); // 1.0

    wBc1.resize((Nx + 2 * NxGhost) * (Ny + 2 * NyGhost));
    wBc2.resize((Nx + 2 * NxGhost) * (Ny + 2 * NyGhost));

    for (int k = 0; k < Nx + 2 * NxGhost; k++) {
        for (int l = 0; l < Ny + 2 * NyGhost; l++) {
            wBc1[k + l*(Nx + 2 * NxGhost)] = wbc1(k, l);
            wBc2[k + l*(Nx + 2 * NxGhost)] = wbc2(k, l);
        }
    }

    if (fabs(phival) > ds) phival = std::copysign(ds, phival);

    interfacePos[0] = xCellCentres[i-NxGhost + (j-NyGhost)*Nx] - phival * normal[0];
    interfacePos[1] = yCellCentres[i-NxGhost + (j-NyGhost)*Nx] - phival * normal[1];

    sampleL[0] = interfacePos[0] + ds * normal[0];
    sampleL[1] = interfacePos[1] + ds * normal[1];

    sampleR[0] = interfacePos[0] - ds * normal[0];
    sampleR[1] = interfacePos[1] - ds * normal[1];

    LeftState = grid_bilinear_interpolate<A4Double>(dx, dy, xCellCentres, yCellCentres, xLeftDomain, yBottomDomain, Nx, nVar, wBc1, sampleR);
    RightState = grid_bilinear_interpolate<A4Double>(dx, dy, xCellCentres, yCellCentres, xLeftDomain, yBottomDomain, Nx, nVar, wBc2, sampleL);
}

/// Exact RP Solver to get intermediate states
void GhostFluidMethod::ComputeIntermediateStates(A4Double& wlinter, A4Double& wrinter) {

    double rholstar, alstar, Shl, Stl, Sl;
    double rhorstar, arstar, Shr, Str, Sr;
    double pstar, ustar;

    pstar = exactRPSolver.ComputePStar();
    ustar = exactRPSolver.ComputeUStar(pstar);

    // Wave pattern at interface
    int Pattern;
    exactRPSolver.FindPattern(pstar, Pattern);
    exactRPSolver.ComputeRhoStar(pstar, ustar, Pattern, rholstar, alstar, Shl, Stl, Sl, rhorstar, arstar, Shr, Str, Sr);

    wlinter[0] = rholstar;
    wlinter[1] = ustar;
    wlinter[2] = 0;
    wlinter[3] = pstar;

    wrinter[0] = rhorstar;
    wrinter[1] = ustar;
    wrinter[2] = 0;
    wrinter[3] = pstar;
}

/// Apply boundary conditions and update level set first order in space and time
void GhostFluidMethod::UpdateLevelSet(const VectA4& u1, const VectA4& u2, const double& dt, VDouble& phi, VDouble& phibc) {

    double phipx, phimx, phipy, phimy, tempPhix, tempPhiy;
    A4Double w1, w2;

    /// Euler upwind first order update
    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {
            phipx = (phibc[i+1 + j*(Nx + 2 * NxGhost)] - phibc[i + j*(Nx + 2 * NxGhost)])/dx;
            phimx = (phibc[i + j*(Nx + 2 * NxGhost)] - phibc[i-1 + j*(Nx + 2 * NxGhost)])/dx;
            phipy = (phibc[i + (j+1)*(Nx + 2 * NxGhost)] - phibc[i + j*(Nx + 2 * NxGhost)])/dy;
            phimy = (phibc[i + j*(Nx + 2 * NxGhost)] - phibc[i + (j-1)*(Nx + 2 * NxGhost)])/dy;

            if (phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) {  // updating with star state instead, change to < 0 for Mach test and 2D, with reinitialisation
                w1 = convertVars1.conservativeToprimitive(u1(i, j));
                if (w1[1] <= 0.0) tempPhix = phipx;
                else tempPhix = phimx;
                if (w1[2] <= 0.0) tempPhiy = phipy;
                else tempPhiy = phimy;
                phi[i-NxGhost + (j-NyGhost)*Nx] = phibc[i + j*(Nx + 2 * NxGhost)] - dt * (w1[1] * tempPhix + w1[2] * tempPhiy);
            } else {
                w2 = convertVars2.conservativeToprimitive(u2(i, j));
                if (w2[1] <= 0.0) tempPhix = phipx;
                else tempPhix = phimx;
                if (w2[2] <= 0.0) tempPhiy = phipy;
                else tempPhiy = phimy;
                phi[i-NxGhost + (j-NyGhost)*Nx] = phibc[i + j*(Nx + 2 * NxGhost)] - dt * (w2[1] * tempPhix + w2[2] * tempPhiy);
            }
        }
    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {
            phibc[i + j*(Nx + 2 * NxGhost)] = phi[i-NxGhost + (j-NyGhost)*Nx];
        }
    }
}

/// Populating ghost cells using fast sweeping
void GhostFluidMethod::PopulatingGhostCells(const VectPair& interfaceInt, const VDouble& phibc, VectA4& Q1, VectA4& Q2) const {

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            if (! isInterfaceCell) {
                if (phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) {  // Material 2 is ghost here
                    for (int n = 0; n < nVar; n++) {
                        Q2(i, j, n) = std::copysign(1e100, Q2(i, j, n));
                    }
                } else {  // Material 1 is ghost here
                    for (int n = 0; n < nVar; n++) {
                        Q1(i, j, n) = std::copysign(1e100, Q1(i, j, n));
                    }
                }
            }

        }
    }

    std::array<double, 2> normal{0, 0};

    for (int i = NxGhost; i < Nx + NxGhost; i++) {  // x sweep positive direction

        for (int j = NyGhost; j < Ny + NyGhost; j++) {

            bool isInterfaceCell = false;
            bool negativegradient = false;

            if (phibc[i + j*(Nx + 2 * NxGhost)] - phibc[i-1 + j*(Nx + 2 * NxGhost)] < 0.0) negativegradient = true;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            normal = GetNormalVector(phibc, i ,j);

            if (! negativegradient) {
                if (! isInterfaceCell) {
                    if (phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) Q1(i, j) = FastSweepingQPos(Q1, phibc, i, j, normal);  // Q1 is ghost here
                    else Q2(i, j) = FastSweepingQNeg(Q2, phibc, i, j, normal);  // Q2 is ghost here
                }
            }
        }
    }

    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {

        for (int j = NyGhost; j < Ny+NyGhost; j++) {  // y sweep in positive direction

            bool isInterfaceCell = false;
            bool negativegradient = false;

            if (phibc[i + (j)*(Nx + 2 * NxGhost)] - phibc[i + (j-1)*(Nx + 2 * NxGhost)] < 0.0) negativegradient = true;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            normal = GetNormalVector(phibc, i ,j);

            if (! negativegradient) {
                if (! isInterfaceCell) {
                    if (phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) Q1(i, j) = FastSweepingQPos(Q1, phibc, i, j, normal);  // Q1 is ghost here
                    else Q2(i, j) = FastSweepingQNeg(Q2, phibc, i, j, normal);  // Q2 is ghost here
                }
            }
        }
    }

    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {  // x sweep in negative direction

        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {

            bool isInterfaceCell = false;
            bool negativegradient = false;

            if (phibc[(i) + j*(Nx + 2 * NxGhost)] - phibc[i+1 + j*(Nx + 2 * NxGhost)] < 0.0) negativegradient = true;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            normal = GetNormalVector(phibc, i , j);

            if (! negativegradient) {
                if (! isInterfaceCell) {
                    if (phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) Q1(i, j) = FastSweepingQPos(Q1, phibc, i, j, normal);  // Q1 is ghost here
                    else Q2(i, j) = FastSweepingQNeg(Q2, phibc, i, j, normal);  // Q2 is ghost here
                }
            }
        }
    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {

        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {  // y sweep negative direction

            bool isInterfaceCell = false;
            bool negativegradient = false;

            if (phibc[i + (j)*(Nx + 2 * NxGhost)] - phibc[i + (j+1)*(Nx + 2 * NxGhost)] < 0.0) negativegradient = true;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            normal = GetNormalVector(phibc, i ,j);

            if (! negativegradient) {
                if (! isInterfaceCell) {
                    if (phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) Q1(i, j) = FastSweepingQPos(Q1, phibc, i, j, normal);  // Q1 is ghost here
                    else Q2(i, j) = FastSweepingQNeg(Q2, phibc, i, j, normal);  // Q2 is ghost here
                }
            }
        }
    }
}

/// Fast sweeping Q routine for phi > 0
A4Double GhostFluidMethod::FastSweepingQPos(const VectA4& Q, const VDouble& phibc, int i, int j, const std::array<double, 2>& normal) const {

    A4Double Qtmp;
    double a, b;

    for (int n = 0; n < nVar; n++) {

        if (phibc[i-1 + j*(Nx + 2 * NxGhost)] < phibc[i+1 + j*(Nx + 2 * NxGhost)]) a = Q(i-1, j, n); // Qx
        else a =  Q(i+1, j, n);
        if (phibc[i + (j-1)*(Nx + 2 * NxGhost)] < phibc[i + (j+1)*(Nx + 2 * NxGhost)]) b = Q(i, j-1, n); // Qy
        else b = Q(i, j+1, n);

        Qtmp[n] = (a * dy*normal[0] + b * dx*normal[1]) / (dy*normal[0] + dx*normal[1]);
    }

    return Qtmp;
}

/// Fast sweeping Q routine for phi < 0
A4Double GhostFluidMethod::FastSweepingQNeg(const VectA4& Q, const VDouble& phibc, int i, int j, const std::array<double, 2>& normal) const {

    A4Double Qtmp;
    double a, b;

    for (int n = 0; n < nVar; n++) {

        if (phibc[i-1 + j*(Nx + 2 * NxGhost)] > phibc[i+1 + j*(Nx + 2 * NxGhost)]) a = Q(i-1, j, n); // Qx
        else a =  Q(i+1, j, n);
        if (phibc[i + (j-1)*(Nx + 2 * NxGhost)] > phibc[i + (j+1)*(Nx + 2 * NxGhost)]) b = Q(i, j-1, n); // Qy
        else b = Q(i, j+1, n);

        Qtmp[n] = (a * dy*normal[0] + b * dx*normal[1]) / (dy*normal[0] + dx*normal[1]);
    }

    return Qtmp;
}

/// Reinitialisation of the level set function using the fast sweeping method
void GhostFluidMethod::Reinitialisation(const VectPair& interfaceInt, VDouble& phi, VDouble& phibc) const {

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            if (! isInterfaceCell) {
                phibc[i + j*(Nx + 2 * NxGhost)] = std::copysign(1e100, phibc[i + j*(Nx + 2 * NxGhost)]);
            }

        }
    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            // if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
            if (! isInterfaceCell) {
                if (phibc[i + j * (Nx + 2 * NxGhost)] > 0.0) phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
                else phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
            }
        }
    }

    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {
        for (int j = NyGhost; j < Ny+NyGhost; j++) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            //if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
            if (! isInterfaceCell) {
                if (phibc[i + j * (Nx + 2 * NxGhost)] > 0.0) phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
                else phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
            }

        }
    }

    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {
        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            // if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
            if (! isInterfaceCell) {
                if (phibc[i + j * (Nx + 2 * NxGhost)] > 0.0) phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
                else phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
            }

        }
    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {

            bool isInterfaceCell = false;

            for (const auto& cell : interfaceInt) {
                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
            }

            // if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] > 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
            if (! isInterfaceCell) {
                if (phibc[i + j * (Nx + 2 * NxGhost)] > 0.0) phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodPositive(phibc, i, j);
                else phibc[i + j * (Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
            }

        }
    }

//    for (int i = NxGhost; i < Nx + NxGhost; i++) {
//        for (int j = NyGhost; j < Ny + NyGhost; j++) {
//
//            bool isInterfaceCell = false;
//
//            for (const auto& cell : interfaceInt) {
//                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
//            }
//
//            if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
//
//        }
//    }
//
//    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {
//        for (int j = NyGhost; j < Ny+NyGhost; j++) {
//
//            bool isInterfaceCell = false;
//
//            for (const auto& cell : interfaceInt) {
//                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
//            }
//
//            if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
//
//        }
//    }
//
//    for (int i = Nx+NxGhost-1; i >= NxGhost; i--) {
//        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {
//
//            bool isInterfaceCell = false;
//
//            for (const auto& cell : interfaceInt) {
//                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
//            }
//
//            if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
//
//        }
//    }
//
//    for (int i = NxGhost; i < Nx + NxGhost; i++) {
//        for (int j = Ny+NyGhost-1; j >= NyGhost; j--) {
//
//            bool isInterfaceCell = false;
//
//            for (const auto& cell : interfaceInt) {
//                if (i-NxGhost == cell.first && j-NyGhost == cell.second) isInterfaceCell = true;
//            }
//
//            if (! isInterfaceCell && phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) phibc[i + j*(Nx + 2 * NxGhost)] = FastSweepingMethodNegative(phibc, i, j);
//
//        }
//    }

    for (int i = NxGhost; i < Nx + NxGhost; i++) {
        for (int j = NyGhost; j < Ny + NyGhost; j++) {
            phi[i-NxGhost + (j-NyGhost)*Nx] = phibc[i + j*(Nx + 2 * NxGhost)];
        }
    }

}

/// Fast sweeping method for phi > 0
double GhostFluidMethod::FastSweepingMethodPositive(VDouble& phibc, int i, int j) const {

    double phival = phibc[i + j*(Nx + 2 * NxGhost)];

    double a = std::min(phibc[i-1 + j*(Nx + 2 * NxGhost)], phibc[i+1 + j*(Nx + 2 * NxGhost)]);
    double b = std::min(phibc[i + (j-1)*(Nx + 2 * NxGhost)], phibc[i + (j+1)*(Nx + 2 * NxGhost)]);
    double phiTmp;

    if (a - b <= -dx) phiTmp = dx + a;
    else if (a - b >= dy) phiTmp = dy + b;
    else {
        double dxdy = dx * dy;
        double dxsq = dx * dx;
        double dysq = dy * dy;

        phiTmp = (a * dysq + b * dxsq + dxdy * sqrt(dxsq + dysq - (a - b) * (a - b))) / (dxsq + dysq);
    }

    return fabs(phival) < fabs(phiTmp) ? phival : phiTmp;
}

/// Fast sweeping method for phi < 0
double GhostFluidMethod::FastSweepingMethodNegative(VDouble& phibc, int i, int j) const {

    double phival = phibc[i + j*(Nx + 2 * NxGhost)];

    double a = std::max(phibc[i-1 + j*(Nx + 2 * NxGhost)], phibc[i+1 + j*(Nx + 2 * NxGhost)]);
    double b = std::max(phibc[i + (j-1)*(Nx + 2 * NxGhost)], phibc[i + (j+1)*(Nx + 2 * NxGhost)]);
    double phiTmp;

    if (b - a <= -dx) phiTmp = a - dx;
    else if (b - a >= dy) phiTmp = b - dy;
    else {
        double dxdy = dx * dy;
        double dxsq = dx * dx;
        double dysq = dy * dy;

        phiTmp = (a * dysq + b * dxsq - dxdy * sqrt(dxsq + dysq - (a - b) * (a - b))) / (dxsq + dysq);
    }

    return fabs(phival) < fabs(phiTmp) ? phival : phiTmp;;
}

/// Velocity Getters
// Get normal vector based on upwind/centred direction of velocity
std::array<double, 2> GhostFluidMethod::GetNormalVector(const VDouble& phibc, int i, int j) const {

    double dphix, dphiy, mag;
    std::array<double, 2> n{0, 0};

    /// Upwind
//    if (w[1] < 0) phix = (phibc[i+1 + j*(Nx + 2 * NxGhost)] - phibc[i + j*(Nx + 2 * NxGhost)])/dx;
//    else phix = (phibc[i + j*(Nx + 2 * NxGhost)] - phibc[i-1 + j*(Nx + 2 * NxGhost)])/dx;

//    if (w[2] < 0) phiy = (phibc[i + (j+1)*(Nx + 2 * NxGhost)] - phibc[i + j*(Nx + 2 * NxGhost)])/dy;
//    else phiy = (phibc[i + j*(Nx + 2 * NxGhost)] - phibc[i + (j-1)*(Nx + 2 * NxGhost)])/dy;

    /// Centred
    dphix = (phibc[i+1 + j*(Nx + 2 * NxGhost)] - phibc[i-1 + j*(Nx + 2 * NxGhost)])/(2*dx);
    dphiy = (phibc[i + (j+1)*(Nx + 2 * NxGhost)] - phibc[i + (j-1)*(Nx + 2 * NxGhost)])/(2*dy);
    mag = sqrt(dphix*dphix + dphiy*dphiy);

    if (mag < 1e-12) {
        n[0] = 1.0;
        n[1] = 0.0;
    } else {
        n[0] = dphix/mag;
        n[1] = dphiy/mag;
    }

    n[0] = dphix;
    n[1] = dphiy;

    return n;
}

//double GhostFluidMethod::GetNormalVelocity(const A4Double& w, const std::array<double, 2>& n) {
//
//    return w[1]*n[0] + w[2]*n[1];
//
//}

//std::array<double, 2> GhostFluidMethod::GetTangentialVelocity(const A4Double& w, const std::array<double, 2>& n, const double& vn) {
//
//    return std::array<double, 2> {w[1] - vn*n[0], w[2] - vn*n[1]};
//
//}

//void GhostFluidMethod::SetBoundaryConditions(const VectA4& uIC, VectA4& ubc) const {
//
//    /// Copy the valid state
//    for (int i = 0; i < Nx; i++) {
//        for (int j = 0; j < Ny; j++) {
//            ubc(i + NxGhost, j + NyGhost) = uIC(i, j);
//        }
//    }
//
//    switch(xBC)
//    {
//        case(Neumann):
//        {
//            for (int j = 0; j < Ny; j++) {
//                for (int i = 0; i < NxGhost; i++) {
//                    ubc(i, j + NyGhost) = uIC(0, j);
//                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(Nx - 1, j);
//                }
//            }
//            break;
//        }
//        case(Periodic):
//        {
//            /// Periodic BC in x-direction
//            for (int j = 0; j < Ny; j++) {
//                for (int i = 0; i < NxGhost; i++) {
//                    ubc(i, j + NyGhost) = uIC(Nx - NxGhost + i, j);
//                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(NxGhost + i, j);
//                }
//            }
//            break;
//        }
//        default:
//        {
//            std::cerr << "Boundary conditions not recognized.";
//            std::exit(EXIT_FAILURE);
//        }
//    }
//
//    switch(yBC)
//    {
//        case(Neumann):
//        {
//            /// Neumann boundary conditions
//            for (int i = 0; i < Nx; i++) {
//                for (int j = 0; j < NyGhost; j++) {
//                    ubc(i + NxGhost, j) = uIC(i, 0);
//                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, Ny - 1);
//                }
//            }
//            break;
//        }
//        case(Periodic):
//        {
//            /// Implement here other boundary conditions ... (reflective for example)
//            /// Periodic BC in y-direction
//            for (int i = 0; i < Nx; i++) {
//                for (int j = 0; j < NyGhost; j++) {
//                    ubc(i + NxGhost, j) = uIC(i, Ny - NyGhost + j);
//                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, NyGhost + j);
//                }
//            }
//            break;
//        }
//        case(Reflective):
//        {
//            /// Reflective BC in y-direction
//            for (int i = 0; i < Nx; i++) {
//                for (int j = 0; j < NyGhost; j++) {
//
//                    for (int k = 0; k < nVar; k++) {
//                        /// Reflecting u_y and B_y
//                        if (k == 2 || k == 6) {
//                            ubc(i + NxGhost, j, k) = -uIC(i, NyGhost - 1 - j, k);
//                            ubc(i + NxGhost, Ny + NyGhost + j, k) = -uIC(i, Ny - 1 - j, k);
//                        } else {
//                            ubc(i + NxGhost, j, k) = uIC(i, NyGhost - 1 - j, k);
//                            ubc(i + NxGhost, Ny + NyGhost + j, k) = uIC(i, Ny - 1 - j, k);
//                        }
//                    }
//
//                }
//            }
//            break;
//        }
//        default:
//        {
//            std::cerr << "Boundary conditions not recognized.";
//            std::exit(EXIT_FAILURE);
//        }
//    }
//}

/// 1D Implementation
//void GhostFluidMethod::ComputeGhostFluidBoundaries1D(const unsigned int& interfaceInt, const unsigned int& count, VectA4& u1, VectA4& u2) {
//
//    A4Double w1, w2, tmpW;
//
//    if (count == 0) {
//
//        w1 = convertVars1.conservativeToprimitive(u1(interfaceInt, 0));
//        double p1 = w1[3];
//        double rho1 = w1[0];
//
//        w2 = convertVars2.conservativeToprimitive(u2(interfaceInt+1, 0));
//        double p2 = w2[3];
//        double rho2 = w2[0];
//
//        double C1 = p1/pow(rho1, GammaL);
//        double C2 = p2/pow(rho2, GammaR);
//
//        for (int i = interfaceInt - 2; i <= interfaceInt; i++) {
//            tmpW = convertVars1.conservativeToprimitive(u1(i, 0));
//            tmpW[0] = exp(log(tmpW[3]/C2)/GammaR);
//            u2(i, 0) = convertVars2.primitiveToconservative(tmpW);
//        }
//
//        for (int i = interfaceInt+1; i <= interfaceInt+3; i++) {
//            tmpW = convertVars2.conservativeToprimitive(u2(i, 0));
//            tmpW[0] = exp(log(tmpW[3]/C1)/GammaL);
//            u1(i, 0) = convertVars1.primitiveToconservative(tmpW);
//        }
//
//    } else if (count == 1) {
//
//        w1 = convertVars1.conservativeToprimitive(u1(interfaceInt+1, 0));
//        double p1 = w1[3];
//        double rho1 = w1[0];
//
//        w2 = convertVars2.conservativeToprimitive(u2(interfaceInt, 0));
//        double p2 = w2[3];
//        double rho2 = w2[0];
//
//        double C1 = p1/pow(rho1, GammaL);
//        double C2 = p2/pow(rho2, GammaR);
//
//        for (int i = interfaceInt - 2; i <= interfaceInt; i++) {
//            tmpW = convertVars2.conservativeToprimitive(u2(i, 0));
//            tmpW[0] = exp(log(tmpW[3]/C1)/GammaL);
//            u1(i, 0) = convertVars1.primitiveToconservative(tmpW);
//        }
//
//        for (int i = interfaceInt+1; i <= interfaceInt+3; i++) {
//            tmpW = convertVars1.conservativeToprimitive(u1(i, 0));
//            tmpW[0] = exp(log(tmpW[3]/C2)/GammaR);
//            u2(i, 0) = convertVars2.primitiveToconservative(tmpW);
//        }
//
//    }
//}

//void GhostFluidMethod::ComputeRiemannGhostFluidBoundaries1D(const VectPair& interfaceInt, VectA4& u1, VectA4& u2) const {
//
//    for (int i = 0; i <= interfaceInt[0].first + NxGhost; i++) {
//        u2(i, NyGhost) = u2(interfaceInt[0].first + NxGhost + 1, NyGhost);
//    }
//
//    for (int i = interfaceInt[0].first + NxGhost + 1; i < Nx + 2 * NxGhost; i++) {
//        u1(i, NyGhost) = u1(interfaceInt[0].first + NxGhost, NyGhost);
//    }
//
//}

/// Other codes
//    double phix, phiy, sqrttrm, phinew;
//    double a, b, c;
//
//    for (int i = 0; i < Nx; i++) {
//        for (int j = 0; j < Ny; j++) {
//
//            bool isInterfaceCell = false;
//
//            for (const auto& cell : interfaceInt) {
//                if (i == cell.first && j == cell.second) isInterfaceCell = true;
//            }
//
//            if (!isInterfaceCell) {
//                if (phi[i + j*Nx] > 0) {
//                    phi[i + j*Nx] = 1e100;
//                    phix = std::min<double>(phi[i+1 + j*Nx], phi[i-1 + j*Nx]);
//                    phiy = std::min<double>(phi[i + (j+1)*Nx], phi[i + (j-1)*Nx]);
//                    a = dx*dx + dy*dy;
//                    b = 2*phix*dy*dy + 2*phiy*dx*dx;
//                    c = phix*phix*dy*dy + phiy*phiy*dx*dx - dx*dx*dy*dy;
//                    sqrttrm = b*b - 4*a*c;
//                    if (sqrttrm < 0) {
//                        b = 2*phix*dy*dy;
//                        c = phix*phix*dy*dy - dx*dx*dy*dy;
//                        sqrttrm = b*b - 4*a*c;
//                    }
//                    phinew = (-b + sqrt(sqrttrm))/(2*a);
//                    if (phinew < phi[i + j*Nx]) phi[i + j*Nx] = phinew;
//                } else if (phi[i + j*Nx] < 0) {
//                    phi[i + j*Nx] = -1e100;
//                    phix = std::max<double>(phi[i+1 + j*Nx], phi[i-1 + j*Nx]);
//                    phiy = std::max<double>(phi[i + (j+1)*Nx], phi[i + (j-1)*Nx]);
//                    a = dx*dx + dy*dy;
//                    b = 2*phix*dy*dy + 2*phiy*dx*dx;
//                    c = phix*phix*dy*dy +phiy*phiy*dx*dx - dx*dx*dy*dy;
//                    sqrttrm = b*b - 4*a*c;
//                    if (sqrttrm < 0) {
//                        b = 2*phix*dy*dy;
//                        c = phix*phix*dy*dy - dx*dx*dy*dy;
//                        sqrttrm = b*b - 4*a*c;
//                    }
//                    phinew = (-b - sqrt(sqrttrm))/(2*a);
//                    if (phinew > phi[i + j*Nx]) phi[i + j*Nx] = phinew;
//                }
//            }
//        }
//    }

//void GhostFluidMethod::IterativeMethod(const int& N, const VDouble& phibc, VectA4& w1, VectA4& w2, VectA4& wbc1, VectA4& wbc2) {
//
//    A4Double deltax, deltay;
//    double dt = 0.5 * std::min(dx, dy);
//    std::array<double, 2> normal{0, 0};
//
//    for (int i = NxGhost; i < Nx + NxGhost; i++) {
//        for (int j = NyGhost; j < Ny + NyGhost; j++) {
//            w1(i-NxGhost, j-NyGhost) = wbc1(i, j);
//            w2(i-NxGhost, j-NyGhost) = wbc2(i, j);
//        }
//    }
//
//    /// Treat wbc as current and w as updated
//
//    for (int sweep = 0; sweep < N; sweep++) {
//
//        for (int i = NxGhost; i < Nx + NxGhost; i++) {
//            for (int j = NyGhost; j < Ny + NyGhost; j++) {
//                wbc1(i, j) = w1(i-NxGhost, j-NyGhost);
//                wbc2(i, j) = w2(i-NxGhost, j-NyGhost);
//            }
//        }
//
//        for (int i = NxGhost; i < Nx + NxGhost; i++) {
//            for (int j = NyGhost; j < Ny + NyGhost; j++) {
//
//                normal = GetNormalVector(phibc, i ,j);
//
//                double mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
//
//                if (mag < 1e-12) {
//                    normal[0] = 1.0;
//                    normal[1] = 0.0;
//                } else {
//                    normal[0] /= mag;
//                    normal[1] /= mag;
//                }
//
//                if (phibc[i + j*(Nx + 2 * NxGhost)] < 0.0) {  // Material 2 is ghost
//
//                    if (normal[0] > 0.0) deltax = LoopVariables(wbc2(i+1, j), wbc2(i, j));
//                    else deltax = LoopVariables(wbc2(i, j), wbc2(i-1, j));
//                    if (normal[1] > 0.0) deltay = LoopVariables(wbc2(i, j+1), wbc2(i, j));
//                    else deltay = LoopVariables(wbc2(i, j), wbc2(i, j-1));
//
//                    for (int n = 0; n < nVar; n++) {
//                        w2(i - NxGhost, j - NyGhost, n) = wbc2(i, j, n) + (dt/dx) * normal[0] * deltax[n] + (dt/dy) * normal[1] * deltay[n];
//                    }
//
//                } else {
//
//                    if (normal[0] > 0.0) deltax = LoopVariables(wbc1(i, j), wbc1(i-1, j));
//                    else deltax = LoopVariables(wbc1(i+1, j), wbc1(i, j));
//                    if (normal[1] > 0.0) deltay = LoopVariables(wbc1(i, j), wbc1(i, j-1));
//                    else deltay = LoopVariables(wbc1(i, j+1), wbc1(i, j));
//
//                    for (int n = 0; n < nVar; n++) {
//                        w1(i - NxGhost, j - NyGhost, n) = wbc1(i, j, n) - (dt/dx) * normal[0] * deltax[n] - (dt/dy) * normal[1] * deltay[n];
//                    }
//
//                }
//            }
//        }
//
//        SetBoundaryConditions(w1, wbc1);
//        SetBoundaryConditions(w2, wbc2);
//    }
//
//    for (int i = NxGhost; i < Nx + NxGhost; i++) {
//        for (int j = NyGhost; j < Ny + NyGhost; j++) {
//            wbc1(i, j) = w1(i-NxGhost, j-NyGhost);
//            wbc2(i, j) = w2(i-NxGhost, j-NyGhost);
//        }
//    }
//
//}
//
// Function to loop through conservative/primitive variable array to get difference
//A4Double GhostFluidMethod::LoopVariables(const A4Double& x1, const A4Double& x2) const {
//
//    A4Double tmp;
//
//    for (int n = 0; n < nVar; n++) {
//        tmp[n] = x1[n] - x2[n];
//    }
//
//    return tmp;
//}