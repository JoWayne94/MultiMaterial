#ifndef EULER_H
#define EULER_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "ConvertVars.h"
#include "FluxUpdate.h"
#include "TypeDefs.h"
#include "ExactRPSolver.h"
#include "GhostFluidMethod.h"
#include "EquationofStates.h"

using namespace std;

class Euler {

public:

    Euler();
    ~Euler();

    void SetDomainSize(double xlength, double ylength, double xLeft, double xRight, double yBottom, double yTop);
    void SetNumCells(unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost);
    void SetCourantNumber(double cfl);
    void SetMeshSize(double deltax, double deltay);
    void SetnVar(unsigned int nvar);
    void SetTestNum(unsigned int testNum);
    void SetFileName(std::string Name);
    // void SetBCType(BCtype xbc, BCtype ybc);

    /// Main solver class functions
    void Initialise();
    void Update();

    /// Static member data
    double Cfl, dx, dy;
    int Nx, Ny, NxGhost, NyGhost, N;
    unsigned int nVar, testNumber;
    double Lx, Ly, xLeftDomain, xRightDomain, yBottomDomain, yTopDomain;
    double GammaL, GammaR, PinfL, PinfR;

private:

    std::string name;
    BCtype xBC, yBC;

    VectA4 *u1 = nullptr;
    VectA4 *uBc1 = nullptr;
    VectA4 *u2 = nullptr;
    VectA4 *uBc2 = nullptr;
    VectA4 *wExact = nullptr;
    VDouble *xCells = nullptr;
    VDouble *yCells = nullptr;
    VDouble *phi = nullptr;
    VDouble *phiBc = nullptr;

    /// To be defined parameters
    double dt, dt1, dt2, T;
    double RP_xDisc, RP_yDisc;
    A4Double RP_Shock, RP_LeftState, RP_RightState;
    VectPair interfaceList;

    /// Private functions to set initial and boundary conditions
    void ConstructFVDomain(VDouble& xCellCentres, VDouble& yCellCentres) const ;
    void SetBoundaryConditions(const VectA4& uIC, VectA4& ubc) const ;
    void SetInitialConditions1(const VDouble& xCellCentres, const VDouble& yCellCentres, VectA4& uIC, double& finalt, double& xDiscontinuity, double& yDiscontinuity,
                               double& gammaL, double& gammaR, double& pinfL, double& pinfR) ;
    void SetInitialConditions2(const VDouble& xCellCentres, const VDouble& yCellCentres, VectA4& uIC, VDouble& phiIC) ;
    void GenerateData(const VDouble& xCellCentres, const VDouble& yCellCentres, const VectA4& uOut1, const VectA4& uOut2, const VDouble& phiOut, const VectA4& wExactOut) ;

    /// Classes required
    ConvertVars convertVarsL;
    ConvertVars convertVarsR;
    ExactRPSolver exactRPSolver;
    FluxUpdate fluxUpdate1;
    FluxUpdate fluxUpdate2;
    GhostFluidMethod ghostFluidMethod;
    EquationofStates equationofStatesL;
    EquationofStates equationofStatesR;
};

#endif