#ifndef GHOSTFLUIDMETHOD_H
#define GHOSTFLUIDMETHOD_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "TypeDefs.h"
#include "ConvertVars.h"
#include "ExactRPSolver.h"
#include "interpolation.hpp"

using namespace std;

class GhostFluidMethod {

public:

    GhostFluidMethod();
    ~GhostFluidMethod();

    void SetVariables(unsigned int nvar, unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost, double gammaL, double gammaR,
                      double pinfL, double pinfR, double deltax, double deltay, double xLeft, double yBot, BCtype xbc, BCtype ybc);
    VectPair FindInterfaceLocation(const VDouble& phi) const;
    void UpdateLevelSet(const VectA4& u1, const VectA4& u2, const double& dt, VDouble& phi, VDouble& phibc);
    void Reinitialisation(const VectPair& interfaceInt, VDouble& phi, VDouble& phibc) const;
    void SetGhostFluid(VectA4& u1, VectA4& u2, VectA4& ubc1, VectA4& ubc2, const VDouble& phi, VDouble& phibc, const VDouble& xCellCentres, const VDouble& yCellCentres,
                       const VectPair& interfaceInt); // phi Neumann BC here

private:

    unsigned int nVar;
    int Nx, Ny, NxGhost, NyGhost;
    double GammaL, GammaR, PinfL, PinfR, dx, dy, xLeftDomain, yBottomDomain;
    BCtype xBC, yBC;

    ConvertVars convertVars1;
    ConvertVars convertVars2;
    ExactRPSolver exactRPSolver;

    /// 1D code
    // void ComputeGhostFluidBoundaries1D(const unsigned int& interfaceInt, const unsigned int& count, VectA4& u1, VectA4& u2);
    // void ComputeRiemannGhostFluidBoundaries1D(const VectPair& interfaceInt, VectA4& u1, VectA4& u2) const;

    /// 2D Implementation
    std::array<double, 2> GetNormalVector(const VDouble& phi, int i, int j) const;
//    double GetNormalVelocity(const A4Double& w, const std::array<double, 2>& n);
//    std::array<double, 2> GetTangentialVelocity(const A4Double& w, const std::array<double, 2>& n, const double& vn);

    /// GFM
    void GetInterpolatedMixedRPStates(const VDouble& xCellCentres, const VDouble& yCellCentres, double phival, const std::array<double, 2>& normal, const VectA4& wbc1, const VectA4& wbc2, int i, int j, A4Double& LeftState, A4Double& RightState) const;
    void ComputeIntermediateStates(A4Double& wlinter, A4Double& wrinter);

    /// Fast Sweeping Method
    double FastSweepingMethodPositive(VDouble& phibc, int i, int j) const;
    double FastSweepingMethodNegative(VDouble& phibc, int i, int j) const;

    /// Extrapolating Variables
    void PopulatingGhostCells(const VectPair& interfaceInt, const VDouble& phibc, VectA4& Q1, VectA4& Q2) const ;
    A4Double FastSweepingQPos(const VectA4& Q, const VDouble& phibc, int i, int j, const std::array<double, 2>& normal) const ;
    A4Double FastSweepingQNeg(const VectA4& Q, const VDouble& phibc, int i, int j, const std::array<double, 2>& normal) const ;

//    void IterativeMethod(const int& N, const VDouble& phibc, VectA4& w1, VectA4& w2, VectA4& wbc1, VectA4& wbc2);
//    void SetBoundaryConditions(const VectA4& uIC, VectA4& ubc) const ;
//    A4Double LoopVariables(const A4Double& x1, const A4Double& x2) const;

};

#endif