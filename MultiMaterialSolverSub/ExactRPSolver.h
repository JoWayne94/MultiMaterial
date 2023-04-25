#ifndef EXACTRPSOLVER_H
#define EXACTRPSOLVER_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "TypeDefs.h"
#include "EquationofStates.h"

using namespace std;

class ExactRPSolver {

public:

    ExactRPSolver();
    ~ExactRPSolver();

    void SetVariables(const A4Double& wl, const A4Double& wr, double gammaL, double gammaR, double pinfL, double pinfR);
    void ComputeExactSolution(const double& time, const VDouble& x, const double& xd, VectA4& uExact);

    double ComputePStar();
    double ComputeUStar(const double& pstar);
    void ComputeRhoStar(const double& pstar, const double& ustar, const int& Pattern, double& rholstar, double& alstar, double& Shl, double& Stl, double& Sl, double& rhorstar, double& arstar, double& Shr, double& Str, double& Sr) const;
    void FindPattern(const double& pstar, int& Pattern) const;

private:

    unsigned int nVar;
    double rhol, rhor, ul, ur, pl, pr, al, ar;
    double GammaL, GammaR, PinfL, PinfR;
    double G1,G2,G3,G4,G5,G6,G7,G9,G10,G11,G12,G13,G14,G15;

    EquationofStates equationofStatesL, equationofStatesR;

    double GuessPressure() const;
    double fl_NR(const double& p) const;
    double fr_NR(const double& p) const;
    double f_NR(const double& p) const;
    double df_NR(const double& p) const;
    A4Double Sample(const double& pstar, const double& ustar, const int& Pattern, const double& csi) const;

};

#endif