#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "ConvertVars.h"
#include "EquationofStates.h"
#include "TypeDefs.h"

using namespace std;

class NumericalFlux {

public:

    NumericalFlux();
    ~NumericalFlux();

    void SetVariables(const A4Double& Ql, const A4Double& Qr, unsigned int nvar, double gamma, double pinf);
    void WaveEstimates(double& Sl, double& Sr, int dim) const;
    void HLLFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fhll, int dim);
    void HLLCFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fhllc, int dim);
    void FORCEFlux(const A4Double& Ql, const A4Double& Qr, A4Double& Fforce, const double& dx, const double& dt, int dim);

private:

    /// Classes required
    ConvertVars convertVars;
    EquationofStates equationofStates;

    unsigned int nVar;
    double Sl, Sr, Gamma, Pinf;
    A4Double wl, wr;
    double rhol, rhor, uxl, uxr, uyl, uyr, pl, pr, al, ar;

};

#endif