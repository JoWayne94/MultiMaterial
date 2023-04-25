#ifndef FLUXUPDATE_H
#define FLUXUPDATE_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "NumericalFlux.h"
#include "ConvertVars.h"
#include "EquationofStates.h"
#include "TypeDefs.h"
#include <ctime>

using namespace std;

class FluxUpdate {

public:

    FluxUpdate();
    ~FluxUpdate();

    void SetVariables(unsigned int nx, unsigned int ny, double cfl, unsigned int nxGhost, unsigned int nyGhost, unsigned int nvar, double finalt, double gamma, double pinf);
    void Initialise();
    void ComputeDt(const VectA4& ubc, const double& deltax, const double& deltay, const double& time, double& deltat);
    void UpdatewithFluxes(const VectA4& ubcOld, VectA4& unew, const double& deltax, const double& deltat, int dim, PairDouble& timerPair);

private:

    /// Classes required
    NumericalFlux numericalFlux;
    ConvertVars convertVars;

    /// Timing parameters
    clock_t start, end;
    double elapsed_timehts, elapsed_timeflux;

    unsigned int Nx, Ny;
    unsigned int NxGhost, NyGhost;
    unsigned int nVar;
    double Cfl, Pinf, Gamma, T;
    VectA4 *Ql = nullptr;
    VectA4 *Qr = nullptr;

    void HalfTimeStepUpdate(const VectA4& ubc, const double& deltax, const double& deltat, A4Double& qhtsl, A4Double& qhtsr, int k, int l, int dim);

};

#endif