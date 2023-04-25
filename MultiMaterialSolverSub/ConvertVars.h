#ifndef CONVERTVARS_H
#define CONVERTVARS_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "EquationofStates.h"
#include "TypeDefs.h"

using namespace std;

class ConvertVars {

public:

    ConvertVars();
    ~ConvertVars();

    void SetVariables(unsigned int nvar, double gamma, double pinf);
    A4Double conservativeToprimitive(const A4Double& U);
    A4Double primitiveToconservative(const A4Double& W);
    A4Double EulerFlux(const A4Double& U, int dim);

private:

    /// Classes required
    EquationofStates equationofStates;

    unsigned int nVar;
    double Gamma;
    double Pinf;

    A4Double Q;
    A4Double w;
    A4Double F;
    A4Double Wtemp;

};

#endif