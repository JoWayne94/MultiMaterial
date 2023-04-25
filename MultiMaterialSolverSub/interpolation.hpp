#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#pragma once
#include "TypeDefs.h"
#include <vector>
#include <array>

using namespace std;

template <typename T>
inline T linear_interpolate
(
	const T& fL,
	const T& fR,
	const double alpha,
    const unsigned int& nVar
)
{
    T tmp;

    for (int n = 0; n < nVar; n++) {
        tmp[n] = fL[n] + alpha * (fR[n] - fL[n]);
    }

	return tmp;
}

template <typename T>
inline T bilinear_interpolate
(
	const T& fBL,
	const T& fBR,
	const T& fTL,
	const T& fTR,
	const double alpha,
	const double beta,
    const unsigned int& nVar
)
{
	T fB = linear_interpolate<T>(fBL, fBR, alpha, nVar);
	T fT = linear_interpolate<T>(fTL, fTR, alpha, nVar);
	return linear_interpolate<T>(fB, fT, beta, nVar);
}

template <typename T>
inline T grid_bilinear_interpolate
(
	const double& dx, const double& dy, const VDouble& xCellCentres, const VDouble& yCellCentres, const double& xLeftDomain, const double& yBottomDomain,
    const unsigned int& Nx, const unsigned int& nVar, const std::vector<T>& grid, const std::array<double, 2>& pos
)
{
	std::array<double, 2> displaced_pos{0, 0};

	displaced_pos[0] = pos[0] + 0.5 * dx;
	displaced_pos[1] = pos[1] + 0.5 * dy;

	int i = std::floor((displaced_pos[0] - xLeftDomain) / dx) + 2;  // NGhost
	int j = std::floor((displaced_pos[1] - yBottomDomain) / dy) + 2;
    double alpha = (pos[0] - xCellCentres[i-3 + (j-3)*Nx]) / dx;
	double beta = (pos[1] - yCellCentres[i-3 + (j-3)*Nx]) / dy;

	return bilinear_interpolate<T>(grid[i-1 + (j-1)*(Nx + 2 * 2)], grid[i + (j-1)*(Nx + 2 * 2)], grid[(i-1) + j*(Nx + 2 * 2)], grid[i + j*(Nx + 2 * 2)], alpha, beta, nVar);
}

#endif