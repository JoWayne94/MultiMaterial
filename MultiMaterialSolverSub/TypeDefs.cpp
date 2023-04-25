#include "TypeDefs.h"

VectA4::VectA4() {
}

void VectA4::SetSize(unsigned int nx, unsigned int ny)
{
    Nx = nx;
    Ny = ny;
    N = Nx * Ny;
    m_data.resize(N);
}

A4Double& VectA4::operator()(int i, int j)
{
    return m_data[i + j*Nx];
}

const A4Double& VectA4::operator()(int i, int j) const
{
    return m_data[i + j*Nx];
}

double& VectA4::operator()(int i, int j, int k)
{
    return m_data[i + j*Nx][k];
}

const double& VectA4::operator()(int i, int j, int k) const
{
    return m_data[i + j*Nx][k];
}