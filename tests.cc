
#define _USE_MATH_DEFINES
#include <cmath>
#include "PeriodicArray.h"
#include "RK4Register.h"
#include <iostream>
#include <string>

#define NX 10
#define NY 20
#define NZ 30
#define LX 15.0
#define LY 20.0
#define LZ 25.0
#define DT 0.1

int main(int argc, char **argv)
{
  int i, j, k;

  // Array tests

  PeriodicArray<int, double> * arr = new PeriodicArray<int, double>(NX, NY, NZ, LX, LY, LZ);

  for(i=0; i<NX; ++i)
    for(j=0; j<NY; ++j)
      for(k=0; k<NZ; ++k)
      {
        (*arr)(i,j,k) = std::sin(2.0*M_PI*j/NY);
      }

  std::cout << "Array parameters nx, ny, nz, lx, ly, lz, dx, dy, dz are:"
    << std::endl << "  "
    << arr->nx << ", " << arr->ny << ", " << arr->nz << ", "
    << arr->lx << ", " << arr->ly << ", " << arr->lz << ", "
    << arr->dx << ", " << arr->dy << ", " << arr->dz << std::endl;

  std::cout << "Array values along y-axis are:" << std::endl;
  std::cout << "  [" << (*arr)(0,0,0);
  for(j=1; j<NY; ++j)
    std::cout << ", " << (*arr)(0,j,0);
  std::cout << "]" << std::endl;


  // RK4 register tests

  RK4Register<int, double> * rk4 = new RK4Register<int, double>(NX, NY, NZ, LX, LY, LZ, DT);

  for(i=0; i<NX; ++i)
    for(j=0; j<NY; ++j)
      for(k=0; k<NZ; ++k)
      {
        (*rk4)(i,j,k) = std::sin(2.0*M_PI*j/NY);
        rk4->_p(i,j,k) = std::sin(2.0*M_PI*j/NY);
        rk4->_a(i,j,k) = std::sin(2.0*M_PI*j/NY);
        rk4->_c(i,j,k) = std::sin(2.0*M_PI*j/NY);
        rk4->_f(i,j,k) = std::sin(2.0*M_PI*j/NY);
      }

  std::cout << "RK4 values along y-axis are:" << std::endl;
  std::cout << "  [" << (*rk4)(0,0,0);
  for(j=1; j<NY; ++j)
    std::cout << ", " << (*rk4)(0,j,0);
  std::cout << "]" << std::endl;

  std::cout << "An interpolated value between " << (*rk4)(0,1,0) << " and "
    << (*rk4)(0,2,0) << " is: " <<
    rk4->getInterpolatedValue(0.0, 1.5, 0.0) << std::endl;

  std::cout << "Running RK routines...";
  rk4->stepInit();
  rk4->K1Finalize();
  rk4->K2Finalize();
  rk4->K3Finalize();
  rk4->K4Finalize();
  std::cout << " done." << std::endl;

  return 0;
}
