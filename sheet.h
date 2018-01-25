/** @file sheet.h
 * @brief Class implementing a "sheet" / interpolated phase-space
 * method for solving the Vlasov equation in Newtonian gravity.
 *  Vocabulary:
 *   - Phase-space: Reference to 6-d (x, y, z, vx, vy, vz) space in which
 *     matter lives.
 *   - Sheet: 3-dimensional slice through 6-d phase-space on which all matter
 *      content lives.
 *   - Sheet coordinates: coordinates parameterizing the 3-d slice (s1, s2, s3)
 *   - Metric coordinates: x, y, z
 *   - Position fields: \vec{x}(\vec{s})
 *
 *  Two coordinate systems:
 *  - Sheet coordinates
 *    position and velocity fields live on these coordinates
 *  - Metric coordinates
 *    Density and metric (Newtonian potential)
 */

#ifndef SHEET
#define SHEET

#include "RK4Register.h"
#include "PeriodicArray.h"
#include "Fourier.h"
#include "Timer.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include "zlib.h"
#include <limits>
#include <sstream>
#include <iomanip>
#include <fstream>

/**
 * Class used to run a sheet sim.
 */
class SheetSimulation
{
public:

  // internal types
  typedef long long idx_t;
  typedef double real_t;
  typedef RK4Register<idx_t, real_t> RK4_t;
  typedef PeriodicArray<idx_t, real_t> array_t;
  typedef Fourier<idx_t, real_t> fourier_t;

  TimerManager _timer;
  idx_t step;

  enum Verbosity { none, minimal, verbose, debug };
  Verbosity verbosity;

  enum depositScheme { CIC, PCS };
  enum carrierCountScheme { per_dx, per_ds };
  enum initializationType { uniform1d, uniform1dv, overdensity1d };

  // Simulation information
  struct Specs
  {
    idx_t nx, ny, nz; ///< Metric grid resolution
    idx_t ns1, ns2, ns3; ///< Phase-space sheet resolution
    real_t lx, ly, lz; ///< Metric grid physical dimensions

    idx_t carriers_per_dx,
          carriers_per_dy,
          carriers_per_dz;

    real_t dt; ///< timestep

    depositScheme deposit;
    carrierCountScheme carrier_count_scheme;
    initializationType initialization_type;
  };

  Specs specs = {0};
  idx_t s1_dbg = 0, s2_dbg = 0, s3_dbg = 0;

  // Phase-space fields, parametrized by sheet coordinates
  RK4_t *Dx, *Dy, *Dz; ///< Phase-space displacement fields
  RK4_t *vx, *vy, *vz; ///< Phase-space velocity fields

  // Metric-space fields
  array_t *rho; ///< Metric-space density
  array_t *dx_phi, *dy_phi, *dz_phi; ///< Metric-space derivatives of metric
                                     ///< (Newtonian potential)
  array_t *phi; ///< Gravitational potential

  fourier_t * fourierX; ///< FFT in metric space

  /**
   * Wrapper to print TimerManager information
   */
  void printTimingInformation()
  {
    std::cout << _timer;
  }

  /**
   * Functions to convert s-indices to non-displaced coordinates
   */
  real_t _S1IDXtoX0(idx_t s1) { return s1*specs.lx/specs.ns1; }
  real_t _S2IDXtoY0(idx_t s2) { return s2*specs.ly/specs.ns2; }
  real_t _S3IDXtoZ0(idx_t s3) { return s3*specs.lz/specs.ns3; }
  real_t _S1IDXtoX0(real_t s1) { return s1*specs.lx/specs.ns1; }
  real_t _S2IDXtoY0(real_t s2) { return s2*specs.ly/specs.ns2; }
  real_t _S3IDXtoZ0(real_t s3) { return s3*specs.lz/specs.ns3; }

  void _MassDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
    bool announce)
  {
    switch (specs.deposit)
    {
      case PCS:
        _PCSDeposit(weight, x_idx, y_idx, z_idx, announce);
        break;
      case CIC:
      default:
        _CICDeposit(weight, x_idx, y_idx, z_idx, announce);
        break;
    }
  }

  void _CICDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
    bool announce)
  {
    idx_t ix, iy, iz,
          ixp, iyp, izp;

    real_t x_f, y_f, z_f,
           x_h, y_h, z_h;

    // gridpoint index "left" of x_idx
    ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );
    ixp = ix + 1;
    x_f = x_idx - (real_t) ix;
    x_h = 1.0 - x_f;

    iy = (y_idx < 0 ? (idx_t) y_idx - 1 : (idx_t) y_idx );
    iyp = iy + 1;
    y_f = y_idx - (real_t) iy;
    y_h = 1.0 - y_f;

    iz = (z_idx < 0 ? (idx_t) z_idx - 1 : (idx_t) z_idx );
    izp = iz + 1;
    z_f = z_idx - (real_t) iz;
    z_h = 1.0 - z_f;

    if(announce)
    {
      std::cout << "  depositing near " << ix << "," << iy << "," << iz << "; ";
      std::cout << "weight = " << weight << ", x_h=" << x_h << ", x_f=" << x_f << "; sum="
      << x_h*y_h*z_h
         +x_h*y_h*z_f
         +x_h*y_f*z_h
         +x_h*y_f*z_f
         +x_f*y_h*z_h
         +x_f*y_h*z_f
         +x_f*y_f*z_h
         +x_f*y_f*z_f << std::endl;
    }

    (*rho)(ix, iy, iz) += x_h*y_h*z_h*weight;
    (*rho)(ix, iy, izp) += x_h*y_h*z_f*weight;
    (*rho)(ix, iyp, iz) += x_h*y_f*z_h*weight;
    (*rho)(ix, iyp, izp) += x_h*y_f*z_f*weight;

    (*rho)(ixp, iy, iz) += x_f*y_h*z_h*weight;
    (*rho)(ixp, iy, izp) += x_f*y_h*z_f*weight;
    (*rho)(ixp, iyp, iz) += x_f*y_f*z_h*weight;
    (*rho)(ixp, iyp, izp) += x_f*y_f*z_f*weight;
  }

  void _PCSDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
    bool announce)
  {
    idx_t ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );
    idx_t iy = (y_idx < 0 ? (idx_t) y_idx - 1 : (idx_t) y_idx );
    idx_t iz = (z_idx < 0 ? (idx_t) z_idx - 1 : (idx_t) z_idx );

    real_t norm = 0.0;
    for(idx_t i=-1; i<=2; ++i)
      for(idx_t j=-1; j<=2; ++j)
        for(idx_t k=-1; k<=2; ++k)
        {
          real_t s = std::sqrt( std::pow(ix+i-x_idx, 2)
            + std::pow(iy+j-y_idx, 2)
            + std::pow(iz+k-z_idx, 2) );
          
          if(s<1.0)
          {
            norm += (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
          }
          else if(s<2.0 && s>=1.0)
          {
            norm += std::pow(2.0 - std::abs(s), 3)/6.0;
          }
        }

    real_t pcs;
    for(idx_t i=-1; i<=2; ++i)
      for(idx_t j=-1; j<=2; ++j)
        for(idx_t k=-1; k<=2; ++k)
        {
          real_t s = std::sqrt( std::pow(ix+i-x_idx, 2)
            + std::pow(iy+j-y_idx, 2)
            + std::pow(iz+k-z_idx, 2) );
          
          if(s<1.0)
          {
            pcs = (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
            (*rho)(ix+i, iy+j, iz+k) += pcs*weight/norm;
          }
          else if(s<2.0 && s>=1.0)
          {
            pcs = std::pow(2.0 - std::abs(s), 3)/6.0;
            (*rho)(ix+i, iy+j, iz+k) += pcs*weight/norm;
          }
        }
  }

  /**
   * Get Min/max x/y/z coordinates at voxel corners
   */
  real_t _getXRangeInSVoxel(RK4_t * DX, idx_t s1_idx, idx_t s2_idx,
    idx_t s3_idx, real_t X0_lower, real_t X0_upper)
  {
    std::vector<real_t> X_vals_on_bounds = {
      (*DX)(s1_idx,s2_idx,s3_idx) + X0_lower,
      (*DX)(s1_idx+1,s2_idx,s3_idx) + X0_upper,
      (*DX)(s1_idx,s2_idx+1,s3_idx) + X0_lower,
      (*DX)(s1_idx,s2_idx,s3_idx+1) + X0_lower,
      (*DX)(s1_idx+1,s2_idx+1,s3_idx) + X0_upper,
      (*DX)(s1_idx+1,s2_idx,s3_idx+1) + X0_upper,
      (*DX)(s1_idx,s2_idx+1,s3_idx+1) + X0_lower,
      (*DX)(s1_idx+1,s2_idx+1,s3_idx+1) + X0_upper
    };

    real_t X_min = std::numeric_limits<real_t>::max();
    real_t X_max = std::numeric_limits<real_t>::min();
    for(int i=0; i<8; ++i)
    {
      if(X_vals_on_bounds[i] > X_max) X_max = X_vals_on_bounds[i];
      if(X_vals_on_bounds[i] < X_min) X_min = X_vals_on_bounds[i];
    }

    return X_max - X_min;
  }

  /**
   * Compute conribution to rho(x) from data in a phase-space
   * sheet voxel and add to rho(x) grid. Do so via 1501.01959 mass
   * deposition scheme.
   * TODO: improve; consider higher-order or analytic versions of this
   */
  void _pushSheetMassToRho(idx_t s1, idx_t s2, idx_t s3)
  {
    _timer["_pushSheetMassToRho"].start();

    idx_t num_x_carriers, num_y_carriers, num_z_carriers;

    switch(specs.carrier_count_scheme)
    {
      case per_dx:
        num_x_carriers = specs.carriers_per_dx == 0 ? 1 : specs.carriers_per_dx*(
          (idx_t) (0.5 + _getXRangeInSVoxel(Dx, s1, s2, s3, _S1IDXtoX0(s1), _S1IDXtoX0(s1+1)) / rho->dx ) );
        num_y_carriers = specs.carriers_per_dy == 0 ? 1 : specs.carriers_per_dy*(
          (idx_t) (0.5 + _getXRangeInSVoxel(Dy, s1, s2, s3, _S2IDXtoY0(s2), _S2IDXtoY0(s2+1)) / rho->dy ) );
        num_z_carriers = specs.carriers_per_dz == 0 ? 1 : specs.carriers_per_dz*(
          (idx_t) (0.5 + _getXRangeInSVoxel(Dz, s1, s2, s3, _S3IDXtoZ0(s3), _S3IDXtoZ0(s3+1)) / rho->dz ) );
        break;
      case per_ds:
      default:
        num_x_carriers = specs.carriers_per_dx;
        num_y_carriers = specs.carriers_per_dy;
        num_z_carriers = specs.carriers_per_dz;
        break;
    }

    if(num_x_carriers <= 0) num_x_carriers = 1;
    if(num_y_carriers <= 0) num_y_carriers = 1;
    if(num_z_carriers <= 0) num_z_carriers = 1;

    idx_t num_carriers = num_x_carriers*num_y_carriers*num_z_carriers;
    if(verbosity == debug && s1==s1_dbg && s2==s2_dbg && s3==s3_dbg)
    {
      std::cout << "  Performing projection for debug point using "
        << num_carriers << " carriers." << std::endl;
    }

    real_t weight = 1.0 / (real_t) num_carriers
      * specs.nx*specs.ny*specs.nz / specs.ns1/specs.ns2/specs.ns3;

    // distribute mass from all carriers
    idx_t i, j, k;
    for(i=0; i<num_x_carriers; ++i)
      for(j=0; j<num_y_carriers; ++j)
        for(k=0; k<num_z_carriers; ++k)
        {
          real_t carrier_s1 = (real_t) s1 + (real_t) i / (real_t) num_x_carriers;
          real_t carrier_s2 = (real_t) s2 + (real_t) j / (real_t) num_y_carriers;
          real_t carrier_s3 = (real_t) s3 + (real_t) k / (real_t) num_z_carriers;

          real_t carrier_x_idx = ( _S1IDXtoX0(carrier_s1)
            + Dx->getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3) ) / rho->dx;
          real_t carrier_y_idx = ( _S2IDXtoY0(carrier_s2)
            + Dy->getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3) ) / rho->dy;
          real_t carrier_z_idx = ( _S3IDXtoZ0(carrier_s3)
            + Dz->getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3) ) / rho->dz;

          bool announce = false;
          if(verbosity == debug && /*s1==s1_dbg &&*/ s2==s2_dbg && s3==s3_dbg)
          {
            std::cout << "  depositing carrier from (" << carrier_s1 << ","
              << carrier_s2 << "," << carrier_s3
              << ") at (" << carrier_x_idx << ","
              << carrier_y_idx << "," << carrier_z_idx << ")." << std::endl;
            announce = true;
          }

          _MassDeposit(weight, carrier_x_idx, carrier_y_idx, carrier_z_idx, announce);
        }

    _timer["_pushSheetMassToRho"].stop();
  }

  /**
   * Set metric potentials (or really just derivatives thereof)
   */
  void _setMetricPotentials()
  {
    _timer["_setMetricPotentials"].start();
    if(verbosity == debug)
    {
      std::cout << "  Setting Metric potentials...";
    }

    // Just set \phi'(x)  ~ 1/k \rho
    // TODO: consider large-scale correction terms?
    *dx_phi = *rho;
    *dy_phi = *rho;
    *dz_phi = *rho;
    *phi = *rho;

    fourierX->inverseGradient(&(*dx_phi)[0], 1); // &(*dx_phi)[0] is pointer to dx_phi internal array
    fourierX->inverseGradient(&(*dy_phi)[0], 2);
    fourierX->inverseGradient(&(*dz_phi)[0], 3);
    fourierX->inverseLaplacian(&(*phi)[0]);
    if(verbosity == debug)
    {
      std::cout << " done." << std::endl;
    }
    _timer["_setMetricPotentials"].stop();
  }

  /**
   * Intermediate RK4 calculations for phase space fields
   */
  void _RK4Calc()
  {
    _timer["_RK4Calc"].start();
    idx_t i ,j, k;

    // density projection
    // clear rho
    if(verbosity == debug) std::cout << "  Performing density projection..." << std::flush;
    for(idx_t p=0; p<rho->nx*rho->ny*rho->nz; ++p)
      (*rho)[p] = 0;
    for(i=0; i<specs.ns1; ++i)
      for(j=0; j<specs.ns2; ++j)
        for(k=0; k<specs.ns3; ++k)
          _pushSheetMassToRho(i, j, k);
    if(verbosity == debug) std::cout << " done." << std::endl << std::flush;

    // Metric potential from density
    _setMetricPotentials();

    // RK4 calculation
    if(verbosity == debug) std::cout << "  Performing RK4 calculation..." << std::flush;
    for(i=0; i<specs.ns1; ++i)
      for(j=0; j<specs.ns2; ++j)
        for(k=0; k<specs.ns3; ++k)
        {
          Dx->_c(i,j,k) = vx->_a(i,j,k);
          Dy->_c(i,j,k) = vy->_a(i,j,k);
          Dz->_c(i,j,k) = vz->_a(i,j,k);

          real_t x_pt = Dx->_a(i,j,k) + _S1IDXtoX0(i);
          real_t y_pt = Dy->_a(i,j,k) + _S2IDXtoY0(j);
          real_t z_pt = Dz->_a(i,j,k) + _S3IDXtoZ0(k);
          vx->_c(i,j,k) = -1.0*dx_phi->getTriCubicInterpolatedValue(
            x_pt/dx_phi->dx, y_pt/dx_phi->dy, z_pt/dx_phi->dz);
          vy->_c(i,j,k) = -1.0*dy_phi->getTriCubicInterpolatedValue(
            x_pt/dy_phi->dx, y_pt/dy_phi->dy, z_pt/dy_phi->dz);
          vz->_c(i,j,k) = -1.0*dz_phi->getTriCubicInterpolatedValue(
            x_pt/dz_phi->dx, y_pt/dz_phi->dy, z_pt/dz_phi->dz);
        }
    if(verbosity == debug) std::cout << " done." << std::endl << std::flush;
    _timer["_RK4Calc"].stop();
  }

  void _stepInit()
  {
    Dx->stepInit();
    Dy->stepInit();
    Dz->stepInit();
    vx->stepInit();
    vy->stepInit();
    vz->stepInit();
  }

  void _K1Finalize()
  {
    Dx->K1Finalize();
    Dy->K1Finalize();
    Dz->K1Finalize();
    vx->K1Finalize();
    vy->K1Finalize();
    vz->K1Finalize();
  }

  void _K2Finalize()
  {
    Dx->K2Finalize();
    Dy->K2Finalize();
    Dz->K2Finalize();
    vx->K2Finalize();
    vy->K2Finalize();
    vz->K2Finalize();
  }

  void _K3Finalize()
  {
    Dx->K3Finalize();
    Dy->K3Finalize();
    Dz->K3Finalize();
    vx->K3Finalize();
    vy->K3Finalize();
    vz->K3Finalize();
  }

  void _K4Finalize()
  {
    Dx->K4Finalize();
    Dy->K4Finalize();
    Dz->K4Finalize();
    vx->K4Finalize();
    vy->K4Finalize();
    vz->K4Finalize();
  }

  void RKStep()
  {
    if(verbosity == debug) std::cout << "Initializing RK step." << std::endl << std::flush;
    _stepInit();
    
    if(verbosity == debug) std::cout << "Performing K1 step." << std::endl << std::flush;
    _RK4Calc();
    if(verbosity == debug) std::cout << "  Finalizing." << std::endl << std::flush;
    _K1Finalize();

    if(verbosity == debug) std::cout << "Performing K2 step." << std::endl << std::flush;
    _RK4Calc();
    if(verbosity == debug) std::cout << "  Finalizing." << std::endl << std::flush;
    _K2Finalize();

    if(verbosity == debug) std::cout << "Performing K3 step." << std::endl << std::flush;
    _RK4Calc();
    if(verbosity == debug) std::cout << "  Finalizing." << std::endl << std::flush;
    _K3Finalize();

    if(verbosity == debug) std::cout << "Performing K4 step." << std::endl << std::flush;
    _RK4Calc();
    if(verbosity == debug) std::cout << "  Finalizing." << std::endl << std::flush;
    _K4Finalize();

    step++;
  }

  void initialize(const Specs specs_in)
  {
    specs = specs_in;

    // Phase-space variables
    Dx = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);
    Dy = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);
    Dz = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);
    vx = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);
    vy = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);
    vz = new RK4_t(specs.ns1, specs.ns2, specs.ns3,
      specs.lx, specs.ly, specs.lz, specs.dt);

    // Metric-space fields
    rho = new array_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz);
    dx_phi = new array_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz);
    dy_phi = new array_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz);
    dz_phi = new array_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz);
    phi = new array_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz);

    // For computing Fourier transforms
    fourierX = new fourier_t(specs.nx, specs.ny, specs.nz,
      specs.lx, specs.ly, specs.lz, &(*rho)[0]);
  }

  SheetSimulation(const Specs specs_in)
  {
    verbosity = none;
    step = 0;
    initialize(specs_in);
  }

  SheetSimulation(const Specs specs_in, Verbosity verbosity_in)
  {
    verbosity = verbosity_in;
    initialize(specs_in);
  }

  ~SheetSimulation()
  {
    delete Dx;
    delete Dy;
    delete Dz;

    delete vx;
    delete vy;
    delete vz;

    delete rho;
    delete dx_phi;
    delete dy_phi;
    delete dz_phi;
    delete phi;
    delete fourierX;
  }

  void initializeFields()
  {
    if(verbosity > none)
      std::cout << "Initializing fields..." << std::endl;

    switch(specs.initialization_type)
    {
      case uniform1d :
        _initialize1DUniform();
        break;

      case uniform1dv :
        _initialize1DUniformMoving();
        break;

      case overdensity1d :
      default :
        _initialize1DOverdensity();
        break;
    }
    
    _initialize1DUniform();

    // copy to _a registers as well,
    // set metric for output
    _stepInit();
    _RK4Calc();
  }

  void _initialize1DUniform()
  {
    // position fields should just be coordinates
    for(idx_t i=0; i<specs.ns1; i++)
      for(idx_t j=0; j<specs.ns2; j++)
        for(idx_t k=0; k<specs.ns3; k++)
        {
          Dx->_p(i,j,k) = 1.0/specs.nx/5.0;
        }
  }

  void _initialize1DUniformMoving()
  {
    // position fields should just be coordinates
    for(idx_t i=0; i<specs.ns1; i++)
      for(idx_t j=0; j<specs.ns2; j++)
        for(idx_t k=0; k<specs.ns3; k++)
        {
          Dx->_p(i,j,k) = 1.0/specs.nx/5.0;
          vx->_p(i,j,k) = 0.01;
        }
  }

  void _initialize1DOverdensity()
  {
    // position fields should just be coordinates
    for(idx_t i=0; i<specs.ns1; i++)
      for(idx_t j=0; j<specs.ns2; j++)
        for(idx_t k=0; k<specs.ns3; k++)
        {
          real_t x_frac = ((real_t) i)/specs.ns1 - 0.5;
          Dx->_p(i,j,k) = 1.0/specs.nx/5.0 - 0.5*std::exp(-1.0*std::pow(x_frac/0.05, 2))*x_frac;
            //+0.3*std::exp(-1.0*std::pow((x_frac-.2)/0.05, 2.0))*x_frac;
          Dy->_p(i,j,k) = 0.0;
          Dz->_p(i,j,k) = 0.0;
        }
  }

  std::string toStr(real_t val)
  {
    std::ostringstream out;
    out << std::setprecision(17) << val;
    return out.str();
  }

  void writeDensity(std::string filename)
  {
    _timer["writeDensity"].start();
    filename = filename + ".gz";
    gzFile fi = gzopen(filename.c_str(), "wb");

    idx_t points = specs.nx*specs.ny*specs.nz;
    for(idx_t p=0; p<points; ++p)
    {
      std::string str = toStr((*rho)[p]) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);
    _timer["writeDensity"].stop();
  }

  void writedphis(std::string filename)
  {
    _timer["writedphis"].start();
    filename = filename + ".gz";
    gzFile fi = gzopen(filename.c_str(), "wb");

    idx_t points = specs.nx*specs.ny*specs.nz;
    for(idx_t p=0; p<points; ++p)
    {
      std::string str = toStr((*dx_phi)[p]) + "\t"
        + toStr((*dy_phi)[p]) + "\t"
        + toStr((*dz_phi)[p]) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);
    _timer["writedphis"].stop();
  }

  void writePositions(std::string filename)
  {
    _timer["writePositions"].start();
    filename = filename + ".gz";
    gzFile fi = gzopen(filename.c_str(), "wb");

    idx_t points = specs.ns1*specs.ns2*specs.ns3;
    for(idx_t p=0; p<points; ++p)
    {
      std::string str = toStr((*Dx)[p]) + "\t"
        + toStr((*Dy)[p]) + "\t"
        + toStr((*Dz)[p]) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);
    _timer["writePositions"].stop();
  }

  void writeStrips(std::string filename_base)
  {
    _timer["writeStrips"].start();
    std::string filename = filename_base + "_density.strip.gz";
    gzFile fi = gzopen(filename.c_str(), "ab");
    for(idx_t p=0; p<specs.nx; ++p)
    {
      std::string str = toStr((*rho)(p,1,1)) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);

    filename = filename_base + "_x.strip.gz";
    fi = gzopen(filename.c_str(), "ab");
    for(idx_t p=0; p<specs.ns1; ++p)
    {
      std::string str = toStr((*Dx)(p,1,1) + _S1IDXtoX0(p)) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);

    filename = filename_base + "_vx.strip.gz";
    fi = gzopen(filename.c_str(), "ab");
    for(idx_t p=0; p<specs.ns1; ++p)
    {
      std::string str = toStr((*vx)(p,1,1)) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);

    filename = filename_base + "_dx_phi.strip.gz";
    fi = gzopen(filename.c_str(), "ab");
    for(idx_t p=0; p<specs.nx; ++p)
    {
      std::string str = toStr((*dx_phi)(p,1,1)) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);

    filename = filename_base + "_phi.strip.gz";
    fi = gzopen(filename.c_str(), "ab");
    for(idx_t p=0; p<specs.nx; ++p)
    {
      std::string str = toStr((*phi)(p,1,1)) + "\t";
      gzwrite(fi, str.c_str(), str.length());
    }
    gzclose(fi);
    _timer["writeStrips"].stop();
  }

  void writeInfo(std::string directory)
  {
    std::string filename (directory + "/SheetSimulation.info");

    std::ofstream ofs;
    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);

    ofs << "nx = " << specs.nx << std::endl;
    ofs << "ny = " << specs.ny << std::endl;
    ofs << "nz = " << specs.nz << std::endl;
    ofs << "ns1 = " << specs.ns1 << std::endl;
    ofs << "ns2 = " << specs.ns2 << std::endl;
    ofs << "ns3 = " << specs.ns3 << std::endl;
    ofs << "lx = " << specs.lx << std::endl;
    ofs << "ly = " << specs.ly << std::endl;
    ofs << "lz = " << specs.lz << std::endl;
    ofs << "carriers_per_dx = " << specs.carriers_per_dx << std::endl;
    ofs << "carriers_per_dy = " << specs.carriers_per_dy << std::endl;
    ofs << "carriers_per_dz = " << specs.carriers_per_dz << std::endl;
    ofs << "dt = " << specs.dt << std::endl;

    if(specs.deposit == CIC)
      ofs << "deposit = CIC " << std::endl;
    if(specs.deposit == PCS)
      ofs << "deposit = PCS " << std::endl;

    ofs.close();
  }

};

#endif // include guard
