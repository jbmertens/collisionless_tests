#ifndef RK4REGISTER_H
#define RK4REGISTER_H

#include <string>
#include <utility>
#include <iostream>
#include "PeriodicArray.h"

/**
 * @brief RK4 Class for integration
 * 
 * @tparam IT Index type
 * @tparam RT Real type
 */
template<typename IT, typename RT>
class RK4Register : public PeriodicArray<IT, RT>
{
private:
  RT sim_dt; ///< Simulation timestep

  RT * _array_p; ///< "_p" register: contains data from _p_revious step
  RT * _array_a; ///< "_a" register: containes _a_ctive data needed for _c_omputations
  RT * _array_c; ///< "_c" register: contains _c_omputed values
  RT * _array_f; ///< "_f" register: containes final value of RK4 step

public:
  /**
   * @brief Initialize class variables
   * @details Set "dt" for class instance; grid dimensions
   * 
   * @param nx_in num. grid points in x-direction
   * @param ny_in num. grid points in y-direction
   * @param nz_in num. grid points in z-direction
   * @param lx_in physical length in x-direction
   * @param ly_in physical length in y-direction
   * @param lz_in physical length in z-direction
   * @param sim_dt_in initial timestep
   */
  RK4Register(IT nx_in, IT ny_in, IT nz_in, RT lx_in, RT ly_in, RT lz_in,
    RT sim_dt_in) :
    PeriodicArray<IT, RT>(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in)
  {
    // "active" array references PeriodicArray array
    _array_a = this->_array;

    // allocate storage for additional registers
    // TODO: vectors
    _array_p = new RT[this->_pts];
    _array_c = new RT[this->_pts];
    _array_f = new RT[this->_pts];
    // TODO: OMP simd
    // read https://stackoverflow.com/questions/14674049/parallel-for-vs-omp-simd-when-to-use-each
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      _array_p[i] = 0.0;
      _array_c[i] = 0.0;
      _array_f[i] = 0.0;
    }

    setDt(sim_dt_in);
  }

  /**
   * @brief Set "dt" for RK4Register instance
   */
  void setDt(RT sim_dt_in)
  {
    sim_dt = sim_dt_in;
  }

  ~RK4Register()
  {
    delete [] _array_p;
    delete [] _array_a;
    delete [] _array_c;
    delete [] _array_f;
  }

  /**
   * @brief Swap array (pointers) for _a and _c registers
   */
  void swap_a_c()
  {
    std::swap(_array_a, _array_c);
    this->_array = _array_a;
  }

  /**
   * @brief Swap array (pointers) for _p and _f registers
   */
  void swap_p_f()
  {
    std::swap(_array_p, _array_f);
  }

  void stepInit()
  {
    IT i;
#ifdef USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for(i=0; i<this->_pts; ++i)
    {
      _array_a[i] = _array_p[i];
      _array_f[i] = 0;
    }
  }

  void K1Finalize()
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      _array_f[i] += sim_dt*_array_c[i]/6.0;
      _array_c[i] = _array_p[i] + sim_dt*_array_c[i]/2.0;
    }

    swap_a_c();
  }

  void K2Finalize()
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      _array_f[i] += sim_dt*_array_c[i]/3.0;
      _array_c[i] = _array_p[i] + sim_dt*_array_c[i]/2.0;
    }

    swap_a_c();
  }

  void K3Finalize()
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      _array_f[i] += sim_dt*_array_c[i]/3.0;
      _array_c[i] = _array_p[i] + sim_dt*_array_c[i];
    }
    
    swap_a_c();
  }

  void K4Finalize()
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      _array_f[i] += sim_dt*_array_c[i]/6.0 + _array_p[i];
      _array_p[i] = _array_f[i];
    }
  }

  RT& _p(const IT & i, const IT & j, const IT & k) { return _array_p[this->idx(i, j, k)]; }
  RT& _a(const IT & i, const IT & j, const IT & k) { return _array_a[this->idx(i, j, k)]; }
  RT& _c(const IT & i, const IT & j, const IT & k) { return _array_c[this->idx(i, j, k)]; }
  RT& _f(const IT & i, const IT & j, const IT & k) { return _array_f[this->idx(i, j, k)]; }

  RT& operator()(const IT & i, const IT & j, const IT & k)
  {
    return _array_a[this->idx(i, j, k)];
  }

  RT& operator[](IT idx)
  {
    return _array_a[idx];
  }

};

#endif
