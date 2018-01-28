#ifndef PERIODICARRAY_H
#define PERIODICARRAY_H

#include <limits>
#include <stdexcept>
#include "TriCubicInterpolator.h"

template<typename IT, typename RT>
class PeriodicArray
{
protected:

  IT _nx, _ny, _nz;
  RT _lx, _ly, _lz;
  RT _dx, _dy, _dz;
  IT _pts;

  // TODO: __restrict__
  RT * _array;

private:

  void _init(IT nx_in, IT ny_in, IT nz_in,
    RT lx_in, RT ly_in, RT lz_in)
  {
    _nx = nx_in;
    _ny = ny_in;
    _nz = nz_in;

    _lx = lx_in;
    _ly = ly_in;
    _lz = lz_in;

    _dx = lx_in / (RT) nx_in;
    _dy = ly_in / (RT) ny_in;
    _dz = lz_in / (RT) nz_in;

    _pts = _nx*_ny*_nz;

    _array = new RT[_pts];
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<_pts; ++i)
    {
      _array[i] = 0.0;
    }
  }

  IT _IT_mod(IT n, IT d) const
  {
    // TODO: check and see if having conditional here is faster
    if ( n >= 0 && n < d )
      return n;

    IT mod = n % d;
    if(mod < 0)
      mod += d;
    return mod;
  }

public:

  const IT &nx, &ny, &nz;
  const RT &lx, &ly, &lz;
  const RT &dx, &dy, &dz;

  PeriodicArray(IT n_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(n_in, n_in, n_in, 1.0, 1.0, 1.0);
  }

  PeriodicArray(IT nx_in, IT ny_in, IT nz_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(nx_in, ny_in, nz_in, 1.0, 1.0, 1.0);
  }

  PeriodicArray(IT nx_in, IT ny_in, IT nz_in, RT lx_in, RT ly_in, RT lz_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in);
  }

  ~PeriodicArray()
  {
    delete [] _array;
  }

  IT idx(IT i_in, IT j_in, IT k_in) const // TODO
  {
    IT i=_IT_mod(i_in, _nx), j=_IT_mod(j_in, _ny), k=_IT_mod(k_in, _nz);
    return i*_ny*_nz + j*_nz + k;
  }

  RT& operator()(const IT & i, const IT & j, const IT & k)
  {
    IT p = idx(i, j, k);
    return _array[p];
  }

  PeriodicArray& operator=(const PeriodicArray& other) {
    // check for self-assignment
    if(&other == this)
      return *this;

    RT this_pts = this->_nx*this->_ny*this->_nz;
    RT other_pts = other.nx*other.ny*other.nz;
    if(this_pts != other_pts)
      throw std::runtime_error ("Incompatible array shapes.");

    this->_nx = other.nx;
    this->_nz = other.ny;
    this->_ny = other.nz;
    this->_lx = other.lx;
    this->_lz = other.ly;
    this->_ly = other.lz;
    this->_dx = other.dx;
    this->_dz = other.dy;
    this->_dy = other.dz;
    this->_pts = this->_nx*this->_ny*this->_nz;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i=0; i<this->_pts; ++i)
    {
      this->_array[i] = other._array[i];
    }

    return *this;
  }

  RT& operator[](IT idx)
  {
    return _array[idx];
  }

  RT sum()
  {
    RT res = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:res)
#endif
    for(IT i=0; i<_pts; ++i)
    {
      res += _array[i];
    }
    return res;

  }
  
  RT avg()
  {
    return sum() / (RT) _pts;
  }

  RT min()
  {
    RT min_res = std::numeric_limits<RT>::max();
#ifdef USE_OPENMP
#pragma omp parallel for 
#endif
    for(IT i = 0; i<_pts; i++)
    {
#ifdef USE_OPENMP
#pragma omp critical
#endif
      if(_array[i] < min_res)
        min_res = _array[i];
    }
    return min_res;
  }

  RT max()
  {
    RT max_res = std::numeric_limits<RT>::min();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(IT i = 0; i<_pts; i++)
    {
#ifdef USE_OPENMP
#pragma omp critical
#endif
      {
        if(_array[i] > max_res)
          max_res = _array[i];
      }
    }

    return max_res;
  }

  // Weighted averaging / trilinear interpolation via
  // https://en.wikipedia.org/wiki/Trilinear_interpolation#Method
  RT getInterpolatedValue(RT i_in, RT j_in, RT k_in)
  {
    IT il = i_in < 0 ? (IT) i_in - 1 : (IT) i_in; // Index "left" of i
    RT id = i_in - il; // fractional difference
    IT jl = j_in < 0 ? (IT) j_in - 1 : (IT) j_in; // same as ^ but j
    RT jd = j_in - jl;
    IT kl = k_in < 0 ? (IT) k_in - 1 : (IT) k_in; // same as ^ but k
    RT kd = k_in - kl;

    RT c00 = _array[idx(il, jl, kl)]*(1-id) + _array[idx(il+1, jl, kl)]*id;
    RT c01 = _array[idx(il, jl, kl+1)]*(1-id) + _array[idx(il+1, jl, kl+1)]*id;
    RT c10 = _array[idx(il, jl+1, kl)]*(1-id) + _array[idx(il+1, jl+1, kl)]*id;
    RT c11 = _array[idx(il, jl+1, kl+1)]*(1-id) + _array[idx(il+1, jl+1, kl+1)]*id;
    RT c0 = c00*(1-jd) + c10*jd;
    RT c1 = c01*(1-jd) + c11*jd;

    return c0*(1-kd) + c1*kd;
  }

  RT getTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in)
  {
    // TODO: consider higher order interpolation methods
    return getCINTTriCubicInterpolatedValue(i_in, j_in, k_in);
    return getLMTriCubicInterpolatedValue(i_in, j_in, k_in);
  }

  // Lekien-Marsden cubic spline
  RT getLMTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in)
  {
    IT il = i_in < 0 ? (IT) i_in - 1 : (IT) i_in; // Index "left" of i
    RT id = i_in - il; // fractional difference
    IT jl = j_in < 0 ? (IT) j_in - 1 : (IT) j_in; // same as ^ but j
    RT jd = j_in - jl;
    IT kl = k_in < 0 ? (IT) k_in - 1 : (IT) k_in; // same as ^ but k
    RT kd = k_in - kl;

    RT a[64];
    RT f[64];

    for(IT i=0; i<4; ++i)
      for(IT j=0; j<4; ++j)
        for(IT k=0; k<4; ++k)
        {
          f[i*16 + j*4 + k] = _array[idx(il+i-1, jl+j-1, kl-1)];
        }

    compute_tricubic_coeffs(a, f);
    return evaluate_interpolation(a, id, jd, kd);
  }

  // Catmull-Rom cubic spline
  RT CINT(RT u, RT p0, RT p1, RT p2, RT p3)
  {
    return 0.5*(
          (u*u*(2.0 - u) - u)*p0
        + (u*u*(3.0*u - 5.0) + 2)*p1
        + (u*u*(4.0 - 3.0*u) + u)*p2
        + u*u*(u - 1.0)*p3
      );
  }

  RT getCINTTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in)
  {
    // TODO: consider higher order interpolation methods
    IT il = i_in < 0 ? (IT) i_in - 1 : (IT) i_in; // Index "left" of i
    RT id = i_in - il; // fractional difference
    IT jl = j_in < 0 ? (IT) j_in - 1 : (IT) j_in; // same as ^ but j
    RT jd = j_in - jl;
    IT kl = k_in < 0 ? (IT) k_in - 1 : (IT) k_in; // same as ^ but k
    RT kd = k_in - kl;

    // interpolated value at (i*, j*, k_in)
    RT * F_i_j_kd = new RT[16];
    for(IT i=0; i<4; ++i)
      for(IT j=0; j<4; ++j)
        F_i_j_kd[i*4+j] = CINT(kd,
          _array[idx(il+i-1, jl+j-1, kl-1)], _array[idx(il+i-1, jl+j-1, kl+0)],
          _array[idx(il+i-1, jl+j-1, kl+1)], _array[idx(il+i-1, jl+j-1, kl+2)]);

    // interpolated value at (i*, j_in, k_in)
    RT * F_i_jd_kd = new RT[4];
    for(IT i=0; i<4; ++i)
      F_i_jd_kd[i] = CINT(jd, F_i_j_kd[i*4+0], F_i_j_kd[i*4+1], F_i_j_kd[i*4+2], F_i_j_kd[i*4+3]);

    // interpolated value at (i_in, j_in, k_in)
    RT Fijk = CINT(id, F_i_jd_kd[0], F_i_jd_kd[1], F_i_jd_kd[2], F_i_jd_kd[3]);

    delete [] F_i_j_kd;
    delete [] F_i_jd_kd;

    return Fijk;
  }

  RT xDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      + 1.0/12.0*_array[idx(i_in - 2, j_in, k_in)]
      - 2.0/3.0*_array[idx(i_in - 1, j_in, k_in)]
      + 2.0/3.0*_array[idx(i_in + 1, j_in, k_in)]
      - 1.0/12.0*_array[idx(i_in + 2, j_in, k_in)]
    ) / _dx;
  }

  RT yDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      + 1.0/12.0*_array[idx(i_in, j_in - 2, k_in)]
      - 2.0/3.0*_array[idx(i_in, j_in - 1, k_in)]
      + 2.0/3.0*_array[idx(i_in, j_in + 1, k_in)]
      - 1.0/12.0*_array[idx(i_in, j_in + 2, k_in)]
    ) / _dy;
  }

  RT zDer(RT i_in, RT j_in, RT k_in)
  {
    return (
     + 1.0/12.0*_array[idx(i_in, j_in, k_in - 2)]
     - 2.0/3.0*_array[idx(i_in, j_in, k_in - 1)]
     + 2.0/3.0*_array[idx(i_in, j_in, k_in + 1)]
     - 1.0/12.0*_array[idx(i_in, j_in, k_in + 2)]
    ) / _dz;
  }

  RT xxDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      - 1.0/12.0*_array[idx(i_in + 2, j_in, k_in)]
      + 4.0/3.0*_array[idx(i_in + 1, j_in, k_in)]
      - 5.0/2.0*_array[idx(i_in + 0, j_in, k_in)]
      + 4.0/3.0*_array[idx(i_in - 1, j_in, k_in)]
      - 1.0/12.0*_array[idx(i_in - 2, j_in, k_in)]
    ) / _dx / _dx;
  }

  RT yyDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      - 1.0/12.0*_array[idx(i_in, j_in + 2, k_in)]
      + 4.0/3.0*_array[idx(i_in, j_in + 1, k_in)]
      - 5.0/2.0*_array[idx(i_in, j_in + 0, k_in)]
      + 4.0/3.0*_array[idx(i_in, j_in - 1, k_in)]
      - 1.0/12.0*_array[idx(i_in, j_in - 2, k_in)]
    ) / _dy / _dy;
  }

  RT zzDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      - 1.0/12.0*_array[idx(i_in, j_in, k_in + 2)]
      + 4.0/3.0*_array[idx(i_in, j_in, k_in + 1)]
      - 5.0/2.0*_array[idx(i_in, j_in, k_in + 0)]
      + 4.0/3.0*_array[idx(i_in, j_in, k_in - 1)]
      - 1.0/12.0*_array[idx(i_in, j_in, k_in - 2)]
    ) / _dz / _dz;
  }

  RT xyDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      _array[idx(i_in + 1, j_in + 1, k_in)]
      + _array[idx(i_in - 1, j_in - 1, k_in)]
      - _array[idx(i_in - 1, j_in + 1, k_in)]
      - _array[idx(i_in + 1, j_in - 1, k_in)]
      - 1.0/16.0*(
        _array[idx(i_in + 2, j_in + 2, k_in)]
        + _array[idx(i_in - 2, j_in - 2, k_in)]
        - _array[idx(i_in - 2, j_in + 2, k_in)]
        - _array[idx(i_in + 2, j_in - 2, k_in)]
      )
    ) / 3.0 / _dx / _dy;
  }

  RT xzDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      _array[idx(i_in + 1, j_in, k_in + 1)]
      + _array[idx(i_in - 1, j_in, k_in - 1)]
      - _array[idx(i_in - 1, j_in, k_in + 1)]
      - _array[idx(i_in + 1, j_in, k_in - 1)]
      - 1.0/16.0*(
        _array[idx(i_in + 2, j_in, k_in + 2)]
        + _array[idx(i_in - 2, j_in, k_in - 2)]
        - _array[idx(i_in - 2, j_in, k_in + 2)]
        - _array[idx(i_in + 2, j_in, k_in - 2)]
      )
    ) / 3.0 / _dx / _dz;
  }

  RT yzDer(RT i_in, RT j_in, RT k_in)
  {
    return (
      _array[idx(i_in, j_in + 1, k_in + 1)]
      + _array[idx(i_in, j_in - 1, k_in - 1)]
      - _array[idx(i_in, j_in - 1, k_in + 1)]
      - _array[idx(i_in, j_in + 1, k_in - 1)]
      - 1.0/16.0*(
        _array[idx(i_in, j_in + 2, k_in + 2)]
        + _array[idx(i_in, j_in - 2, k_in - 2)]
        - _array[idx(i_in, j_in - 2, k_in + 2)]
        - _array[idx(i_in, j_in + 2, k_in - 2)]
      )
    ) / 3.0 / _dy / _dz;
  }

};



#endif
