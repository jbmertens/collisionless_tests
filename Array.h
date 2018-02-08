#ifndef ARRAY_H
#define ARRAY_H

#include <limits>
#include <stdexcept>
#include "TriCubicInterpolator.h"

#define ARRAY_LOOP(i, pts) \
  for(IT i=0; i<pts; ++i)

#define PARALLEL_ARRAY_LOOP(i, pts) \
  _Pragma("omp parallel for") \
  ARRAY_LOOP(i, pts)

// number of ghost cells
#define NG 2

enum arrayInterpolationType { Trilinear, LMTricubic, CINTTricubic };

template<typename IT, typename RT>
class Array
{
protected:

  IT _nx, _ny, _nz; ///< array dims
  IT _nx_inc, _ny_inc, _nz_inc; ///< array dims incl. ghost zones
  RT _lx, _ly, _lz; ///< physical dims
  RT _dx, _dy, _dz;
  IT _pts_inc; ///< # points (incl. ghost)
  IT _pts_exc; ///< # points (excl. ghost)

  RT * _array;
  
  arrayInterpolationType interpolation_type;

private:

  void _init(IT nx_in, IT ny_in, IT nz_in,
    RT lx_in, RT ly_in, RT lz_in)
  {
    _nx = nx_in;
    _ny = ny_in;
    _nz = nz_in;

    _nx_inc = nx_in + 2*NG;
    _ny_inc = ny_in + 2*NG;
    _nz_inc = nz_in + 2*NG;

    _lx = lx_in;
    _ly = ly_in;
    _lz = lz_in;

    _dx = lx_in / (RT) nx_in;
    _dy = ly_in / (RT) ny_in;
    _dz = lz_in / (RT) nz_in;

    _pts_inc = _nx_inc*_ny_inc*_nz_inc;
    _pts_exc = _nx*_ny*_nz;

    _array = new RT[_pts_inc];

    PARALLEL_ARRAY_LOOP(i, _pts_inc)
    {
      _array[i] = 0.0;
    }

    interpolation_type = CINTTricubic;
  }

public:

  const IT &nx, &ny, &nz;
  const RT &lx, &ly, &lz;
  const RT &dx, &dy, &dz;

  Array(IT n_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(n_in, n_in, n_in, 1.0, 1.0, 1.0);
  }

  Array(IT nx_in, IT ny_in, IT nz_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(nx_in, ny_in, nz_in, 1.0, 1.0, 1.0);
  }

  Array(IT nx_in, IT ny_in, IT nz_in,
    RT lx_in, RT ly_in, RT lz_in) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    _init(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in);
  }

  Array(IT nx_in, IT ny_in, IT nz_in,
    RT lx_in, RT ly_in, RT lz_in, arrayInterpolationType type) :
    nx(_nx), ny(_ny), nz(_nz),
    lx(_lx), ly(_ly), lz(_lz),
    dx(_dx), dy(_dy), dz(_dz)
  {
    interpolation_type = type;
    _init(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in);
  }

  ~Array()
  {
    delete [] _array;
  }

  void setInterpolationType(arrayInterpolationType type)
  {
    interpolation_type = type;
  }

  IT idx(const IT & i_in, const IT & j_in, const IT & k_in) const
  {
    return (i_in+NG)*_ny_inc*_nz_inc + (j_in+NG)*_nz_inc + (k_in+NG);
  }

  RT& operator[](IT idx)
  {
    return _array[idx];
  }

  RT& operator()(const IT & i, const IT & j, const IT & k)
  {
    return _array[idx(i, j, k)];
  }

  Array& operator=(const Array& other) {
    // check for self-assignment
    if(&other == this)
      return *this;

    if(this->_nx != other.nx || this->_ny != other.ny || this->_nz != other.nz)
      throw std::runtime_error ("Incompatible array shapes.");

    PARALLEL_ARRAY_LOOP(i, this->_pts_inc)
    {
      this->_array[i] = other._array[i];
    }

    return *this;
  }


  void addPeriodicGhostContributions()
  {
    // contributions from 6 faces

    // x faces
#pragma omp parallel for
    for(IT j=-NG; j<_ny+NG; ++j)
      for(IT k=-NG; k<_nz+NG; ++k)
      {
        for(IT i=-NG; i<0; ++i)
        {
          _array[idx(i+_nx, j, k)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
        for(IT i=_nx; i<_nx+NG; ++i)
        {
          _array[idx(i-_nx, j, k)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
      }

    // y faces
#pragma omp parallel for
    for(IT i=-NG; i<_nx+NG; ++i)
      for(IT k=-NG; k<_nz+NG; ++k)
      {
        for(IT j=-NG; j<0; ++j)
        {
          _array[idx(i, j+_ny, k)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
        for(IT j=_ny; j<_ny+NG; ++j)
        {
          _array[idx(i, j-_ny, k)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
      }

    // z faces
#pragma omp parallel for
    for(IT i=-NG; i<_nx+NG; ++i)
      for(IT j=-NG; j<_ny+NG; ++j)
      {
        for(IT k=-NG; k<0; ++k)
        {
          _array[idx(i, j, k+_nz)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
        for(IT k=_nz; k<_nz+NG; ++k)
        {
          _array[idx(i, j, k-_nz)] += _array[idx(i, j, k)];
          _array[idx(i, j, k)] = 0; // contribution already added; zero ghost cell value.
        }
      }
  }

  void applyPeriodicBoundaryConditions()
  {
    // contributions from 6 faces

    // x faces
    for(IT j=0; j<_ny; ++j)
      for(IT k=0; k<_nz; ++k)
      {
        for(IT i=-NG; i<0; ++i)
          _array[idx(i, j, k)] = _array[idx(i+_nx, j, k)];
        for(IT i=_nx; i<_nx+NG; ++i)
          _array[idx(i, j, k)] = _array[idx(i-_nx, j, k)];
      }

    // y faces (plus x ghosts)
    for(IT i=-NG; i<_nx+NG; ++i)
      for(IT k=0; k<_nz; ++k)
      {
        for(IT j=-NG; j<0; ++j)
          _array[idx(i, j, k)] = _array[idx(i, j+_ny, k)];
        for(IT j=_ny; j<_ny+NG; ++j)
          _array[idx(i, j, k)] = _array[idx(i, j-_ny, k)];
      }

    // z faces (plus x and y ghosts)
    for(IT i=-NG; i<_nx+NG; ++i)
      for(IT j=-NG; j<_ny+NG; ++j)
      {
        for(IT k=-NG; k<0; ++k)
          _array[idx(i, j, k)] = _array[idx(i, j, k+_nz)];
        for(IT k=_nz; k<_nz+NG; ++k)
          _array[idx(i, j, k)] = _array[idx(i, j, k-_nz)];
      }
  }

  RT getInterpolatedValueAtModX(RT i_in, RT j_in, RT k_in) const
  {
    RT i_mod = i_in - _nx*std::floor(i_in/_nx);
    RT j_mod = j_in - _ny*std::floor(j_in/_ny);
    RT k_mod = k_in - _nz*std::floor(k_in/_nz);
    return getInterpolatedValue(i_mod, j_mod, k_mod);
  }

  RT getInterpolatedValue(RT i_in, RT j_in, RT k_in) const
  {
    switch(interpolation_type)
    {
      case Trilinear:
      default:
        return getTriLinearInterpolatedValue(i_in, j_in, k_in);
        break;

      case LMTricubic:
        return getLMTriCubicInterpolatedValue(i_in, j_in, k_in);
        break;

      case CINTTricubic:
        return getCINTTriCubicInterpolatedValue(i_in, j_in, k_in);
        break;
    }
  }

  // Weighted averaging / trilinear interpolation via
  // https://en.wikipedia.org/wiki/Trilinear_interpolation#Method
  RT getTriLinearInterpolatedValue(RT i_in, RT j_in, RT k_in) const
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

  // Lekien-Marsden cubic spline
  RT getLMTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in) const
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
  RT CINT(RT u, RT p0, RT p1, RT p2, RT p3) const
  {
    return 0.5*(
          (u*u*(2.0 - u) - u)*p0
        + (u*u*(3.0*u - 5.0) + 2)*p1
        + (u*u*(4.0 - 3.0*u) + u)*p2
        + u*u*(u - 1.0)*p3
      );
  }

  RT getCINTTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in) const
  {
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

  RT xDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
      + 1.0/12.0*_array[idx(i_in - 2, j_in, k_in)]
      - 2.0/3.0*_array[idx(i_in - 1, j_in, k_in)]
      + 2.0/3.0*_array[idx(i_in + 1, j_in, k_in)]
      - 1.0/12.0*_array[idx(i_in + 2, j_in, k_in)]
    ) / _dx;
  }

  RT yDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
      + 1.0/12.0*_array[idx(i_in, j_in - 2, k_in)]
      - 2.0/3.0*_array[idx(i_in, j_in - 1, k_in)]
      + 2.0/3.0*_array[idx(i_in, j_in + 1, k_in)]
      - 1.0/12.0*_array[idx(i_in, j_in + 2, k_in)]
    ) / _dy;
  }

  RT zDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
     + 1.0/12.0*_array[idx(i_in, j_in, k_in - 2)]
     - 2.0/3.0*_array[idx(i_in, j_in, k_in - 1)]
     + 2.0/3.0*_array[idx(i_in, j_in, k_in + 1)]
     - 1.0/12.0*_array[idx(i_in, j_in, k_in + 2)]
    ) / _dz;
  }

  RT xxDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
      - 1.0/12.0*_array[idx(i_in + 2, j_in, k_in)]
      + 4.0/3.0*_array[idx(i_in + 1, j_in, k_in)]
      - 5.0/2.0*_array[idx(i_in + 0, j_in, k_in)]
      + 4.0/3.0*_array[idx(i_in - 1, j_in, k_in)]
      - 1.0/12.0*_array[idx(i_in - 2, j_in, k_in)]
    ) / _dx / _dx;
  }

  RT yyDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
      - 1.0/12.0*_array[idx(i_in, j_in + 2, k_in)]
      + 4.0/3.0*_array[idx(i_in, j_in + 1, k_in)]
      - 5.0/2.0*_array[idx(i_in, j_in + 0, k_in)]
      + 4.0/3.0*_array[idx(i_in, j_in - 1, k_in)]
      - 1.0/12.0*_array[idx(i_in, j_in - 2, k_in)]
    ) / _dy / _dy;
  }

  RT zzDer(RT i_in, RT j_in, RT k_in) const
  {
    return (
      - 1.0/12.0*_array[idx(i_in, j_in, k_in + 2)]
      + 4.0/3.0*_array[idx(i_in, j_in, k_in + 1)]
      - 5.0/2.0*_array[idx(i_in, j_in, k_in + 0)]
      + 4.0/3.0*_array[idx(i_in, j_in, k_in - 1)]
      - 1.0/12.0*_array[idx(i_in, j_in, k_in - 2)]
    ) / _dz / _dz;
  }

  RT xyDer(RT i_in, RT j_in, RT k_in) const
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

  RT xzDer(RT i_in, RT j_in, RT k_in) const
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

  RT yzDer(RT i_in, RT j_in, RT k_in) const
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
