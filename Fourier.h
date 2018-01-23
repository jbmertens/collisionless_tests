#ifndef FOURIER
#define FOURIER

#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <cmath>

template<typename IT, typename RT>
class Fourier
{
  // FFT field
  fftw_complex *f_field;

  IT nx, ny, nz; ///< # Samples in each dimension
  RT lx, ly, lz; ///< Physical domain lengths
  RT dx, dy, dz; ///< Length elements

  // plans for taking FFTs
  fftw_plan p_c2r;
  fftw_plan p_r2c;

  IT _F_IDX(IT i, IT j, IT k)
  {
    return i*ny*nz + j*nz + k;
  }

  IT _FFT_IDX(IT i, IT j, IT k)
  {
    return i*ny*(nz/2+1) + j*(nz/2+1) + k;
  }

public:

  /**
   * @brief Initialize a fourier class instance
   * @details Create fftw plans, allocate memory
   * 
   * @param nx points in x-direction
   * @param ny points in y-direction
   * @param nz points in z-direction
   * @param field any nx*ny*nz grid for planning
   */
  Fourier(IT nx_in, IT ny_in, IT nz_in, RT lx_in, RT ly_in, RT lz_in,
    RT *tmp_field)
  {
    nx = nx_in; ny = ny_in; nz = nz_in;
    lx = lx_in; ly = ly_in; lz = lz_in;
    dx = lx/nx; dy = ly/ny; dz = lz/nz;

    //fftw_malloc
    f_field = (fftw_complex *) fftw_malloc(nx*ny*(nz/2+1)
                                           *((long long) sizeof(fftw_complex)));

    // create plans
    p_r2c = fftw_plan_dft_r2c_3d(nx, ny, nz,
                                 tmp_field, f_field,
                                 FFTW_MEASURE);
    p_c2r = fftw_plan_dft_c2r_3d(nx, ny, nz,
                                 f_field, tmp_field,
                                 FFTW_MEASURE);
  }

  ~Fourier()
  {
    // dealloc
    fftw_free(f_field);
    fftw_destroy_plan(p_r2c);
    fftw_destroy_plan(p_c2r);
  }

  void inverseLaplacian(RT *field)
  {
    IT i, j, k;
    RT px, py, pz, pmag;

    fftw_execute_dft_r2c(p_r2c, field, f_field);

    for(i=0; i<nx; i++)
    {
      px = (RT) (i<=nx/2 ? i : i-nx);
      for(j=0; j<ny; j++)
      {
        py = (RT) (j<=ny/2 ? j : j-ny);
        for(k=0; k<nz/2+1; k++)
        {
          pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          pmag = std::sqrt( (px/lx)*(px/lx) + (py/ly)*(py/ly) + (pz/lz)*(pz/lz) )*2.0*M_PI;

          f_field[fft_index][0] /= -pmag*pmag*nx*ny*nz;
          f_field[fft_index][1] /= -pmag*pmag*nx*ny*nz;
        }
      }
    }
    // zero mode?
    f_field[0][0] = 0;
    f_field[0][1] = 0;

    fftw_execute_dft_c2r(p_c2r, f_field, field);
  }

  /**
   * Compute derivative of an "offset periodic" function,
   * Ie. a function f(s) that is periodic in f(s) - s.
   * Then, f'(s) = d_s ( f(s) - s ) + 0 or 1 if s = s
   * And d_s ( f(s) - s ) can be computed via FFT.
   */
  void offsetPeriodicGradient(RT *field, IT offset_dir, IT grad_dir)
  {
    RT is_i_offset = offset_dir == 1 ? 1 : 0;
    RT is_j_offset = offset_dir == 2 ? 1 : 0;
    RT is_k_offset = offset_dir == 3 ? 1 : 0;

    // subtract off "s" gradient
    for(IT i=0; i<nx; ++i)
      for(IT j=0; j<ny; ++j)
        for(IT k=0; k<nz; ++k)
        {
          field[_F_IDX(i,j,k)] -= i*dx*is_i_offset
            + j*dy*is_j_offset + k*dz*is_k_offset;
        }

    // compute "normal" gradient
    periodicGradient(field, grad_dir);

    // add 1 back in if grad_dir is offset_dir
    if(offset_dir == grad_dir)
      for(IT i=0; i<nx; ++i)
        for(IT j=0; j<ny; ++j)
          for(IT k=0; k<nz; ++k)
          {
            field[_F_IDX(i,j,k)] += 1;
          }

  }

  void periodicGradient(RT *field, IT grad_dir)
  {
    IT i, j, k;
    RT px, py, pz, p;

    fftw_execute_dft_r2c(p_r2c, field, f_field);

    for(i=0; i<nx; i++)
    {
      px = (RT) (i<=nx/2 ? i : i-nx);
      for(j=0; j<ny; j++)
      {
        py = (RT) (j<=ny/2 ? j : j-ny);
        for(k=0; k<nz/2+1; k++)
        {
          pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          p = 2.0*M_PI * ( grad_dir==1 ? px/lx : (grad_dir==2 ? py/ly : pz/lz) );
          if(grad_dir == 1 && i == nx/2) p = 0;
          if(grad_dir == 2 && j == ny/2) p = 0;
          if(grad_dir == 3 && k == nz/2) p = 0;

          RT f_re = f_field[fft_index][0];
          RT f_im = f_field[fft_index][1];
          RT norm = 1.0/nx/ny/nz;

          f_field[fft_index][0] = -p*f_im*norm;
          f_field[fft_index][1] = p*f_re*norm;
        }
      }
    }
    // zero mode?
    f_field[0][0] = 0;
    f_field[0][1] = 0;

    fftw_execute_dft_c2r(p_c2r, f_field, field);
  }

  void inverseGradient(RT *field, IT grad_dir)
  {
    IT i, j, k;
    RT px, py, pz, p;

    fftw_execute_dft_r2c(p_r2c, field, f_field);

    for(i=0; i<nx; i++)
    {
      px = (RT) (i<=nx/2 ? i : i-nx);
      for(j=0; j<ny; j++)
      {
        py = (RT) (j<=ny/2 ? j : j-ny);
        for(k=0; k<nz/2+1; k++)
        {
          pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          p = 2.0*M_PI * ( grad_dir==1 ? px/lx : (grad_dir==2 ? py/ly : pz/lz) );
          if(grad_dir == 1 && i == nx/2) p = 0;
          if(grad_dir == 2 && j == ny/2) p = 0;
          if(grad_dir == 3 && k == nz/2) p = 0;

          RT f_re = f_field[fft_index][0];
          RT f_im = f_field[fft_index][1];
          RT norm = 1.0/nx/ny/nz;
          
          if(p == 0)
          {
            f_field[fft_index][0] = 0;
            f_field[fft_index][1] = 0;
          }
          else
          {
            f_field[fft_index][0] = f_im/p*norm;
            f_field[fft_index][1] = -f_re/p*norm;
          }
        }
      }
    }

    fftw_execute_dft_c2r(p_c2r, f_field, field);
  }

};

#endif