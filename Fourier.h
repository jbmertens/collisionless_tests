#ifndef FOURIER
#define FOURIER

#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <iostream>
#include <omp.h>

#define ROUND_2_IT(f) ((IT)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

template<typename IT, typename RT>
class Fourier
{
  // FFT field
  fftwf_complex * f_field;
  // Array for FFTs
  RT * unghosted_field;

  IT nx, ny, nz; ///< # Samples in each dimension
  RT lx, ly, lz; ///< Physical domain lengths
  RT dx, dy, dz; ///< Length elements

  // plans for taking FFTs
  fftwf_plan p_c2r;
  fftwf_plan p_r2c;

  IT _F_IDX(IT i, IT j, IT k)
  {
    return i*ny*nz + j*nz + k;
  }

  IT _F_IDX_G(IT i, IT j, IT k, IT ng)
  {
    return (i+ng)*(ny+2*ng)*(nz+2*ng) + (j+ng)*(nz+2*ng) + k+ng;
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

    std::cout << "Initializing FFTW class with OpenMP support, using "
      << omp_get_max_threads() << " threads." << std::endl;
    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());

    //fftwf_malloc
    f_field = (fftwf_complex *) fftwf_malloc(nx*ny*(nz/2+1)
                                           *((long long) sizeof(fftwf_complex)));
    unghosted_field = new RT[nx*ny*nz];

    // create plans
    p_r2c = fftwf_plan_dft_r2c_3d(nx, ny, nz,
                                 tmp_field, f_field,
                                 FFTW_ESTIMATE);
    p_c2r = fftwf_plan_dft_c2r_3d(nx, ny, nz,
                                 f_field, tmp_field,
                                 FFTW_ESTIMATE);
  }

  ~Fourier()
  {
    // dealloc
    delete [] unghosted_field;
    fftwf_free(f_field);
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2r);
    fftwf_cleanup_threads();
  }

  void inverseLaplacian(RT *field, IT window_power, IT num_ghost_cells)
  {
    IT i, j, k;

    // copy non-ghost cells to unghosted_field
#pragma omp parallel for collapse(3)
    for(i=0; i<nx; ++i)
      for(j=0; j<ny; ++j)
        for(k=0; k<nz; ++k)
        {
          unghosted_field[_F_IDX(i,j,k)]
            = field[_F_IDX_G(i,j,k,num_ghost_cells)];
        }

    // inverse laplacian on unghosted_field
    inverseLaplacian(unghosted_field, window_power);

    // copy back to non-ghost cells
#pragma omp parallel for collapse(3)
    for(i=0; i<nx; ++i)
      for(j=0; j<ny; ++j)
        for(k=0; k<nz; ++k)
        {
          field[_F_IDX_G(i,j,k,num_ghost_cells)]
           = unghosted_field[_F_IDX(i,j,k)];
        }
  }

  void inverseLaplacian(RT *field)
  {
    inverseLaplacian(field, 0);
  }

  void inverseLaplacian(RT *field, IT window_power)
  {
    IT i, j, k;

    fftwf_execute_dft_r2c(p_r2c, field, f_field);

#pragma omp parallel for collapse(3)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz/2+1; k++)
        {
          RT px = (RT) (i<=nx/2 ? i : i-nx);
          RT py = (RT) (j<=ny/2 ? j : j-ny);
          RT pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          RT pmag = std::sqrt( (px/lx)*(px/lx)
            + (py/ly)*(py/ly) + (pz/lz)*(pz/lz) )*2.0*M_PI;

          RT window_corr = 1.0;
          if(window_power > 0)
          {
            RT sinc_x = px == 0 ? 1.0 : std::sin(M_PI*px/nx)/(M_PI*px/nx);
            RT sinc_y = py == 0 ? 1.0 : std::sin(M_PI*py/ny)/(M_PI*py/ny);
            RT sinc_z = pz == 0 ? 1.0 : std::sin(M_PI*pz/nz)/(M_PI*pz/nz);
            window_corr = std::pow(sinc_x*sinc_y*sinc_z, window_power);
          }

          f_field[fft_index][0] /= -pmag*pmag*nx*ny*nz / window_corr;
          f_field[fft_index][1] /= -pmag*pmag*nx*ny*nz / window_corr;
        }

    // zero mode?
    f_field[0][0] = 0;
    f_field[0][1] = 0;

    fftwf_execute_dft_c2r(p_c2r, f_field, field);
  }

  void periodicGradient(RT *field, IT grad_dir)
  {
    IT i, j, k;

    fftwf_execute_dft_r2c(p_r2c, field, f_field);

#pragma omp parallel for collapse(3)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz/2+1; k++)
        {
          RT px = (RT) (i<=nx/2 ? i : i-nx);
          RT py = (RT) (j<=ny/2 ? j : j-ny);
          RT pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          RT p = 2.0*M_PI * ( grad_dir==1 ? px/lx : (grad_dir==2 ? py/ly : pz/lz) );
          if(grad_dir == 1 && i == nx/2) p = 0;
          if(grad_dir == 2 && j == ny/2) p = 0;
          if(grad_dir == 3 && k == nz/2) p = 0;

          RT f_re = f_field[fft_index][0];
          RT f_im = f_field[fft_index][1];
          RT norm = 1.0/nx/ny/nz;

          f_field[fft_index][0] = -p*f_im*norm;
          f_field[fft_index][1] = p*f_re*norm;
        }

    // zero mode?
    f_field[0][0] = 0;
    f_field[0][1] = 0;

    fftwf_execute_dft_c2r(p_c2r, f_field, field);
  }

  void inverseGradient(RT *field, IT grad_dir)
  {
    IT i, j, k;

    fftwf_execute_dft_r2c(p_r2c, field, f_field);

#pragma omp parallel for collapse(3)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz/2+1; k++)
        {
          RT px = (RT) (i<=nx/2 ? i : i-nx);
          RT py = (RT) (j<=ny/2 ? j : j-ny);
          RT pz = (RT) k;

          IT fft_index = _FFT_IDX(i,j,k);

          RT p = 2.0*M_PI * ( grad_dir==1 ? px/lx : (grad_dir==2 ? py/ly : pz/lz) );
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

    fftwf_execute_dft_c2r(p_c2r, f_field, field);
  }

  /**
   * @brief Gaussian random field ICs
   */
  void gaussianRandomRealization(RT * field, IT num_ghost_cells)
  {
    // GRF in internal array
    gaussianRandomRealization(unghosted_field);
    // copy to field w/ ghosts
    IT i, j, k;
#pragma omp parallel for collapse(3)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz; k++)
        {
          field[_F_IDX_G(i,j,k,num_ghost_cells)]
            = unghosted_field[_F_IDX(i,j,k)];
        }
  }

  /**
   * @brief Gaussian random field ICs
   */
  void gaussianRandomRealization(RT * field)
  {
    IT i, j, k;
    RT scale;

    // don't expect to run a simulation larger than this
    // Otherwise, ICs will be inconsistent between resolutions
    IT MAX_N = 4096;

#pragma omp parallel for collapse(3)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz/2+1; k++)
        {
if(std::abs(i) < 4 && std::abs(j) < 4 && std::abs(k) < 4)
{ // PS cutoff
          RT px = (RT) (i<=nx/2 ? i : i-nx);
          RT py = (RT) (j<=ny/2 ? j : j-ny);
          RT pz = (RT) k;

          RT pmag = std::sqrt( (px/lx)*(px/lx) + (py/ly)*(py/ly)
            + (pz/lz)*(pz/lz) )*2.0*M_PI;

          IT fft_index = _FFT_IDX(i, j, k);

          IT seed = MAX_N*MAX_N*i + MAX_N*j + k;

          // initialize rng using index as a seed
          std::random_device rd;
          std::mt19937 gen(seed);
          std::normal_distribution<RT> gaussian_distribution;
          std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*M_PI);

          RT rand_mag = gaussian_distribution(gen);
          RT rand_phase = angular_distribution(gen);

          // Scale by power spectrum
          scale = 1.0*std::pow(pmag, -3);
          f_field[fft_index][0] = scale*rand_mag*std::cos(rand_phase);
          f_field[fft_index][1] = scale*rand_mag*std::sin(rand_phase);
}
        }

    // No zero-mode
    f_field[_FFT_IDX(0,0,0)][0] = 0;
    f_field[_FFT_IDX(0,0,0)][1] = 0;

    // FFT back; 'field' array should now be populated
    // with a gaussian random field
    fftwf_execute_dft_c2r(p_c2r, f_field, field);
  }

};

#endif
