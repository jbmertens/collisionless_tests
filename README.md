# Newtonian sheet method.

Code here is a not-quite-complete implementation of the sheet method.

For now, compile by doing something like:
```
g++ main.cc -std=c++11 -lfftw3 -O3 -ffast-math -flto -Wall -lz -march=native
```

TODO:

 - Figure out units
 - Work on lightcone / raytracing aspect
 - Random power spectrum realization
   - Using 2LPT for ICs? https://arxiv.org/pdf/0910.0258.pdf
   - Dx_i = D_1 (d_i phi) + ...
   - v_i = -D_1 f_1 H (d_i phi)
   - ... For some phi.
 - Window function correction, TSC, etc?
 - Parallelize
