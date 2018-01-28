# Newtonian sheet method.

Code here is a not-quite-complete implementation of the sheet method.

For now, compile by doing something like:
```
g++ main.cc -std=c++11 -lfftw3_omp -lfftw3 -lm -flto -O3 -Wall -lz -ffast-math -march=native -fopenmp
```

TODO:

 - Figure out units
 - Work on lightcone / raytracing aspect
 - Window function correction, add'l deposition schemes, etc?
 - Parallelize
