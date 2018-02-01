
# Add your custom build info here.

LDFLAGS += -L${HOME}/local/lib

# Set the C++ compiler here. Otherwise, the environment variable $CXX is read.
# That's usually what you want.
#CXX=icpc

# Optimizations
OPTFLAGS += -O3
OPTFLAGS += -ffast-math
OPTFLAGS += -flto
OPTFLAGS += -march=native

# Warnings, etc.
CXXFLAGS += -Wall
CXXFLAGS += -fopenmp
CXXFLAGS += -g

CXXFLAGS += $(OPTFLAGS)

