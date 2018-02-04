#ifndef COSMOLOGY
#define COSMOLOGY

#include <cmath>

/**
 * @file      Standard cosmology functions,
 * assuming a matter-dominated cosmology (omega_m = 1).
 * TODO: generalize.
 */

namespace cosmology
{

/**
 * @brief      Return scale factor at time t, using units
 *  where a = 1 at t = 1.
 */
template<typename RT>
RT a(RT t)
{
  // want a = 1 at t = 1 => a0 = 1
  return std::pow(t, 2.0/3.0);
}

/**
 * @brief      Hubble parameter.
 */
template<typename RT>
RT H(RT t)
{
  return 2.0/3.0/t;
}

/**
 * @brief      Growth function
 */
template<typename RT>
RT D(RT t)
{
  return a<RT>(t);
}

/**
 * @brief      Growth rate
 */
template<typename RT>
RT f(RT t)
{
  return 1.0;
}

/**
 * @brief      2nd-order growth function
 */
template<typename RT>
RT D2(RT t)
{
  return -3.0*std::pow(D<RT>(t),2)/7.0;
}

/**
 * @brief      2nd-order growth rate
 */
template<typename RT>
RT f2(RT t)
{
  return 2.0;
}

}; // namespace

#endif