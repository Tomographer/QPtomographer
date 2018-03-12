
#ifndef DNORMTOMO_USE_SCS
#define DNORMTOMO_USE_SCS

//
// NOTES: INSTRUCTIONS FOR COMPILING SCS: EDIT YOUR scs-X.Y.Z/scs.mk, AND FIND AND CHANGE
// THE FOLLOWING VARIABLES TO READ:
//
// ...
// CTRLC = 0
// ...
// USE_OPENMP = 0
// ...
// USE_LAPACK = 1
// ...
//


namespace SCS {
  extern "C" {
#   include <scs.h>
#   include <linSys.h>
#   include <linsys/amatrix.h>
  }
} // namespace SCS





namespace tomo_internal {
inline SCS::scs_int ltrilinindex(SCS::scs_int i, SCS::scs_int j, SCS::scs_int d)
{
  assert(i >= j && i < d);
  // Drafts & Calculations Vol VI ~60% 11/29/2016
  //return j*(d-j+1) + j*(j-1)/2 + i - j;
  return d*j - (j*j+j)/2 + i;
}
}




#endif
