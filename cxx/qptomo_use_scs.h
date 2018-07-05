
#ifndef QPTOMO_USE_SCS
#define QPTOMO_USE_SCS

//
// NOTE: INSTRUCTIONS FOR COMPILING SCS CAN BE FOUND IN OUR MAIN "README.TXT"
//


//namespace SCS {
  extern "C" {
#   include <scs.h>
#   include <linsys.h>
#   include <linsys/amatrix.h>
  }
//} // namespace SCS




namespace tomo_internal {
inline scs_int ltrilinindex(scs_int i, scs_int j, scs_int d)
{
  assert(i >= j && i < d);
  // Drafts & Calculations Vol VI ~60% 11/29/2016
  //return j*(d-j+1) + j*(j-1)/2 + i - j;
  return d*j - (j*j+j)/2 + i;
}
}




#endif
