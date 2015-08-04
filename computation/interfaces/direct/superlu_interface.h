#ifndef SELDON_FILE_SUPERLU_INTERFACE_H /* allow multiple inclusions */
#define SELDON_FILE_SUPERLU_INTERFACE_H

#ifdef SELDON_WITH_SUPERLU_MT
#include "slu_mt_zdefs.h"
#else
#include "slu_zdefs.h"
#endif

// since slu_zdefs.h and slu_ddefs.h can not be included
// at the same time, needed functions of slu_ddefs are copied here
// version : SuperLU 5.0, SuperLU_MT 3.0
extern "C"
{

  extern void
  dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
                         int *, int *, Stype_t, Dtype_t, Mtype_t);
  
  extern void
  dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
                       Stype_t, Dtype_t, Mtype_t);
    
#ifdef SELDON_WITH_SUPERLU_MT
  extern int_t  superlu_dQuerySpace (int_t, SuperMatrix *, SuperMatrix *, int_t, 
                                     superlu_memusage_t *);
  
  extern void pdgstrf (superlumt_options_t *, SuperMatrix *, int_t *, 
                       SuperMatrix *, SuperMatrix *, Gstat_t *, int_t *);

  extern void pdgstrf_init (int_t, fact_t, trans_t, yes_no_t, int_t, int_t, double, yes_no_t, double,
                            int_t *, int_t *, void *, int_t, SuperMatrix *,
                            SuperMatrix *, superlumt_options_t *, Gstat_t *);
  
  extern void dgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
                      int_t*, int_t*, SuperMatrix*, Gstat_t *, int_t *);  
#else
  extern void    dgstrf (superlu_options_t*, SuperMatrix*,
                         int, int, int*, void *, int, int *, int *, 
                         SuperMatrix *, SuperMatrix *, GlobalLU_t *,
                         SuperLUStat_t*, int *);
  
  extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                         SuperMatrix *, SuperLUStat_t*, int *);
  
  extern int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
#endif
  
}

#endif /* __SUPERLU_INTERFACE */
