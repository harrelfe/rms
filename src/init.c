#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* Claude Sonnet 5 2026-08-30
   Added ormidx registration; updated ormll's argument count: 33 -> 38
   for the new ia, ia2, sgn, ib, nb arguments, then 38 -> 36 after
   removing y, y2 (no longer used in ormll's computation -- see
   ormll.f90 and the new ormidx.f90). Added ormeta registration,
   17 args, later extended to 19 (added s1, s2: per-observation score
   contributions credited separately to ia(i)/ia2(i), needed for the
   sparse missing-information correction to the final information
   matrix -- see orm.rfit.r's sparseInfoMatrix). */

/* .Fortran calls */
extern void F77_NAME(lrmll)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(matinv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(robcovf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormidx)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormeta)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormll)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"lrmll",   (DL_FUNC) &F77_NAME(lrmll),   19},
    {"matinv",  (DL_FUNC) &F77_NAME(matinv),  11},
    {"robcovf", (DL_FUNC) &F77_NAME(robcovf),  8},
    {"ormidx",  (DL_FUNC) &F77_NAME(ormidx),  10},
    {"ormeta",  (DL_FUNC) &F77_NAME(ormeta),  19},
    {"ormll",   (DL_FUNC) &F77_NAME(ormll),   36},
    {NULL, NULL, 0}
};

void R_init_rms(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
