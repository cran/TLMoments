#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP TLMoments_pseudo_C(SEXP, SEXP);
extern SEXP TLMoments_pwm_C(SEXP, SEXP);
extern SEXP TLMoments_PWM_to_TLMoments(SEXP, SEXP, SEXP);
extern SEXP TLMoments_TLMoment_direct(SEXP, SEXP, SEXP, SEXP);
extern SEXP TLMoments_TLMoment_PWM(SEXP, SEXP, SEXP, SEXP);
extern SEXP TLMoments_TLMoments_recurrence(SEXP, SEXP, SEXP, SEXP);
extern SEXP TLMoments_TLMoments_recursive(SEXP, SEXP, SEXP, SEXP);
extern SEXP TLMoments_z_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP TLMoments_Z_C(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"TLMoments_pseudo_C",             (DL_FUNC) &TLMoments_pseudo_C,             2},
  {"TLMoments_pwm_C",                (DL_FUNC) &TLMoments_pwm_C,                2},
  {"TLMoments_PWM_to_TLMoments",     (DL_FUNC) &TLMoments_PWM_to_TLMoments,     3},
  {"TLMoments_TLMoment_direct",      (DL_FUNC) &TLMoments_TLMoment_direct,      4},
  {"TLMoments_TLMoment_PWM",         (DL_FUNC) &TLMoments_TLMoment_PWM,         4},
  {"TLMoments_TLMoments_recurrence", (DL_FUNC) &TLMoments_TLMoments_recurrence, 4},
  {"TLMoments_TLMoments_recursive",  (DL_FUNC) &TLMoments_TLMoments_recursive,  4},
  {"TLMoments_z_C",                  (DL_FUNC) &TLMoments_z_C,                  4},
  {"TLMoments_Z_C",                  (DL_FUNC) &TLMoments_Z_C,                  3},
  {NULL, NULL, 0}
};

void R_init_TLMoments(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
