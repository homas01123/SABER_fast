#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP c_load_pure_water(SEXP, SEXP);
extern SEXP c_load_a0_a1(SEXP, SEXP, SEXP);
extern SEXP c_pure_water_iop(SEXP);
extern SEXP c_iop_from_oac(SEXP, SEXP);
extern SEXP c_forward_am03(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_load_r_rs_b(SEXP, SEXP);
extern SEXP c_compute_r_rs_b_lmm(SEXP);
extern SEXP c_snell_law(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"c_load_pure_water", (DL_FUNC) &c_load_pure_water, 2},
  {"c_load_a0_a1", (DL_FUNC) &c_load_a0_a1, 3},
  {"c_pure_water_iop", (DL_FUNC) &c_pure_water_iop, 1},
  {"c_iop_from_oac", (DL_FUNC) &c_iop_from_oac, 2},
  {"c_forward_am03", (DL_FUNC) &c_forward_am03, 7},
  {"c_load_r_rs_b", (DL_FUNC) &c_load_r_rs_b, 2},
  {"c_compute_r_rs_b_lmm", (DL_FUNC) &c_compute_r_rs_b_lmm, 1},
  {"c_snell_law", (DL_FUNC) &c_snell_law, 2},
  {NULL, NULL, 0}
};

void R_init_SABER(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

