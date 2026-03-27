#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

// External accessors from load_data.c
extern SEXP get_cached_a_w(void);
extern SEXP get_cached_bb_w(void);
extern SEXP get_cached_a0(void);
extern SEXP get_cached_a1(void);
extern SEXP get_cached_r_rs_b(void);
extern SEXP get_cached_r_rs_b_colnames(void);
extern int wavelength_match(SEXP wavelengths);
extern SEXP c_build_cache(SEXP wavelengths);

// [[register]]
SEXP c_get_stan_data(SEXP wavelength_sexp) {
  if (!Rf_isReal(wavelength_sexp))
    Rf_error("wavelength must be a numeric vector");

  // Build cache if needed (same as c_pure_water_iop does)
  if (!wavelength_match(wavelength_sexp)) {
    c_build_cache(wavelength_sexp);
  }

  // Retrieve all cached data
  SEXP a_w = get_cached_a_w();
  SEXP bb_w = get_cached_bb_w();
  SEXP a0 = get_cached_a0();
  SEXP a1 = get_cached_a1();
  SEXP r_rs_b = get_cached_r_rs_b();
  SEXP colnames = get_cached_r_rs_b_colnames();

  if (a_w == NULL || a0 == NULL || r_rs_b == NULL)
    Rf_error("Cache not initialized. This should not happen.");

  // Create list output
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, a_w);
  SET_VECTOR_ELT(out, 1, bb_w);
  SET_VECTOR_ELT(out, 2, a0);
  SET_VECTOR_ELT(out, 3, a1);
  SET_VECTOR_ELT(out, 4, r_rs_b);
  SET_VECTOR_ELT(out, 5, colnames);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(names, 0, Rf_mkChar("a_w"));
  SET_STRING_ELT(names, 1, Rf_mkChar("bb_w"));
  SET_STRING_ELT(names, 2, Rf_mkChar("a0"));
  SET_STRING_ELT(names, 3, Rf_mkChar("a1"));
  SET_STRING_ELT(names, 4, Rf_mkChar("r_b_matrix"));
  SET_STRING_ELT(names, 5, Rf_mkChar("r_b_colnames"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}
