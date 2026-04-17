#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// externs
extern SEXP get_cached_a_w(void);
extern SEXP get_cached_bb_w(void);
extern SEXP get_cached_a0(void);
extern SEXP get_cached_a1(void);
extern int wavelength_match(SEXP wavelengths);
extern SEXP c_build_cache(SEXP wavelengths);

// --- Helper function ---
double get_named_value(SEXP par, const char* name, int* found) {
  SEXP names = Rf_getAttrib(par, R_NamesSymbol);
  int n = Rf_length(par);

  for (int i = 0; i < n; i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      *found = 1;
      return REAL(par)[i];
    }
  }
  *found = 0;
  return 0.0;
}

// --- Main function ---
// [[register]]
SEXP c_iop_from_oac(SEXP wavelength_sexp, SEXP par_sexp) {
  int nprotect = 0;

  if (!Rf_isReal(wavelength_sexp))
    Rf_error("wavelength must be a numeric vector\n");
  if (!Rf_isReal(par_sexp) || Rf_isNull(Rf_getAttrib(par_sexp, R_NamesSymbol)))
    Rf_error("par must be a named numeric vector\n");

  // Coerce integer wavelength to real
  if (TYPEOF(wavelength_sexp) == INTSXP) {
    wavelength_sexp = PROTECT(Rf_coerceVector(wavelength_sexp, REALSXP));
    nprotect++;
  }

  // Make sure cache is up to date
  if (!wavelength_match(wavelength_sexp)) {
    c_build_cache(wavelength_sexp);
  }

  int n = Rf_length(wavelength_sexp);
  double* wl = REAL(wavelength_sexp);

  // Pointers to cached tables
  double* a0_ptr = REAL(get_cached_a0());
  double* a1_ptr = REAL(get_cached_a1());
  double* aw_ptr = REAL(get_cached_a_w());
  double* bb_w_ptr = REAL(get_cached_bb_w());

  // Prepare output
  SEXP a_out     = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  SEXP bb_out    = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  SEXP a_phy_out = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  SEXP a_dg_out  = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  double* a_res     = REAL(a_out);
  double* bb_res    = REAL(bb_out);
  double* a_phy_res = REAL(a_phy_out);
  double* a_dg_res  = REAL(a_dg_out);

  // Parse parameters
  int found;
  double chl = get_named_value(par_sexp, "chl", &found);
  int has_chl = found;

  double a_dg_440 = get_named_value(par_sexp, "a_dg_440", &found);
  int has_a_dg_440 = found;

  double a_nap_440 = get_named_value(par_sexp, "a_nap_440", &found);
  int has_a_nap_440 = found;

  double bb_p_550 = get_named_value(par_sexp, "bb_p_550", &found);
  int has_bb_p_550 = found;

  double a_dg_s = get_named_value(par_sexp, "a_dg_s", &found);
  int has_a_dg_s = found;

  double a_nap_s_d = get_named_value(par_sexp, "a_nap_s_d", &found);
  int has_a_nap_slope = found;

  double bb_p_gamma = get_named_value(par_sexp, "bb_p_gamma", &found);
  int has_bb_p_gamma = found;

  // Compute IOPs
  for (int i = 0; i < n; i++) {
    double wavelength = wl[i];

    // Phytoplankton absorption
    double a_phy = 0.0;
    if (has_chl) {
      double aph_440 = 0.06 * pow(chl, 0.65);
      double a0 = a0_ptr[i];
      double a1 = a1_ptr[i];
      a_phy = (a0 + a1 * log(aph_440)) * aph_440;
      if (a_phy < 0.0) a_phy = 0.0;
    }

    // CDOM absorption
    double a_g = 0.0;
    if (has_a_dg_440) {
      double slope = has_a_dg_s ? a_dg_s : 0.017;
      a_g = a_dg_440 * exp(-slope * (wavelength - 440.0));
    }

    // NAP absorption
    double a_nap = 0.0;
    if (has_a_nap_440) {
      double slope = has_a_nap_slope ? a_nap_s_d : 0.0116;
      a_nap = a_nap_440 * exp(-slope * (wavelength - 440.0));
    }

    // Particle backscattering
    double bb_p = 0.0;
    if (has_bb_p_550) {
      double gamma = has_bb_p_gamma ? bb_p_gamma : 0.46;
      bb_p = bb_p_550 * pow(wavelength / 550.0, -gamma);
    }

    // Total
    a_res[i]     = aw_ptr[i] + a_phy + a_g + a_nap;
    bb_res[i]    = bb_w_ptr[i] + bb_p;
    a_phy_res[i] = a_phy;
    a_dg_res[i]  = a_g + a_nap;
  }

  // Create output list (a, bb, a_phy, a_dg)
  // NOTE: existing callers use $a and $bb by name — adding elements is safe.
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4)); nprotect++;
  SET_VECTOR_ELT(out, 0, a_out);
  SET_VECTOR_ELT(out, 1, bb_out);
  SET_VECTOR_ELT(out, 2, a_phy_out);
  SET_VECTOR_ELT(out, 3, a_dg_out);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4)); nprotect++;
  SET_STRING_ELT(names, 0, Rf_mkChar("a"));
  SET_STRING_ELT(names, 1, Rf_mkChar("bb"));
  SET_STRING_ELT(names, 2, Rf_mkChar("a_phy"));
  SET_STRING_ELT(names, 3, Rf_mkChar("a_dg"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}
