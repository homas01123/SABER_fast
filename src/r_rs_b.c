#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <string.h>

extern SEXP get_cached_r_rs_b_colnames(void);
extern SEXP get_cached_r_rs_b(void);

// [[register]]
SEXP c_compute_r_rs_b_lmm(SEXP fractions) {
  int nprotect = 0;

  if (!Rf_isReal(fractions) || Rf_isNull(Rf_getAttrib(fractions, R_NamesSymbol)))
    Rf_error("fractions must be a named numeric vector\n");

  SEXP colnames = get_cached_r_rs_b_colnames();
  SEXP interp = get_cached_r_rs_b();
  if (interp == NULL)
    Rf_error("Cached bottom reflectance not initialized. Run c_build_cache() first.\n");

  int n_wl = Rf_nrows(interp);
  int n_class = Rf_ncols(interp);
  double* mat = REAL(interp);

  if (colnames == NULL || Rf_length(colnames) != n_class)
    Rf_error("Bottom class names missing or inconsistent\n");

  SEXP names = Rf_getAttrib(fractions, R_NamesSymbol);
  int n_frac = Rf_length(fractions);
  double* fracs = REAL(fractions);

  // --- NEW: Check if all fractions are NA ---
  int all_na = 1;
  for (int j = 0; j < n_frac; j++) {
    if (!ISNAN(fracs[j])) {
      all_na = 0;
      break;
    }
  }
  if (all_na) {
    Rf_warning("All r_rs_b_* values are NULL (NA); water type is optically deep.");
    return R_NilValue;
  }
  // --- END NEW ---

  // Normalize fractions
  double sum_frac = 0.0;
  for (int j = 0; j < n_frac; j++) {
    if (fracs[j] < 0.0)
      Rf_error("Fraction values must be non-negative\n");
    sum_frac += fracs[j];
  }

  if (sum_frac <= 0.0)
    Rf_error("Sum of fractions must be > 0 for normalization\n");

  SEXP result = PROTECT(Rf_allocVector(REALSXP, n_wl)); nprotect++;
  double* out = REAL(result);

  // Zero initialize
  for (int i = 0; i < n_wl; i++) out[i] = 0.0;

  for (int j = 0; j < n_frac; j++) {
    const char* name = CHAR(STRING_ELT(names, j));
    int matched = -1;

    for (int k = 0; k < n_class; k++) {
      const char* col = CHAR(STRING_ELT(colnames, k));
      if (strcmp(name, col) == 0) {
        matched = k;
        break;
      }
    }

    if (matched < 0) {
      Rf_error("Class name '%s' not found in bottom reflectance matrix\n", name);
    }

    // Normalize the fraction
    double weight = fracs[j] / sum_frac;

    for (int i = 0; i < n_wl; i++) {
      out[i] += weight * mat[i + matched * n_wl];
    }
  }

  UNPROTECT(nprotect);
  return result;
}
