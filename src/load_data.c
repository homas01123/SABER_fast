#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// --- Static global memory ---

// Original spectral tables
static SEXP a_w_wavelength = NULL;
static SEXP a_w_value = NULL;

static SEXP a0_a1_wavelength = NULL;
static SEXP a0_value = NULL;
static SEXP a1_value = NULL;

static SEXP r_rs_b_wavelength = NULL;
static SEXP r_rs_b_matrix = NULL;
static SEXP r_rs_b_colnames = NULL;

// Cached interpolated data
static SEXP cached_wavelengths = NULL;
static SEXP cached_a0 = NULL;
static SEXP cached_a1 = NULL;
static SEXP cached_a_w = NULL;
static SEXP cached_bb_w = NULL;
static SEXP cached_r_rs_b = NULL;

// --- Loader functions ---

// [[register]]
SEXP c_load_pure_water(SEXP wl, SEXP a) {
  if (!Rf_isReal(wl) || !Rf_isReal(a) || Rf_length(wl) != Rf_length(a))
    Rf_error("Invalid pure water data\n");

  if (a_w_wavelength != NULL) R_ReleaseObject(a_w_wavelength);
  if (a_w_value != NULL) R_ReleaseObject(a_w_value);

  a_w_wavelength = wl;
  a_w_value = a;

  R_PreserveObject(a_w_wavelength);
  R_PreserveObject(a_w_value);

  return R_NilValue;
}

// [[register]]
SEXP c_load_a0_a1(SEXP wl, SEXP a0, SEXP a1) {
  if (!Rf_isReal(wl) || !Rf_isReal(a0) || Rf_length(wl) != Rf_length(a0))
    Rf_error("Invalid a0 data\n");
  if (!Rf_isReal(wl) || !Rf_isReal(a1) || Rf_length(wl) != Rf_length(a1))
    Rf_error("Invalid a1 data\n");

  if (a0_a1_wavelength != NULL) R_ReleaseObject(a0_a1_wavelength);
  if (a0_value != NULL) R_ReleaseObject(a0_value);
  if (a1_value != NULL) R_ReleaseObject(a1_value);

  a0_a1_wavelength = wl;
  a0_value = a0;
  a1_value = a1;

  R_PreserveObject(a0_a1_wavelength);
  R_PreserveObject(a0_value);
  R_PreserveObject(a1_value);

  return R_NilValue;
}

// [[register]]
SEXP c_load_r_rs_b(SEXP wavelength, SEXP matrix) {
  if (!Rf_isReal(wavelength) || !Rf_isMatrix(matrix))
    Rf_error("Expecting numeric wavelength and matrix\n");

  int nrow = Rf_nrows(matrix);
  if (Rf_length(wavelength) != nrow)
    Rf_error("wavelength length must match number of rows\n");

  SEXP dimnames = Rf_getAttrib(matrix, R_DimNamesSymbol);
  if (TYPEOF(dimnames) != VECSXP || LENGTH(dimnames) != 2)
    Rf_error("matrix must have dimnames\n");

  SEXP colnames = VECTOR_ELT(dimnames, 1);
  if (!Rf_isString(colnames))
    Rf_error("column names must be character vector\n");

  // --- Free old objects safely ---
  if (r_rs_b_wavelength != NULL) R_ReleaseObject(r_rs_b_wavelength);
  if (r_rs_b_matrix != NULL) R_ReleaseObject(r_rs_b_matrix);
  if (r_rs_b_colnames != NULL) R_ReleaseObject(r_rs_b_colnames);

  // --- Assign and preserve new ones ---
  r_rs_b_wavelength = wavelength;
  r_rs_b_matrix = matrix;
  r_rs_b_colnames = colnames;

  R_PreserveObject(r_rs_b_wavelength);
  R_PreserveObject(r_rs_b_matrix);
  R_PreserveObject(r_rs_b_colnames);

  return R_NilValue;
}


// --- scalar linear interpolation ---
double linear_interpolation(SEXP wl_vec, SEXP val_vec, double wl) {
  int n = Rf_length(wl_vec);
  double* wl_arr = REAL(wl_vec);
  double* val_arr = REAL(val_vec);

  if (wl <= wl_arr[0] || wl >= wl_arr[n-1])
    return 0.0;

  for (int i = 0; i < n - 1; i++) {
    if (wl >= wl_arr[i] && wl <= wl_arr[i+1]) {
      double slope = (val_arr[i+1] - val_arr[i]) / (wl_arr[i+1] - wl_arr[i]);
      return val_arr[i] + slope * (wl - wl_arr[i]);
    }
  }
  Rf_error("Interpolation failed\n");
  return NA_REAL;
}

// --- vector linear interpolation ---
SEXP interpolate_vector(SEXP wl_old, SEXP val_old, SEXP wl_new, SEXP data_name_sexp) {
  if (!Rf_isReal(wl_old) || !Rf_isReal(val_old) || !Rf_isReal(wl_new))
    Rf_error("Inputs must be numeric vectors\n");

  if (!Rf_isString(data_name_sexp) || Rf_length(data_name_sexp) != 1)
    Rf_error("data_name must be a single string\n");

  int n_old = Rf_length(wl_old);
  int n_new = Rf_length(wl_new);

  double* wl_old_ptr = REAL(wl_old);
  double* wl_new_ptr = REAL(wl_new);

  SEXP result = PROTECT(Rf_allocVector(REALSXP, n_new));
  double* val_new_ptr = REAL(result);

  int extrapolated = 0;
  double wl_min = wl_old_ptr[0];
  double wl_max = wl_old_ptr[n_old - 1];

  for (int i = 0; i < n_new; i++) {
    double wl = wl_new_ptr[i];
    val_new_ptr[i] = linear_interpolation(wl_old, val_old, wl);
    if (wl <= wl_min || wl >= wl_max)
      extrapolated = 1;
  }

  if (extrapolated) {
    const char* data_name = CHAR(STRING_ELT(data_name_sexp, 0));
    Rf_warning("Extrapolation detected in %s: requested wavelengths outside [%.1f, %.1f] nm; values set to 0\n",
               data_name, wl_min, wl_max);
  }

  UNPROTECT(1);
  return result;
}

// --- matrix linear interpolation ---
SEXP interpolate_matrix(SEXP wl_old, SEXP mat_old, SEXP wl_new, SEXP data_name_sexp) {
  if (!Rf_isReal(wl_old) || !Rf_isReal(wl_new) || !Rf_isMatrix(mat_old))
    Rf_error("Inputs must be numeric vectors and matrix\n");

  if (!Rf_isString(data_name_sexp) || Rf_length(data_name_sexp) != 1)
    Rf_error("data_name must be a single string\n");

  int n_old = Rf_length(wl_old);
  int n_new = Rf_length(wl_new);
  int n_class = Rf_ncols(mat_old);

  double* wl_old_ptr = REAL(wl_old);
  double* wl_new_ptr = REAL(wl_new);
  double* mat_old_ptr = REAL(mat_old);

  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, n_new, n_class));
  double* mat_new_ptr = REAL(result);

  int extrapolated = 0;
  double wl_min = wl_old_ptr[0];
  double wl_max = wl_old_ptr[n_old - 1];

  // for (int j = 0; j < n_class; j++) {
  //   for (int i = 0; i < n_new; i++) {
  //     double wl = wl_new_ptr[i];
  //
  //     mat_new_ptr[i + j * n_new] = linear_interpolation(wl_old, VECTOR_ELT(mat_old, j), wl);
  //     if (wl <= wl_min || wl >= wl_max)
  //       extrapolated = 1;
  //   }
  // }

  for (int j = 0; j < n_class; j++) {
    for (int i = 0; i < n_new; i++) {
      double wl = wl_new_ptr[i];

      // Extract column slice for interpolation (pointer to column j)
      SEXP temp_col = PROTECT(Rf_allocVector(REALSXP, n_old));
      double* temp_col_ptr = REAL(temp_col);

      for (int k = 0; k < n_old; k++) {
        temp_col_ptr[k] = mat_old_ptr[k + j * n_old];
      }

      mat_new_ptr[i + j * n_new] = linear_interpolation(wl_old, temp_col, wl);
      UNPROTECT(1);

      if (wl <= wl_min || wl >= wl_max)
        extrapolated = 1;
    }
  }

  if (extrapolated) {
    const char* data_name = CHAR(STRING_ELT(data_name_sexp, 0));
    Rf_warning("Extrapolation detected in %s: requested wavelengths outside [%.1f, %.1f] nm; values set to 0\n",
               data_name, wl_min, wl_max);
  }

  UNPROTECT(1);
  return result;
}

// --- Cache Management ---

int wavelength_match(SEXP wl_request) {
  if (cached_wavelengths == NULL) return 0;

  int n_cached = Rf_length(cached_wavelengths);
  int n_req = Rf_length(wl_request);
  if (n_cached != n_req) return 0;

  double* cached = REAL(cached_wavelengths);
  double* req = REAL(wl_request);

  for (int i = 0; i < n_cached; i++) {
    if (fabs(cached[i] - req[i]) > 1e-8)
      return 0;
  }
  return 1;
}

// [[register]]
SEXP c_build_cache(SEXP wavelengths) {
  if (a_w_wavelength == NULL || a_w_value == NULL ||
      a0_a1_wavelength == NULL || a0_value == NULL || a1_value == NULL ||
      r_rs_b_wavelength == NULL || r_rs_b_matrix == NULL) {
    Rf_error("Spectral reference data not fully loaded. Call SABER::.onLoad() to initialize datasets.\n");
  }

  int nprotect = 0;

  if (!Rf_isReal(wavelengths) && !Rf_isInteger(wavelengths))
    Rf_error("wavelength must be numeric or integer vector\n");

  if (TYPEOF(wavelengths) == INTSXP)
    wavelengths = Rf_coerceVector(wavelengths, REALSXP);

  int n = Rf_length(wavelengths);
  int m = Rf_ncols(r_rs_b_matrix);  // number of bottom classes
  double* wl = REAL(wavelengths);

  // --- Release old cache ---
  if (cached_wavelengths != NULL) R_ReleaseObject(cached_wavelengths);
  if (cached_a_w != NULL) R_ReleaseObject(cached_a_w);
  if (cached_bb_w != NULL) R_ReleaseObject(cached_bb_w);
  if (cached_a0 != NULL) R_ReleaseObject(cached_a0);
  if (cached_a1 != NULL) R_ReleaseObject(cached_a1);
  if (cached_r_rs_b != NULL) R_ReleaseObject(cached_r_rs_b);

  // --- Allocate and interpolate ---
  cached_wavelengths = wavelengths;
  R_PreserveObject(cached_wavelengths);

  cached_a_w = PROTECT(interpolate_vector(a_w_wavelength, a_w_value, wavelengths, Rf_mkString("pure water absorption"))); nprotect++;

  cached_a0 = PROTECT(interpolate_vector(a0_a1_wavelength, a0_value, wavelengths, Rf_mkString("phytoplankton a0"))); nprotect++;

  cached_a1 = PROTECT(interpolate_vector(a0_a1_wavelength, a1_value, wavelengths, Rf_mkString("phytoplankton a1"))); nprotect++;

  cached_bb_w = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  double* bb_w_ptr = REAL(cached_bb_w);

  double b1 = 0.00111;
  double lambda1 = 500.0;
  double exponent = -4.32;

  for (int i = 0; i < n; i++) {
    bb_w_ptr[i] = b1 * pow(wl[i] / lambda1, exponent);
  }

  cached_r_rs_b = PROTECT(interpolate_matrix(
    r_rs_b_wavelength, r_rs_b_matrix, wavelengths, Rf_mkString("bottom reflectance"))); nprotect++;

    // --- Preserve all ---
    R_PreserveObject(cached_a_w);
    R_PreserveObject(cached_bb_w);
    R_PreserveObject(cached_a0);
    R_PreserveObject(cached_a1);
    R_PreserveObject(cached_r_rs_b);

    UNPROTECT(nprotect);
    return R_NilValue;
}

// Reset the interpolated cache so the next wavelength-keyed call rebuilds it.
// Called after c_load_r_rs_b replaces the benthic matrix.
SEXP c_reset_cache(void) {
  if (cached_wavelengths != NULL) {
    R_ReleaseObject(cached_wavelengths);
    cached_wavelengths = NULL;
  }
  if (cached_r_rs_b != NULL) {
    R_ReleaseObject(cached_r_rs_b);
    cached_r_rs_b = NULL;
  }
  return R_NilValue;
}

SEXP get_cached_a_w() {
  return cached_a_w;
}

SEXP get_cached_bb_w() {
  return cached_bb_w;
}

SEXP get_cached_a0() {
  return cached_a0;
}

SEXP get_cached_a1() {
  return cached_a1;
}

SEXP get_cached_r_rs_b_colnames() {
  return r_rs_b_colnames;
}

SEXP get_cached_r_rs_b() {
  return cached_r_rs_b;
}
