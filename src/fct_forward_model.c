#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

extern SEXP c_snell_law(SEXP, SEXP);

// Helper: convert degrees to radians
static inline double deg2rad(double deg) {
  return deg * M_PI / 180.0;
}

// [[register]]
SEXP c_forward_am03(SEXP wavelength_sexp, SEXP iop_list, SEXP water_type_sexp,
                    SEXP theta_sun_sexp, SEXP theta_view_sexp,
                    SEXP h_w_sexp, SEXP r_b_sexp) {

  int nprotect = 0;

  if (!Rf_isReal(wavelength_sexp)) Rf_error("wavelength must be numeric");
  if (!Rf_isNewList(iop_list)) Rf_error("iop must be a list with elements 'a' and 'bb'");
  if (!Rf_isReal(theta_sun_sexp) || !Rf_isReal(theta_view_sexp))
    Rf_error("theta_sun and theta_view must be numeric scalars");

  int n = LENGTH(wavelength_sexp);
  double* wl = REAL(wavelength_sexp);

  SEXP a_vec = VECTOR_ELT(iop_list, 0);
  SEXP bb_vec = VECTOR_ELT(iop_list, 1);
  if (LENGTH(a_vec) != n || LENGTH(bb_vec) != n)
    Rf_error("iop$a and iop$bb must be same length as wavelength");

  double* a = REAL(a_vec);
  double* bb = REAL(bb_vec);

  int water_type = INTEGER(water_type_sexp)[0]; // 1 or 2
  double theta_sun = deg2rad(REAL(theta_sun_sexp)[0]);
  double theta_view = deg2rad(REAL(theta_view_sexp)[0]);

  // Optional shallow parameters
  int shallow = (!Rf_isNull(h_w_sexp) && !Rf_isNull(r_b_sexp));
  double* r_b = NULL;
  double h_w = 0.0;

  if (shallow) {
    if (!Rf_isReal(h_w_sexp) || !Rf_isReal(r_b_sexp))
      Rf_error("h_w and r_b must be numeric if provided");
    if (LENGTH(r_b_sexp) != n)
      Rf_error("r_b must match length of wavelength");

    r_b = REAL(r_b_sexp);
    h_w = REAL(h_w_sexp)[0];
  }

  // Output vector
  SEXP rrs_out = PROTECT(Rf_allocVector(REALSXP, n)); nprotect++;
  double* rrs = REAL(rrs_out);

  SEXP geometry = PROTECT(c_snell_law(theta_view_sexp, theta_sun_sexp)); nprotect++;
  double view_w = REAL(VECTOR_ELT(geometry, 0))[0];
  double sun_w  = REAL(VECTOR_ELT(geometry, 1))[0];

  for (int i = 0; i < n; i++) {
    double ext = a[i] + bb[i];
    double omega_b = bb[i] / ext;

    // Fresnel & geometric effects (simplified snell model)
    double f_rs;
    if (water_type == 1) {
      f_rs = 0.095;
    } else {
      f_rs = 0.0512 *
        (1 + (4.6659 * omega_b) +
        (-7.8387 * omega_b * omega_b) +
        (5.4571 * omega_b * omega_b * omega_b)) *
        (1 + (0.1098 / cos(sun_w))) *
        (1 + (0.4021 / cos(view_w)));
    }

    double rrs_deep = f_rs * omega_b;

    if (shallow) {
      double k0 = (water_type == 1) ? 1.0395 : 1.0546;

      double Kd = k0 * (ext / cos(sun_w));
      double kuW = (ext / cos(view_w)) *
        pow(1 + omega_b, 3.5421) *
        (1 - (0.2786 / cos(sun_w)));

      double kuB = (ext / cos(view_w)) *
        pow(1 + omega_b, 2.2658) *
        (1 - (0.0577 / cos(sun_w)));

      double Ars1 = 1.1576;
      double Ars2 = 1.0389;

      rrs[i] = rrs_deep * (1 - (Ars1 * exp(-h_w * (Kd + kuW)))) +
        Ars2 * r_b[i] * exp(-h_w * (Kd + kuB));
    } else {
      rrs[i] = rrs_deep;
    }
  }

  UNPROTECT(nprotect);
  return rrs_out;
}

// ===========================================================================
// Semi-analytical SICF Rrs component (Gilerson et al. 2007 + dual-Gaussian)
//
// Computes the fluorescence contribution to Rrs given a pre-computed
// subsurface downwelling irradiance Ed_0m.  The Ed caching is handled
// on the R side; this function only runs the per-iteration math.
//
// Arguments:
//   wavelength_sexp  – double[n] wavelengths (nm)
//   Ed_0m_sexp       – double[n] subsurface downwelling Ed (W m-2 nm-1)
//   chl_sexp         – double scalar  chlorophyll-a (mg m-3)
//   a_dg_443_sexp    – double scalar  CDOM+NAP absorption at 443 nm (m-1)
//   phi_f_sexp       – double scalar  fluorescence quantum yield
//   coeff_sexp       – double[3]  Gilerson coefficients [c1, c2, c3]
//   scale_sexp       – double scalar  Lf_685 scale factor (default 13.5)
//
// Returns: double[n] Rrs_sicf (sr-1)
// ===========================================================================
SEXP c_sicf_rrs_semi_analytical(SEXP wavelength_sexp, SEXP Ed_0m_sexp,
                                 SEXP chl_sexp,  SEXP a_dg_443_sexp,
                                 SEXP phi_f_sexp, SEXP coeff_sexp,
                                 SEXP scale_sexp) {

  int n      = LENGTH(wavelength_sexp);
  double *wl = REAL(wavelength_sexp);
  double *Ed = REAL(Ed_0m_sexp);

  double chl   = REAL(chl_sexp)[0];
  double adg   = REAL(a_dg_443_sexp)[0];
  double phi_f = REAL(phi_f_sexp)[0];
  double *coef = REAL(coeff_sexp);   /* coef[0]=c1, [1]=c2, [2]=c3 */
  double scale = REAL(scale_sexp)[0];

  /* ---------- Gilerson Lf_685 scalar ---------- */
  double lf_height = phi_f * coef[0] * chl /
                     (1.0 + coef[1] * adg + coef[2] * chl);

  /* Find index closest to 685 nm */
  int idx685 = 0;
  double min_d = fabs(wl[0] - 685.0);
  for (int i = 1; i < n; i++) {
    double d = fabs(wl[i] - 685.0);
    if (d < min_d) { min_d = d; idx685 = i; }
  }
  double Lf_685 = lf_height * Ed[idx685] / scale;

  /* ---------- Dual-Gaussian spectral shape ---------- */
  /* Primary:   centre=685 nm, FWHM=25 nm  → σ² = 25²/(4 ln2) */
  /* Secondary: centre=730 nm, FWHM=50 nm  → σ² = 50²/(4 ln2), amp=0.3 */
  static const double LN2    = 0.6931471805599453;
  double k1 = -4.0 * LN2 / (25.0 * 25.0);   /* = -4ln2/FWHM₁² */
  double k2 = -4.0 * LN2 / (50.0 * 50.0);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *rrs_sicf = REAL(out);

  for (int i = 0; i < n; i++) {
    double d1 = wl[i] - 685.0;
    double d2 = wl[i] - 730.0;
    double shape = exp(k1 * d1 * d1) + 0.3 * exp(k2 * d2 * d2);
    rrs_sicf[i] = Lf_685 * shape / Ed[i];
  }

  UNPROTECT(1);
  return out;
}

// ===========================================================================
// WRF Matrix Discretization
//
// Builds the N×N wavelength-redistribution-function matrix for Chl
// fluorescence.  Exploits the separability of the integrand:
//
//   WRF[i,j] = factors[j] × φ_f × A[i] × B[j]   (i < j only)
//
//   A[i] = Σ_{wex in band_i}  gchl(wex) × wex          excitation sum
//   B[j] = Σ_{wem in band_j}  h(wem) / wem              emission sum
//   h(λ) = dual-Gaussian emission shape (peaks at 685 & 730 nm)
//   gchl(λ) = 1 for 370 ≤ λ ≤ 690 nm, else 0
//
// Arguments:
//   wavelength_sexp – double[N]  wavelength band centres (nm)
//   phi_f_sexp      – double scalar  quantum yield
//
// Returns: double[N,N] WRF matrix (column-major)
// ===========================================================================
SEXP c_discretize_wrf(SEXP wavelength_sexp, SEXP phi_f_sexp) {

  int    N     = LENGTH(wavelength_sexp);
  double *wl   = REAL(wavelength_sexp);
  double phi_f = REAL(phi_f_sexp)[0];

  /* Emission-shape constants matching .wrf_chl() in R */
  static const double SQRT2LN2 = 1.17741002251547;   /* sqrt(2*ln2) */
  static const double PI_      = 3.14159265358979323846;
  double sigma1 = 25.0 / (2.0 * SQRT2LN2);
  double sigma2 = 50.0 / (2.0 * SQRT2LN2);
  double f1  = 1.0 / (sigma1 * sqrt(2.0 * PI_));
  double f2  = 1.0 / (sigma2 * sqrt(2.0 * PI_));
  double kk1 = 0.5 / (sigma1 * sigma1);
  double kk2 = 0.5 / (sigma2 * sigma2);

  static const double deltaw = 1.0;

  /* Temporary arrays (freed at next R GC checkpoint) */
  double *A       = (double *) R_alloc(N, sizeof(double));
  double *B       = (double *) R_alloc(N, sizeof(double));
  double *factors = (double *) R_alloc(N, sizeof(double));

  for (int i = 0; i < N; i++) {
    double bw = (i == 0) ? deltaw : (wl[i] - wl[i - 1]);
    int    ns = (int)(bw / deltaw);
    if (ns < 1) ns = 1;
    factors[i] = (deltaw * deltaw) / bw;

    double wc  = wl[i] + 0.5 * deltaw;   /* band centre */
    double suA = 0.0, suB = 0.0;
    for (int s = 0; s < ns; s++) {
      double w = wc + s * deltaw;
      /* A: excitation weight */
      if (w >= 370.0 && w <= 690.0) suA += w;
      /* B: emission weight  h(w)/w */
      double d1 = w - 685.0, d2 = w - 730.0;
      double h  = f1 * exp(-kk1 * d1 * d1) + 0.3 * f2 * exp(-kk2 * d2 * d2);
      suB += h / w;
    }
    A[i] = suA;
    B[i] = suB;
  }

  /* Allocate and zero output matrix */
  SEXP out  = PROTECT(Rf_allocMatrix(REALSXP, N, N));
  double *W = REAL(out);
  for (int k = 0; k < N * N; k++) W[k] = 0.0;

  /* Fill upper triangle: WRF[i,j] = factors[j]*phi_f*A[i]*B[j]  (i < j) */
  for (int j = 1; j < N; j++) {
    if (B[j] == 0.0) continue;
    double vj  = factors[j] * phi_f * B[j];
    double *col = W + j * N;
    for (int i = 0; i < j; i++) {
      col[i] = vj * A[i];
    }
  }

  UNPROTECT(1);
  return out;
}

// ===========================================================================
// Depth-Integrated SICF Fluorescence Radiance
//
// Computes Lf(λ_em):
//   Lf[j] = (1/4π) Σ_i { E0[i] × a_phy[i] × WRF[i,j] / (Kd[i] + Kd[j]) }
//
// Replaces four sequential N×N R matrix operations with a single fused C
// loop: no temporary matrix allocations, single memory pass, ~50× faster
// than the R vectorised version for N=800 wavelengths.
//
// Arguments:
//   Kd_sexp    – double[N]   diffuse attenuation (m⁻¹)
//   a_phy_sexp – double[N]   phytoplankton absorption (m⁻¹)
//   irrad_sexp – double[N]   scalar irradiance E0 at 0- (W m⁻² nm⁻¹)
//   WRF_sexp   – double[N,N] WRF matrix column-major (output of c_discretize_wrf)
//
// Returns: double[N]  Lf vector (W m⁻² nm⁻¹ sr⁻¹)
// ===========================================================================
SEXP c_sicf_depth_integrated(SEXP Kd_sexp, SEXP a_phy_sexp,
                              SEXP irrad_sexp, SEXP WRF_sexp,
                              SEXP phi_f_sexp) {

  int     N    = LENGTH(Kd_sexp);
  double *Kd   = REAL(Kd_sexp);
  double *aphy = REAL(a_phy_sexp);
  double *E0   = REAL(irrad_sexp);
  double *WRF  = REAL(WRF_sexp);   /* N×N column-major: WRF[i,j] = WRF[i + j*N] */
  double  phi_f = REAL(phi_f_sexp)[0];

  static const double INV4PI = 1.0 / (4.0 * 3.14159265358979323846);
  double  scale = phi_f * INV4PI;

  SEXP Lf_sexp = PROTECT(Rf_allocVector(REALSXP, N));
  double *Lf   = REAL(Lf_sexp);

  for (int j = 0; j < N; j++) {
    double lf_j  = 0.0;
    double kd_j  = Kd[j];
    double *colj = WRF + j * N;    /* stride-1 access down column j */
    for (int i = 0; i < N; i++) {
      double w = colj[i];
      if (w > 0.0)
        lf_j += E0[i] * aphy[i] * w / (Kd[i] + kd_j);
    }
    Lf[j] = lf_j * scale;
  }

  UNPROTECT(1);
  return Lf_sexp;
}
