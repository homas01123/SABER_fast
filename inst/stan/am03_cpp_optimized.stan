// SABER AM03 Stan Model - C++ Optimized (Colleague's Approach)
// Uses external C++ functions from rtm_stan_funcs.hpp for faster performance
// 
// PERFORMANCE: ~3-5x faster than pure Stan due to Eigen optimization
// TRADEOFF: Requires C++ header, slightly more complex setup
//
// Author: Colleague (original), Integrated: 2026-01-15

functions {
  // External C++ function declarations (implementations in rtm_stan_funcs.hpp)
  // These must be declared here so stanc knows about them
  
  // IOP calculation - returns [n_wl, 2] matrix where col1=a, col2=bb
  matrix iop_from_oac_all(vector wavelength, vector a_w, vector a0, vector a1, 
                          vector bb_w, real chl, real a_g_440, real a_nap_440,
                          real a_g_s, real a_nap_s, real bb_p_550, real bb_p_gamma);
  
  // Forward model - returns vector[n_wl] of modeled Rrs
  vector forward_am03_ad(vector wavelength, vector a, vector bb, 
                         int water_type, real theta_sun_deg, real theta_view_deg,
                         int shallow, real h_w, vector r_b);
}

data {
  int<lower=1> n_wl;
  vector[n_wl] wavelength;
  vector[n_wl] rrs_obs;
  vector<lower=1e-12>[n_wl] sigma;

  // Pre-interpolated LUTs on the same wavelength grid
  vector[n_wl] a_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
  vector[n_wl] bb_w;

  // Bottom library on the same wavelength grid
  int<lower=1> n_benthic;
  matrix[n_wl, n_benthic] r_b_matrix;

  // Geometry / flags
  int<lower=1,upper=2> water_type;
  real theta_sun_deg;
  real theta_view_deg;
  int<lower=0,upper=1> shallow;
}

parameters {
  real<lower=0.01, upper=100> chl;
  real<lower=1e-6, upper=10> a_g_440;
  real<lower=1e-6, upper=10> a_nap_440;
  real<lower=1e-6, upper=1> bb_p_550;
  
  real<lower=0.001, upper=0.05> a_g_slope;
  real<lower=0.2, upper=3.0> bb_p_gamma;
  
  real<lower=0.1, upper=50> h_w;
  simplex[n_benthic] r_b_fractions;
  
  real<lower=1e-6> sigma_rrs;
}

transformed parameters {
  vector[n_wl] a;
  vector[n_wl] bb;
  vector[n_wl] r_b;
  vector[n_wl] rrs_model;

  // IOPs from OAC using C++ external function (FAST!)
  {
    real a_nap_s = 0.0116;  // Fixed NAP slope
    matrix[n_wl, 2] iop = iop_from_oac_all(
      wavelength, a_w, a0, a1, bb_w,
      chl, a_g_440, a_nap_440,
      a_g_slope, a_nap_s,
      bb_p_550, bb_p_gamma
    );
    a  = iop[, 1];
    bb = iop[, 2];
  }

  // Bottom reflectance mixture (Stan code - simple enough)
  r_b = r_b_matrix * r_b_fractions;

  // Forward RT model using C++ external function (FAST!)
  rrs_model = forward_am03_ad(
    wavelength, a, bb,
    water_type, theta_sun_deg, theta_view_deg,
    shallow, h_w, r_b
  );
}

model {
  // Priors (same as pure Stan version)
  chl ~ lognormal(log(1.0), 1.5);
  a_g_440 ~ lognormal(log(0.05), 1.5);
  a_nap_440 ~ lognormal(log(0.05), 1.5);
  bb_p_550 ~ lognormal(log(0.01), 1.5);
  
  a_g_slope ~ lognormal(log(0.017), 0.3);
  bb_p_gamma ~ lognormal(log(1.0), 0.5);
  
  h_w ~ lognormal(log(5.0), 1.0);
  
  sigma_rrs ~ normal(0, 0.01);

  // Likelihood
  rrs_obs ~ normal(rrs_model, sigma_rrs);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
  
  // Log-likelihood for model comparison
  real log_lik = normal_lpdf(rrs_obs | rrs_model, sigma_rrs);
}
