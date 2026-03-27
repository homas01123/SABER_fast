// SABER AM03 Stan Model - BEST OF BOTH WORLDS
// Combines:
// - Colleague's C++ IOP calculation (faster, better cache efficiency)
// - Our vectorized forward model (readable, easier to modify)
//
// PERFORMANCE: ~2-3x faster than pure Stan (IOP is main bottleneck)
// FLEXIBILITY: Forward model in Stan = easy to modify/experiment
//
// Author: SABER Team, 2026-01-15

functions {
  // External C++ function for IOP calculation (from colleague)
  matrix iop_from_oac_all(
      vector wavelength,
      vector a_w,
      vector a0,
      vector a1,
      vector bb_w,
      real chl,
      real a_g_440,
      real a_nap_440,
      real a_g_s,
      real a_nap_s,
      real bb_p_550,
      real bb_p_gamma
  );
  
  // Our own Stan forward model (vectorized)
  // This is in Stan language for maximum flexibility
  real snell_law_angle(real theta_air_deg) {
    real n_air = 1.0;
    real n_water = 1.33;
    real theta_air_rad = theta_air_deg * pi() / 180.0;
    real theta_water_rad = asin((n_air / n_water) * sin(theta_air_rad));
    return theta_water_rad;
  }
  
  vector forward_am03_vectorized(
      vector wavelength,
      vector a,
      vector bb,
      int water_type,
      real theta_sun_deg,
      real theta_view_deg,
      int shallow,
      real h_w,
      vector r_b
  ) {
    int N = num_elements(wavelength);
    
    // Refract angles into water
    real theta_sun_w = snell_law_angle(theta_sun_deg);
    real theta_view_w = snell_law_angle(theta_view_deg);
    real cos_sun = cos(theta_sun_w);
    real cos_view = cos(theta_view_w);
    
    // Total extinction and single-scattering albedo (VECTORIZED)
    vector[N] ext = a + bb;
    vector[N] omega_b = bb ./ ext;
    
    // AM03 parametric coefficients (water type 2)
    real Prs1 = 0.0512;
    real Prs2 = 4.6659;
    real Prs3 = -7.8387;
    real Prs4 = 5.4571;
    real Prs5 = 0.1098;
    real Prs7 = 0.4021;
    real kappa0 = 1.0546;
    real Ars1 = 1.1576;
    real Ars2 = 1.0389;
    
    // Deep water Rrs (VECTORIZED)
    vector[N] omega_b2 = omega_b .* omega_b;
    vector[N] omega_b3 = omega_b2 .* omega_b;
    vector[N] f_rs = Prs1 * 
                (1.0 + Prs2 * omega_b + Prs3 * omega_b2 + Prs4 * omega_b3) *
                (1.0 + Prs5 / cos_sun) *
                (1.0 + Prs7 / cos_view);
    vector[N] rrs_deep = f_rs .* omega_b;
    
    if (shallow == 0) {
      return rrs_deep;
    }
    
    // Shallow water attenuation (VECTORIZED)
    vector[N] Kd = kappa0 * (ext / cos_sun);
    vector[N] kuW = (ext / cos_view) .* 
               pow(1.0 + omega_b, 3.5421) * 
               (1.0 - 0.2786 / cos_sun);
    vector[N] kuB = (ext / cos_view) .* 
               pow(1.0 + omega_b, 2.2658) * 
               (1.0 - 0.0577 / cos_sun);
    
    // Combined signal (VECTORIZED)
    vector[N] exp_water = exp(-h_w * (Kd + kuW));
    vector[N] exp_benthic = exp(-h_w * (Kd + kuB));
    vector[N] rrs = rrs_deep .* (1.0 - Ars1 * exp_water) + 
                    Ars2 * r_b .* exp_benthic;
    
    return rrs;
  }
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
    real a_nap_s = 0.0116;
    matrix[n_wl, 2] iop = iop_from_oac_all(
      wavelength, a_w, a0, a1, bb_w,
      chl, a_g_440, a_nap_440,
      a_g_slope, a_nap_s,
      bb_p_550, bb_p_gamma
    );
    a  = iop[, 1];
    bb = iop[, 2];
  }

  // Bottom reflectance mixture
  r_b = r_b_matrix * r_b_fractions;

  // Forward RT model using OUR vectorized Stan code (FLEXIBLE!)
  rrs_model = forward_am03_vectorized(
    wavelength, a, bb,
    water_type, theta_sun_deg, theta_view_deg,
    shallow, h_w, r_b
  );
}

model {
  // Priors
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
  real log_lik = normal_lpdf(rrs_obs | rrs_model, sigma_rrs);
}
