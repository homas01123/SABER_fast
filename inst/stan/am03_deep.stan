// Stan model for SABER AM03 deep water inversion
// Implements Albert & Mobley (2003) forward model for optically deep waters
// No benthic reflectance component
//
// Author: SABER Development Team
// Date: 2026-02-16

functions {
  // Bio-optical model: Phytoplankton absorption from Chl
  vector calculate_a_phy(vector wavelength, real chl, vector a0, vector a1) {
    int N = num_elements(wavelength);
    vector[N] a_phy;
    
    real a_phy_443 = 0.06 * chl^0.65;
    
    for (i in 1:N) {
      real log_term = a0[i] + a1[i] * log(a_phy_443);
      a_phy[i] = fmax(0.0, log_term * a_phy_443);
    }
    
    return a_phy;
  }
  
  // Snell's law: Convert air angles to water angles
  real snell_law_angle(real theta_air_deg) {
    real n_air = 1.0;
    real n_water = 1.33;
    real theta_air_rad = theta_air_deg * pi() / 180.0;
    real theta_water_rad = asin((n_air / n_water) * sin(theta_air_rad));
    return theta_water_rad;
  }
  
  // AM03 forward model for deep water (VECTORIZED)
  vector forward_am03_deep(
      vector wavelength,
      real chl,
      real a_g_440,
      real bb_p_550,
      real a_g_s,
      real bb_p_gamma,
      vector a_w,
      vector bb_w,
      vector a0,
      vector a1,
      real theta_sun_deg,
      real theta_view_deg,
      int water_type
  ) {
    int N = num_elements(wavelength);
    
    // Refract angles into water
    real theta_sun_w = snell_law_angle(theta_sun_deg);
    real theta_view_w = snell_law_angle(theta_view_deg);
    real cos_sun = cos(theta_sun_w);
    real cos_view = cos(theta_view_w);
    
    // Bio-optical conversions (VECTORIZED)
    vector[N] a_phy = calculate_a_phy(wavelength, chl, a0, a1);
    vector[N] a_dg = a_g_440 * exp(-a_g_s * (wavelength - 440.0));  // Combined CDOM + detritus
    vector[N] bb_p = bb_p_550 * pow(wavelength / 550.0, -bb_p_gamma);
    
    // Total IOPs (VECTORIZED)
    vector[N] a_total = a_w + a_phy + a_dg;
    vector[N] bb_total = bb_w + bb_p;
    vector[N] ext = a_total + bb_total;
    vector[N] omega_b = bb_total ./ ext;
    
    // AM03 parametric coefficients (water type 2)
    real Prs1 = 0.0512;
    real Prs2 = 4.6659;
    real Prs3 = -7.8387;
    real Prs4 = 5.4571;
    real Prs5 = 0.1098;
    real Prs7 = 0.4021;
    
    // Deep water Rrs (VECTORIZED)
    vector[N] omega_b2 = omega_b .* omega_b;
    vector[N] omega_b3 = omega_b2 .* omega_b;
    vector[N] f_rs = Prs1 * 
                (1.0 + Prs2 * omega_b + Prs3 * omega_b2 + Prs4 * omega_b3) *
                (1.0 + Prs5 / cos_sun) *
                (1.0 + Prs7 / cos_view);
    vector[N] rrs_deep = f_rs .* omega_b;
    
    return rrs_deep;
  }
}

data {
  // Observations
  int<lower=1> n_wl;
  vector[n_wl] wavelength;
  vector[n_wl] rrs_obs;
  vector<lower=1e-12>[n_wl] sigma;
  
  // Fixed parameters
  int<lower=1,upper=2> water_type;
  real<lower=0,upper=90> theta_sun_deg;
  real<lower=0,upper=90> theta_view_deg;
  
  // Pre-computed lookup tables
  vector[n_wl] a_w;
  vector[n_wl] bb_w;
  vector[n_wl] a0;
  vector[n_wl] a1;
}

parameters {
  // Optically active constituents
  real<lower=0.01, upper=100> chl;
  real<lower=1e-6, upper=10> a_g_440;      // Combined CDOM + detritus at 440nm
  real<lower=1e-6, upper=1> bb_p_550;
  
  // Spectral slopes
  real<lower=0.001, upper=0.10> a_g_s;     // Combined CDOM/detritus slope (increased upper bound)
  real<lower=0.2, upper=3.0> bb_p_gamma;
}

transformed parameters {
  vector[n_wl] rrs_model;
  
  rrs_model = forward_am03_deep(
    wavelength, chl, a_g_440, bb_p_550,
    a_g_s, bb_p_gamma,
    a_w, bb_w, a0, a1,
    theta_sun_deg, theta_view_deg, water_type
  );
}

model {
  // Informative priors based on expected ranges
  chl ~ lognormal(log(4.5), 1.1);          // 95% CI: [0.5, 40]
  a_g_440 ~ lognormal(log(0.13), 0.75);    // 95% CI: [0.03, 0.6]
  bb_p_550 ~ lognormal(log(0.008), 0.5);   // 95% CI: [0.003, 0.02]
  
  a_g_s ~ lognormal(log(0.017), 0.28);     // 95% CI: [0.010, 0.030]
  bb_p_gamma ~ lognormal(log(0.45), 0.4);  // 95% CI: [0.2, 1.0]
  
  // Likelihood (sigma provided in data)
  rrs_obs ~ normal(rrs_model, sigma);
}

generated quantities {
  vector[n_wl] rrs_hat = rrs_model;
  real log_lik = normal_lpdf(rrs_obs | rrs_model, sigma);
}
