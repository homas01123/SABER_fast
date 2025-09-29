---
editor_options: 
  markdown: 
    wrap: 72
---

## Semi-Analytical Bayesian Estimate Retrieval (S.A.B.E.R.) Version 0.1

# Theoretical Background

SABER is a software developed in R to solve forward and inverse problems
of Ocean Color Remote Sensing (OCRS) across different types of oceanic
waters (CASE-I and CASE-II) as well as for inland waters by solving the
Scalar Radiative Transfer Equation (SRTE) using Semi-analytical (SA)
paramterizations.

-   FORWARD MODEL

The forward model calculates Remote-sensing reflectance
($R_\mathrm{rs}$) from user provided spectral IOPs ($a$ &
$b_\mathrm{b}$) using the SA paramterization of SRTE by Albert & Mobley
(2003) as default. The wavelengths for IOPs can be supplied in default
sensor acquisition wavelengths. In absence of IOPs, *in vivo*
Bio-geo-chemical (BGC) quantities; e.g. Chlorophyll concentration
$[chl]$, absorption of Colored Detrital Matter (CDM);
$a_\mathrm{dg}(\lambda = 443)$nm along with particulate backscatter;
$b_\mathrm{bP}(\lambda = 555)$nm can be used along with viewing geometry
for both the cases to simulate $R_\mathrm{rs}(\lambda)$. In the latter
case, *state-of-the-art* Bio-optical relationships are used to determine
BGC component specific $a$ & $b_\mathrm{b}$. In case of shallow waters,
where the benthic substrate has substantial upweling optical signal, the
bottom depth ($H$) and benthic reflectance ($R_\mathrm{B}$) can also be
supplied along with previously mentioned quantities to simulate
$R_\mathrm{rs}(\lambda)$ in shallow optically complex waters. The
forward model can be simply expressed as
$R_\mathrm{rs}(\lambda)=f(a(\lambda),b_\mathrm{b}(\lambda),H, R_\mathrm{B}(\lambda), \theta_\mathrm{v},\Delta\phi_\mathrm{v}, \theta_\mathrm{s})$
where, $a(\lambda)$ and $b_\mathrm{b}(\lambda)$ are absorption and
backscatter varying with wavelength $\lambda$, $H$ is bottom depth,
$R_\mathrm{B}$ is benthic reflectance, $\theta_\mathrm{v}$ is viewing
zenith angle, $\Delta\phi_\mathrm{v}$ is sensor viewing azimuth angle
relative to the solar plane, and $\theta_\mathrm{s}$ is the solar zenith
angle.

The forward model has two options for paramterization to simulate
$R_\mathrm{rs}$, i.e. Albert & Mobley (2003) (AM03 onwards) and Lee et
al. (1998) (L98 onwards) for both optically deep and shallow water
types. However, currently, AM03 is implemented with a much greater
support of bio-optical paramterizations for deriving spectral slopes and
calculation of inelastic scattering (SICF and $f_\mathrm{DOM}$). Thus,
it is recommended to users to use AM03 to perform forward simulations of
SRTE. The AM03 forward model for elastic scattering part can be
expressed as below:

$$R_{r s}^{e, m}(\lambda)=R_{r s, \infty} \cdot[1-\exp {-(K_d+\frac{K_{u, W}}{\cos \theta_v}) H}] + (\sum_{i=1}^n f_{a, i} B_i \cdot R_B) \cdot \exp {-(K_d+\frac{K_{u, B}}{\cos \theta_v}) H}$$

-   Where, $R_{r s, \infty}$ is the contribution from water column only
    (in deep waters, that is the only feasible component) which can be
    expressed as below:

$$R_{r s, \infty}=f^{\uparrow}\left(\omega_b, \theta_s, u, \theta_v\right) \cdot \omega_b=$$

$$f^{\uparrow}\left(\omega_b\right) \cdot f^{\uparrow}\left(\theta_s\right) \cdot f^{\uparrow}(u) \cdot f^{\uparrow}\left(\theta_v\right) \cdot \omega_b= $$

$$f^{\uparrow}\left(\omega_b\right) \cdot f^{\uparrow}\left(\theta_s\right) \cdot f^{\uparrow}(u) \cdot f^{\uparrow}\left(\theta_v\right) \cdot \omega_b= \\$$

$$P_{r s, 1} \cdot\left(1+P_{r s, 2} \cdot \omega_b+P_{r s, 3} \cdot \omega_b^2+P_{r s, 4} \cdot \omega_b^3\right) \cdot\left(1+P_{r s, 5} \cdot \frac{1}{\cos \theta_s}\right) \cdot(1 \\
\left.+P_{r s, 6} \cdot u\right) \cdot\left(1+P_{r s, 7} \cdot \frac{1}{\cos \theta_v}\right) \cdot \omega_b$$

-   Where, $u$ is wind speed at water surface, $\omega_\mathrm{b}$ is
    the single-scattering albedo calculated as
    $\frac{b_\mathrm{b}}{a +b_\mathrm{b}}$ and $P_{r s, 1-7}$ are
    parametric coefficients.

-   $K_\mathrm{d}$ and $K_\mathrm{u,W}$ and $K_\mathrm{u,B}$ are
    downwelling diffused attenuation coefficient and as diffused
    attenuation coefficient for upwelling light from water column and
    bottom respectively. These are paramterized as below:

$$K_d=\kappa_0 \cdot \frac{a+b_b}{\cos \theta_s}, K_u=\left(a+b_b\right) \cdot\left(1+\omega_b\right)^{\kappa_i} \cdot\left(1+\kappa_2 \frac{1}{\cos \theta_s}\right)$$

-   Where, $\kappa_0$, $\kappa_i$ and $\kappa_2$ are parametric
    coefficients. Among rest, $f_\mathrm{a,i}$ is the areal fraction of
    $R_\mathrm{B}$ constrained up to five bottom types and
    $B_\mathrm{i}$ is the Bidirectional Reflectance Distribution
    Function (BRDF) for each bottom type assumed to be lambertian
    ($\pi$).

The $a$ and $b_\mathrm{b}$ values can either be *in situ* IOPs
(processed following the IOCCG/NASA protocols) directly or in absence of
*in situ* IOPs, *in vivo* values of $[chl]$, $a_\mathrm{g}(443)$ and
$a_\mathrm{d}(443)$ can be supplied to the forward model. In that case,
the subsequent IOPs are calculated as following bio-optical
parametrizations:

-   For phytoplankton absorption,
    $$a_\mathrm{\phi}(\lambda) = (a_\mathrm{0}(\lambda) + a_\mathrm{1}(\lambda) \ln[a_\mathrm{\phi}(443)]) \cdot a_\mathrm{\phi}(443)$$
    Where, $a_\mathrm{\phi}(443) = 0.06[chl]^{0.65}$, $a_\mathrm{0}$ and
    $a_\mathrm{1}$ are precomputed parameters from Lee et al. (1994).

-   For coloured detrital matter absorption,
    $$a_\mathrm{dg}(\lambda) = a_\mathrm{dg}(443)\exp {(-S_\mathrm{dg}(\lambda - 443))}$$
    Where, $-S_\mathrm{dg}$ is the spectral slope which can be retrieved
    using QAA in deep waters as following:

$S_{\mathrm{dg}}=0.015+\left(\frac{0.002}{0.6+\frac{r_{\mathrm{rs}}(443)}{r_{r s}(555)}}\right)$
The $r_\mathrm{rs}$ is sub-surface Remote-sensing reflectance which can
be calculated as
$r_\mathrm{rs} = \frac{R_\mathrm{rs}}{0.52 + 1.7 \cdot R_\mathrm{rs}}$.
In case of shallow water, a default value of 0.017 is used which can be
changed by the user as well (For customization See Section
"R-implementation").

-   For particulate backscatter,
    $$b_\mathrm{bp}(\lambda) = b_\mathrm{bp}(555){\left (\frac {\lambda}{555} \right)}^\eta$$
    Where, $\eta$ is the spectral slope which can be retrieved using QAA
    in deep waters as following:

$\eta = 2.0 \left(1 - 1.2 \exp { \left( -0.9 \frac {r_\mathrm{rs}(443)} {r_\mathrm{rs}(555)}\right)}  \right)$.
In case of shallow water, a default value of -0.46 is used which can be
changed by the user as well (For customization See Section
"R-implementation").

The forward model also has support for simulation of inelastic
scattering in water, i.e. Sun Induced Chlorophyll Fluoroscence (SICF)
and CDOM Fluoroscence ($f_\mathrm{DOM}$). The SICF is calculated
semi-analytically following the study by Gilerson et al. (2007) whereas
the $f_\mathrm{DOM}$ is calculated analytically following Mobley (1994).

-   SICF affects the wavelengths between 670-710 nm $\lambda$ with a
    peak at $\sim 685$ nm. The SICF equivalent radiance at the peak SICF
    emission wavelength, 685nm, i.e. $L_\mathrm{\phi}(685)$ is
    calculated as,
    $$L_\mathrm{\phi}(685) = \frac {0.092[chl]} {(1+0.40 a_\mathrm{dg}(443) + 0.078[chl])}$$.
    The spectral shape of SICF is modelled as summation of two Gaussian
    shape with $\mu = 685,730$ and $\sigma = 25,50$ as the Full Wave
    Half Maximum (FWHM) of the Gaussian curve:

$$ hc(\lambda) = \sqrt{\frac{4 \ln 2}{\pi}} \frac{1}{25} \exp \left[-4 \ln 2\left(\frac{\lambda-685}{25}\right)^2\right]+ \sqrt{\frac{4 \ln 2}{\pi}} \frac{1}{50} \exp \left[-4 \ln 2\left(\frac{\lambda-730}{50}\right)^2\right] $$

The SICF equivalent radiance for $670 < \lambda < 710$ can be calculated
as, $L_\mathrm{SICF} = L_\mathrm{\phi}(685) \cdot hc(\lambda)$. Finally
the SICF equivalent $r_\mathrm{rs}$ is calculated as,
$r_\mathrm{rs,SICF}=\frac {L_\mathrm{SICF}} {E_\mathrm{d}(0^-)}$, where,
$E_\mathrm{d}(0^-)$ is the sub-surface downwelling irradiance.The
$E_\mathrm{d}(0^-)$ can either be calculated analytically following
Gregg & Carder (1990) or from a user supplied file for Spectral water
surface irradiance.

-   $f_\mathrm{DOM}$ affects $\lambda$ between 310-600nm and the
    spectral shape can be modelled following a log-normal distribution.
    As the $f_\mathrm{DOM}$ equivalent radiance, $L_\mathrm{fDOM}$ is
    calculated analytically by integrating a wavelength discretized
    source scatter function as

$$\beta_\mathrm{g} (\psi,\lambda^′,\lambda) =   b_\mathrm{g}(\lambda^′)f_\mathrm{g}(\lambda^′, \lambda)\beta_\mathrm{g} (\psi)$$

Where, $b_\mathrm{g}(\lambda^′)$ is loss of photon at excitation
wavelength, $\lambda^′$ which is equivalent to
$a_\mathrm{CDOM}(\lambda^′)$. $f_\mathrm{g}(\lambda^′, \lambda)$ is the
wavelength redistribution function (WRF) calculated as,
$f_\mathrm{g}(\lambda^′, \lambda) = \eta_\mathrm{g}(\lambda^′, \lambda) \frac {\lambda^′} {\lambda}$,
where, $\eta_\mathrm{g}(\lambda^′, \lambda)$ is the shape of quantum
efficiency of $f_\mathrm{DOM}$ modelled after Howes et al. (1992) as

$$\eta_\mathrm{g}(\lambda^′, \lambda) = A_0\left(\lambda^{\prime}\right) \exp \left[-\left(\frac{\frac{1}{\lambda}-\frac{A_1}{\lambda^{\prime}}-B_1}{0.6\left(\frac{A_2}{\lambda^{\prime}}+B_2\right)}\right)^2\right] $$

Where, $A_0(\lambda^{\prime}), A_1, A_2, B_1, B_2$ are coefficients
given by Howes et al. (1992). The phase function is assumed isotropic
thus, $\beta_\mathrm{g} (\psi) = \frac {1} {4\pi}$ is set. Finally, the
$f_\mathrm{DOM}$ equivalent $r_\mathrm{rs}$ is calculated as,
$r_\mathrm{rs,fDOM}=\frac {L_\mathrm{fDOM}} {E_\mathrm{d}(0^-)}$. The
$E_\mathrm{d}(0^-)$ is calculated in a similar way to the same for
$r_\mathrm{rs,SICF}$.

-   INVERSE MODEL

The inverse model is used to retrieve sub-surface
$[chl], a_\mathrm{dg}(443)$ and $H, R_\mathrm{B}$ in case of a shallow
water column from *in situ* or Satelite observations of $R_\mathrm{rs}$.
The values are obtained by an optimization routine that tries to find
the global mimima for the joint distribution of the paramteric space,
$\phi_\mathrm{par} : \{[chl], a_\mathrm{dg}(443), H, R_\mathrm{B}\}^{\mathbb{R}^{\lambda_\mathrm{i} \times D}}$.

In such problems, a cost/objective function, $s_\mathrm{obj}$ is
formulated for $\phi_\mathrm{par}$ which is then minimized in an
optimization problem. S.A.B.E.R. has two built-in $s_\mathrm{obj}$,
i.e. 1. sum-square-of-residuals similar to OLS predictor and 2. sum of
log-likelihood of the forward model where the model error is gaussian in
nature with a known mean $\mu=0$ but standard deviation $\sigma$ is
unknown. For the optimization routines, there are two set of families of
optimizer, i.e. 1. gradient based optimization
(Box-constraint/Levenberg-Marquardt) and 2. Bayesian Markov Chain Monte
Carlo (MCMC) based sampling. The first type calculates optimal set of
parameters by through analytical calculation of the jacobian matrix with
respect to the $s_\mathrm{obj}$ and trying to find the point of
convergence where the change in hessian matrix becomes flat, i.e.
$\frac{\partial{}^2 s_\mathrm{obj}}{\partial{}{[chl]}^2} + \frac{\partial{}^2 s_\mathrm{obj}}{\partial{}{[a_\mathrm{dg}(443)]}^2} + \frac{\partial{}^2 s_\mathrm{obj}}{\partial{}{[H]}^2} + \frac{\partial{}^2 s_\mathrm{obj}}{\partial{}{[R_\mathrm{B}]}^2} = 0$.
However, as the inversion of SRTE in optically complex waters is
*ill-posed* problem, gradient-based solutions may tend to retrieve
mathematically ambiguous sets of parameters, moreover, such solvers are
sensitive to initial values, it can take a long time to reach the global
minima or sometimes, it may fail to reach the global minima. To mitigate
such problems in inversion of SRTE, an additional probabilistic sampling
scheme, MCMC has been integrated in S.A.B.E.R.

The MCMC framework randomly samples parametric sets following a
probabilistic candidate function and compare the quadrature in
$s_\mathrm{obj}$. The convergence here in MCMC is diagnosed by reaching
*stationary distribution* for the parametric space. The posterior
distribution, which is the joint probability distribution of parameter
space, is calculated as:
$$\mathbb{P}(par|R_\mathrm{rs}) \propto \mathbb{P}(par) \mathbb{P}(R_\mathrm{rs}|par)$$

Where $\mathbb{P}(par)$ is the prior probability of parameters.
S.A.B.E.R. is developed with *Weibull* family of distribution for
$[chl], a_\mathrm{dg}(443)$ and $H$ following the assumption that the
values can never be negative in reality. The density function is given
as:
$\{ f(x; \lambda, k) = \frac {k} {\lambda} \left(\frac {x} {\lambda}\right)^{k-1}\exp{(\frac {-x} {\lambda})^k} \} \exists \{ x \geq 0 \}$
where,$\lambda,k$ are the scale and shape parameter of the density
function. For shallow water parameter $R_\mathrm{B}$, which is expressed
as a linear mixture of fractions of possible bottom class, is assumed to
be *Uniform* distribution in nature varying between 0-1. The density
function is given as:
$\{f(x)=\frac {1}{b-a} \} \exists \{ a \leq x \geq \forall a=0, b=1 \}$
where $\mathbb{P}(R_\mathrm{rs}|par)$ is the log-likelihood of the
forward model which actually infers the model noise which is found to be
Gaussian, thus the log-likelihood can be expressed as:

$$
\log \left(\mathcal{L}\left(\mu, \sigma^2\right)\right)=-\frac{n}{2} \log \left(2 \pi \sigma^2\right)-\frac{1}{2 \sigma^2} \sum_{i-1}^n\left(x_i-\mu\right)^2
$$

Where, $x_\mathrm{i}$ is the observed $R_\mathrm{rs}$;
$R_\mathrm{rs, obs}$ that is input to the inverse problem and $\mu$ is
the forward SA.B.E.R. simulated $R_\mathrm{rs}$, The population standard
deviation $\sigma$ is usually unknown and also retrieved as a parameter
from the inversion.

The MCMC sampling can be performed in a variety of stochastic
approaches, S.A.B.E.R. is default with the "Delayed Rejection
Metropolis-Hastings" algorithm and a parallel computation of
log-likelihood values for a large number, i.e. 10,000 number of
iterations. Initial \~2500 iterations are discarded as *"burn-in"*
period prior to ascension to stationarity. Post reaching stationary
distribution, *"Maximum-a-Posterior"* (MAP) estimates;
$\theta_\mathrm{MAP} = \underset{par}{\operatorname{argmax}} \mathbb{P}(par|R_\mathrm{rs})$
are obtained as optimal model solutions of $par$, i.e.
$[chl], a_\mathrm{dg}(443)$, $H$ and $R_\mathrm{B}$ fractions. The
parametric uncertainty is also retrieved from MCMC chains as *"Credible
intervals"*.

Apart from full parametric inversion as described below, an additional
retrieval of spectral bottom reflectance;$R_\mathrm{B}(\lambda)$ using
algebraic simplification of the AM03 forward model has been included in
S.A.B.E.R. However, it is only to be used when the rest of the
parameters, i.e. BGC parameters such as $[chl], a_\mathrm{dg}(443)$ or
the IOPs and $H$ are known. Moreover, the $R_\mathrm{rs,SICF}$ and
$R_\mathrm{rs,fDOM}$ components should be subtracted from the
$R_\mathrm{rs, obs}$ prior to input in the following model where
$R_{r s}^{obs,e}$ is inelastic component contribution excluded
$R_\mathrm{rs, obs}$:

$$ R_B = \frac{R_{r s}^{obs,e}(\lambda) - R_{r s, \infty} \cdot [1 - \exp{(-(K_d + \frac{K_{u,W}}{\cos \theta_v})H)}]}{\exp{(-(K_d + \frac{K_{u,B}}{\cos \theta_v})H)}} $$


# REFERENCES

Albert, A., & Mobley, C. (2003). An analytical model for subsurface
irradiance and remote sensing reflectance in deep and shallow case-2
waters. Optics Express, 11(22), 2873.
<https://doi.org/10.1364/oe.11.002873>

A. Gilerson, J. Zhou, S. Hlaing, I. Ioannou, J. Schalles, B. Gross, F.
Moshary, and S. Ahmed, "Fluorescence component in the reflectance
spectra from coastal waters. Dependence on water composition," Opt.
Express 15, 15702-15721 (2007)

Gregg, Watson W., Carder, K. L., (1990), A simple spectral solar
irradiance model for cloudless maritime atmospheres, Limnology and
Oceanography, 35, doi: 10.4319/lo.1990.35.8.1657.

Lee, Z., Carder, K. L., Hawes, S. K., Steward, R. G., Peacock, T. G., &
Davis, C. O. (1994). Model for the interpretation of hyperspectral
remote-sensing reflectance. Applied Optics, 33(24), 5721.
<https://doi.org/10.1364/AO.33.005721>

Lee, Z., Carder, K. L., Mobley, C. D., Steward, R. G., & Patch, J. S.
(1998). Hyperspectral remote sensing for shallow waters I A
semianalytical model. Applied Optics, 37(27), 6329.
<https://doi.org/10.1364/ao.37.006329>

Lee, ZhongPing, Carder, K. L., & Arnone, R. A. (2002). Deriving inherent
optical properties from water color: a multiband quasi-analytical
algorithm for optically deep waters. Applied Optics, 41(27), 5755.
<https://doi.org/10.1364/ao.41.005755>

Mobley, C. D. (1994). Radiative Transfer: Across the surface. Light and
Water: Radiative Transfer in the Natural Waters, 159–161.
