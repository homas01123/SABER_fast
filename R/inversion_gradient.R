#' Perform inversion based on various optimization methods.
#' Allow for unconstrained inversion of chl [mg/m^-3], ag_440 [m^-1],
#' bbp_550 [m^-1], and in optically shallow waters h_w [m],
#' and bottom reflectance fraction. Optionally, can perform constrained
#' inversion by providing parameters in `fixed_par` instead of `par_init`.
#'
#' @author Soham Mukherjee, Raphael Mabit
#'
#' @param forward_model c("am03", "lee98")
#' @param objective_fct c("log-ll", "SSR", "lee99") SSR is not implemented in the code
#' @param optim_mtd c("nelder-mead", "BFGS", "CG", "L-BFGS-B", "sann", "brent","levenberg-marqardt", "auglag")
#' @param par_init a tibble with the initial parameter to be inverted.
#'  For best results, can be estimated with `pre_fit_inversion`.
#'  \describe{
#'    \item{chl}{chlorophyl-a concentration in [mg/m^3]}
#'    \item{ag_440}{CDOM absorption [m^-1] at 440 nm}
#'    \item{bbp_550}{particulate backscattering [m^-1] at 550 nm}
#'    \item{h_w}{watercolumn height above the bottom [m]}
#'    \item{rb_*}{fraction of end-member bottom reflectance class}
#'    \item{sd}{Optional, standard deviaton of the population for `objective_fct = "log-ll"`}
#'  }
#' @param fixed_par Optional, a tibble with the same columns as par_init.
#'  Parameter defined here will be considered known, hence not retrieved during
#'  optimization. If defined here they must not be defined in `par_init`.
#' @param lower_b Optional, lower boundary for `optim_mtd = "L-BFGS-B"`.
#' If not provided will be calculated in `parse_inverse_parameter`
#' @param upper_b Optional, same as `lower_b`.
#' @param verbose guess
#'
#' @export
inversion_gradient <- function(
    rrs,
    forward_model,
    objective_fct,
    optim_mtd,
    par_inversed,
    par_fixed = NULL,
    lower_b = NULL,
    init_val = NULL,
    upper_b = NULL,
    verbose = F) {
  rlang::inform(paste0("\033[0;33m", "###################################################################", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;39m", "########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;32m", "###################################################################", "\033[0m", "\n"))

  minimization_fct <- objective_factory(
    model = forward_model,
    objective = objective_fct,
    rrs_observed = rrs,
    par_fixed = par_fixed,
    par_inversed = par_inversed
  )

  # Instantiate initial values
  params <- parse_inverse_parameter(
    par_df = data.frame("name" = par_inversed, "value" = init_val),
    optim_mtd=optim_mtd,
    lower_b = lower_b,
    upper_b = upper_b,
    verbose = verbose
  )

  par <- params$par
  names(par) <- params$names
  lower_b <- params$lower
  upper_b <- params$upper
  parscale <- params$parscale

  # Optimization ------------------------------------------------------------

  start.time <- Sys.time()

  if (optim_mtd == "L-BFGS-B") {
    optim_result <- optim(
      par = par,
      fn = minimization_fct,
      method = optim_mtd,
      lower = lower_b,
      upper = upper_b,
      control = list(parscale = parscale)
    )
  }

  if (optim_mtd == "Nelder-Mead" |
    optim_mtd == "SANN" |
    optim_mtd == "Brent") {
    optim_result <- optim(
      par = par,
      fn = minimization_fct,
      method = optim_mtd,
      control = list(parscale = parscale),
      hessian = FALSE
    )
  }

  if (optim_mtd == "levenberg-marqardt") {
    # TODO: Does not support bounds, will break when chl becomes negative !

    lm_result <- marqLevAlg::marqLevAlg(
      b = par,
      fn = minimization_fct,
      print.info = F
    )

    optim_result <- tibble("par" = lm_result$b)
  }

  if (optim_mtd == "auglag") {
    print("Augmented Lagriangian with equality constraints will be used for inversion")

    get_fraction_indices <- function(par_names, prefix = "rb_") {
      which(grepl(paste0("^", prefix), par_names))
    }

    # Identify indices dynamically
    fraction_indices <- get_fraction_indices(names(par))

    # bounds for areal fractions
    fheq <- function(pars) {
      sum(pars[fraction_indices]) - 1 # must sum to 1
    }

    fhin <- function(pars) {
      pars[fraction_indices] # All should be ≥ 0
    }

    optim_result <- alabama::auglag(
      fn = minimization_fct,
      par = par,
      heq = fheq,
      hin = fhin,
      control.outer = list(trace = F, method = "nlminb")
    )
  }

  # Calculate uncertainty ---------------------------------------------------

  # Calculate hessian matrix for var-covar matrix
  if (optim_mtd == "auglag") {
    hessian_inverse <- optim_result$hessian
  } else {
    hessian_inverse <- numDeriv::hessian(
      x = optim_result$par,
      func = minimization_fct # ,
      # data=obsdata
    )
  }

  if (verbose) {
    rownames(hessian_inverse) <- par_inversed
    colnames(hessian_inverse) <- par_inversed
    rlang::inform(paste0("\033[0;32m", "#################### VAR-COV HESSIAN MATRIX #########################", "\033[0m", "\n"))
    prmatrix(hessian_inverse)
  }

  param_estimate <- optim_result$par

  param_sd <- tryCatch({
    sqrt(diag(solve(hessian.inverse)))}, #solve for diagonal elements to get sd
    error = myFun
  )

  # param_sd <- tryCatch(
  #   {
  #     sqrt(diag(solve(hessian_inverse)))
  #   }, # solve for diagonal elements to get sd
  #   objective = NA
  # )

  end.time <- Sys.time()

  if (!is.numeric(param_sd)) {
    rlang::warn(
      paste0("\033[0;31m", "Failed to calculate diagonal of hessian from
             high degree of correlation, coerce to NA", "\033[0m", "\n")
    )
  }

  # Maximum Likelihood Estimates
  mle <- tibble(
    "name" = par_inversed,
    "estimate" = param_estimate,
    "sd" = param_sd
  )

  if (verbose) {
    if (optim_result$convergence == 0) {
      # convergence <- "TRUE"
      rlang::inform(paste0("\033[0;32m", "CONVERGENCE: GLOBAL", "\033[0m", "\n"))
    } else {
      # convergence = "FALSE"
      rlang::inform(paste0("\033[0;34m", "CONVERGENCE: LOCAL", "\033[0m", "\n"))
    }

    time_taken <- end.time - start.time
    rlang::inform(glue::glue("time.elapsed: ", time_taken))
    # return(list(mle, "convergence"= convergence))
  }

  return(mle)
}


#' @export
myFun <- function(x) {
  NA
}

#' pre_fit_inversion
#'

# pre_fit_inversion <- function() {
#   # 5.1 Pre-FIT of initial values
#   if (preFit == TRUE) {
#     pre.Fit <- data.frame(
#       "C_ph" = seq(1, 10, 0.5), # <<USER DEFINED >>
#       "a_cdom.440" = seq(0.5, 5, 0.25),
#       "a.nap.440" = seq(0.01, 0.1, 0.005)
#     )
#
#     pre.Fit.input.LUT <- expand.grid(pre.Fit) # Create pre-Fit parameter space LUT
#
#     preFIT.rrs.forward.LUT <- matrix(
#       nrow = length(pre.Fit.input.LUT$C_ph),
#       ncol = length(wavelength), 0
#     )
#     # Create the Progress Bar
#     pb <- txtProgressBar(
#       min = 0, # Minimum value of the progress bar
#       max = length(pre.Fit.input.LUT$C_ph), # Maximum value of the progress bar
#       style = 3, # Progress bar style (also available style = 1 and style = 2)
#       width = 50, # Progress bar width. Defaults to getOption("width")
#       char = "="
#     )
#
#     reslist <- vector()
#     for (i in 1:length(pre.Fit.input.LUT$C_ph)) { # Create Rrs LUT
#       temp1 <- as.numeric(pre.Fit.input.LUT[i, ])
#       temp2 <- Saber_forward(
#         chl = temp1[1], acdom440 = temp1[2],
#         anap440 = temp1[3], bbp.550 = Fit.input$bbp.550,
#         realdata = obsdata, verbose = F
#       )
#
#       preFIT.rrs.forward.LUT[i, ] <- temp2[[1]]$Rrs
#       reslist[i] <- temp2[[2]]
#       # cat(paste0("\033[0;43m",i," iterations over, ", (nrow(preFIT.rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
#       setTxtProgressBar(pb, i)
#       if (i == length(pre.Fit.input.LUT$C_ph)) {
#         cat(paste0("\033[0;32m", "###############PRE-FIT FINISHED################", "\033[0m", "\n"))
#       }
#     }
#
#     prefit.best <- pre.Fit.input.LUT[which.min(reslist), ] # retrieve best initial values
#     # using C.R.I.S.T.A.L.[minimizing SSR]
#     prefit.best
#
#     rrs.prefit <- Saber_forward(
#       chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440,
#       anap440 = prefit.best$a.nap.440, bbp.550 = Fit.input$bbp.550,
#       realdata = obsdata, verbose = T
#     )[[1]]$Rrs
#
#     # Show prefit spectra (Convert to ggplot2)
#     plot(wavelength, obsdata, type = "l", col = "red", ylim = c(0, max(obsdata)))
#     lines(wavelength, rrs.prefit, col = "green")
#   }
# }
