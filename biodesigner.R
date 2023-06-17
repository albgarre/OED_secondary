
library(tidyverse)
library(FME)
library(MEIGOR)

## Model equation

CPM_model <- function(par, x) {

    num <- (x-par$Xmax)*(x-par$Xmin)^par$n
    den <- (par$Xopt-par$Xmin)^(par$n-1)*( (par$Xopt-par$Xmin)*(x-par$Xopt) - (par$Xopt-par$Xmax)*((par$n-1)*par$Xopt + par$Xmin-par$n*x) )
    gamma <- num/den
    gamma[x < par$Xmin] <- 0
    gamma[x > par$Xmax] <- 0
    
    y <- gamma*par$mu_opt
    
    data.frame(x = x, y = y, sq_y = sqrt(y))
}

Ratkowsky <- function(par, x) {
    
    sq_mu <- par$b*(x - par$Xmin)*(1 - exp(par$c*(x - par$Xmax)))
    
    sq_mu[x < par$Xmin] <- 0
    sq_mu[x > par$Xmax] <- 0
    
    data.frame(x = x, y = sq_mu^2, sq_y = sq_mu)
}

## FIM functions

get_sensitivities <- function(model = "CPM", pars, 
                              tol_design = 1,
                              include = NULL,
                              varscale = 1,
                              parscale = 1,
                              output = "sq_y") {
    
    tgt <- switch(model,
           CPM = CPM_model,
           Ratkowsky = Ratkowsky)
    
    out <- sensFun(tgt, pars, 
            x = seq(pars$Xmin + tol_design, 
                    pars$Xmax - tol_design, 
                    length = 500),
            varscale = 1,
            parscale = 1)  %>% 
        filter(var == output) %>%
        select(-var)
    
    if (!is.null(include)) {
        out <- select(out, x, matches(include))
    }
    
    out
    
    
    
}

calculate_FIM <- function(sensitivities, times) {
    
    old_times <- sensitivities$x
    
    sensitivities <- select(sensitivities, -x) # %>%
    # as.matrix()
    
    interp_sensitivities <- matrix(NA, length(times), ncol(sensitivities))
    
    for (i in 1:ncol(sensitivities)) {
        interp_sensitivities[, i] <- approx(old_times, sensitivities[,i],
                                            times)$y
    }
    
    FIM <- t(interp_sensitivities) %*% interp_sensitivities
    FIM
}

objective_D <- function(times, sensitivities, ...) {
    
    # sensitivities is taken from the global environment
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- -det(FIM)
    
    out
    
}

objective_Emod <- function(times, sensitivities, ...) {
    
    # sensitivities is taken from the global environment
    
    FIM <- calculate_FIM(sensitivities, times)
    
    eig_vals <- eigen(FIM, only.values=TRUE)
    
    min_eig <- min(eig_vals$values)
    
    if (abs(min_eig) < 1e-10) {
        
        return(1e20)
        
    } else {
        
        abs(max(eig_vals$values)/min_eig)
        
    }
    
}

objective_penalty <- function(times, sensitivities, psi, model, pars) {
    
    # sensitivities, psi and objective are taken from the global environment
    
    n_increments <- 6
    
    ## FIM part
    
    FIM_based <- objective_Emod(times, sensitivities)
    
    ## Penalty part
    
    tgt <- switch (model,
        CPM = CPM_model,
        Ratkowsky = Ratkowsky
    )
    
    mus <- tgt(pars, times)$y
    exp_time <- 6/mus
    penalty <- sum(exp_time)*psi  # Works for E-mod only
    
    ## Calculate the objective
    
    FIM_based + penalty
    
}

calculate_OED <- function(model = "CPM", pars,
                          n_points,
                          tol_design = 1,
                          include = NULL,
                          maxeval = 3000,
                          criterion = "D",
                          psi = 1,
                          output = "sq_y") {
    
    ## Get the sensitivities
    
    sens <- get_sensitivities(model = model, pars = pars, 
                              tol_design = tol_design,
                              include = include,
                              varscale = 1,
                              parscale = 1,
                              output = output
    )
    
    ## Options of the optimizer

    opts_global <- list(maxeval=maxeval,  local_solver=0,
                        local_finish="DHC", local_iterprint=1)
    
    ## OED - D criterion
    
    tgt <- switch(criterion,
                  D = objective_D,
                  E = objective_Emod,
                  penalty = objective_penalty)
    
    problem <- list(f=tgt,
                    x_L=rep(pars$Xmin + tol_design, n_points),
                    x_U=rep(pars$Xmax - tol_design, n_points)
    )
    
    MEIGO(problem, opts_global, algorithm="ESS", sensitivities = sens, 
          psi = psi, model = model, pars = pars)
    
}

## Model fitting functions

get_residuals <- function(p, this_data, known, model) {
    
    p <- as.list(c(p, known))
    
    tgt <- switch(model,
                  Ratkowsky = Ratkowsky,
                  CPM = CPM_model)
    
    pred <- tgt(p, this_data$x)$sq_y
    
    this_data$sq_mu_obs - pred
    
}

## Aggregated time

get_total_time <- function(model, pars, times, n_increments = 6) {
    
    tgt <- switch(model,
                  Ratkowsky = Ratkowsky,
                  CPM = CPM_model)
    
    each_time <- n_increments/tgt(pars, times)$y
    sum(each_time)
    
}

## Simulating experiments

simulate_experiments <- function(model, pars_fit, pars_known,
                                 sampled_x, n_repetitions,
                                 sigma_sq_mu, n_exp,
                                 lower = -Inf,
                                 upper = Inf) {
    
    ## Ideal response
    
    tgt <- switch (model,
        CPM = CPM_model,
        Ratkowsky = Ratkowsky
    )
    
    ideal_curve <- tibble(x = sort(sampled_x),
                          sq_mu = tgt(c(pars_fit, pars_known), x)$sq_y
                          ) %>%
        mutate(
            sq_mu = ifelse(sq_mu < 0, 0, sq_mu)
        )
    
    ideal_curve <- bind_rows(replicate(n_repetitions, ideal_curve, simplify = FALSE))
    
    
    ## Simulating experiment ---------------------------------------------------
    
    ## Observed ones
    
    exp <- c(1:n_exp) %>%
        map(.,
            ~ mutate(ideal_curve,
                     sq_mu_obs = sq_mu + rnorm(nrow(ideal_curve), mean = 0, sd = sigma_sq_mu)
            )
        )
    
    ## Fitted models
    
    p0 <- unlist(pars_fit)
    
    # get_residuals(p0, exp_penalty[[1]], known = c(n = 2))
    
    my_models <- exp %>%
        map(.,
            ~ modFit(get_residuals, p0, this_data = ., 
                     known = unlist(pars_known),
                     lower = lower, upper = upper,
                     model = model
                     )
            )
    
    my_models
    
}

## -----------------------------------------------------------------------------

# do_everything <- function(par_values, n_points, tol_design = 1,
#                           n_exp, sigma_sq_mu,
#                           remove_n = FALSE, psi = 1) {
#     
#     ## Sensitivity functions
#     
#     my_sens <- sensFun(Ratkowsky, initial_guess, 
#                        x = seq(guess_min + tol_design, 
#                                guess_max - tol_design, 
#                                length = 200),
#                        varscale = 1,
#                        parscale = 1)  %>% 
#         filter(var == "sq_y") 
#     
#     if (remove_n) {
#         my_sens <- my_sens %>% select(-n)
#     }
#     
#     ## Options of the optimizer
#     
#     max_x <- guess_max - tol_design
#     min_x <- guess_min + tol_design
#     
#     opts_global <- list(maxeval=5000,  local_solver=0,
#                         local_finish="DHC", local_iterprint=1)
#     
#     ## OED - D criterion
#     
#     problem <- list(f=objective_D,
#                     x_L=rep(min_x, n_points),
#                     x_U=rep(max_x, n_points)
#     )
#     
#     OED_Dcrit <- MEIGO(problem, opts_global, algorithm="ESS",
#                        sensitivities = my_sens)
#     
#     ## OED - E criterion
#     
#     problem <- list(f=objective_Emod,
#                     x_L=rep(min_x, n_points),
#                     x_U=rep(max_x, n_points)
#     )
#     
#     OED_Emod <- MEIGO(problem, opts_global, algorithm="ESS",
#                       sensitivities = my_sens)
#     
    # ## Ideal response
    # 
    # ideal_OED <- tibble(x = sort(OED_Emod$xbest),
    #                     sq_mu = Ratkowsky(initial_guess, x)$sq_y
    # ) %>%
    #     bind_rows(., .)  # 2 repetitions
    # 
    # ideal_OED_Dcrit <- tibble(x = sort(OED_Dcrit$xbest),
    #                           sq_mu = Ratkowsky(initial_guess, x)$sq_y
    # ) %>%
    #     bind_rows(., .)  # 2 repetitions
    # 
    # # ideal_penalty <- tibble(x = sort(OED_penalty$xbest),
    # #                         sq_mu = Ratkowsky(initial_guess, x)$sq_y
    # # ) %>%
    # #     bind_rows(., .)  # 2 repetitions
    # 
    # ideal_uniform <- tibble(x = seq(initial_guess$Xmin + tol_design,
    #                                 initial_guess$Xmax - tol_design,
    #                                 length = n_points),
    #                         sq_mu = Ratkowsky(initial_guess, x)$sq_y
    # ) %>%
    #     bind_rows(., .)  # 2 repetitions
    # 
    # ## Simulating experiment ---------------------------------------------------
    # 
    # ## Observed ones
    # 
    # exp_uniform <- c(1:n_exp) %>%
    #     map(.,
    #         ~ mutate(ideal_uniform,
    #                  sq_mu_obs = sq_mu + rnorm(nrow(ideal_uniform), mean = 0, sd = sigma_sq_mu)
    #         )
    #     )
    # 
    # exp_OED <- c(1:n_exp) %>%
    #     map(.,
    #         ~ mutate(ideal_OED,
    #                  sq_mu_obs = sq_mu + rnorm(nrow(ideal_OED), mean = 0, sd = sigma_sq_mu)
    #         )
    #     )
    # 
    # exp_OED_Dcrit <- c(1:n_exp) %>%
    #     map(.,
    #         ~ mutate(ideal_OED_Dcrit,
    #                  sq_mu_obs = sq_mu + rnorm(nrow(ideal_OED_Dcrit), mean = 0, sd = sigma_sq_mu)
    #         )
    #     )
    # 
    # # exp_penalty <- c(1:n_exp) %>%
    # #     map(.,
    # #         ~ mutate(ideal_penalty,
    # #                  sq_mu_obs = sq_mu + rnorm(nrow(ideal_penalty), mean = 0, sd = sigma_sq_mu)
    # #         )
    # #     )
    # 
    # ## Fitted models
    # 
    # p0 <- unlist(initial_guess)
    # 
    # # get_residuals(p0, exp_penalty[[1]], known = c(n = 2))
    # 
    # models_uniform <- exp_uniform %>%
    #     map(.,
    #         ~ modFit(get_residuals, p0, this_data = ., known = c(n = 2),
    #                  lower = c(Xmin = -10, Xmax = 35, b = 0, c = 0),
    #                  upper = c(Xmin = 10, Xmax = 55, b = .25, c = 5)
    #         )
    #     )
    # 
    # models_OED <- exp_OED %>%
    #     map(.,
    #         ~ modFit(get_residuals, p0, this_data = ., known = c(n = 2),
    #                  lower = c(Xmin = -10, Xmax = 35, b = 0, c = 0),
    #                  upper = c(Xmin = 10, Xmax = 55, b = .25, c = 5)
    #         )
    #     )
    # 
    # 
    # models_OED_Dcrit <- exp_OED_Dcrit %>%
    #     map(.,
    #         ~ modFit(get_residuals, p0, this_data = ., known = c(n = 2),
    #                  lower = c(Xmin = -10, Xmax = 35, b = 0, c = 0),
    #                  upper = c(Xmin = 10, Xmax = 55, b = .25, c = 5)
    #         )
    #     )
#     
#     # models_penalty <- exp_penalty %>%
#     #     map(.,
#     #         ~ modFit(get_residuals, p0, this_data = ., known = c(n = 2),
#     #                  lower = c(Xmin = -10, Xmax = 35, b = 0, c = 0),
#     #                  upper = c(Xmin = 10, Xmax = 55, b = .25, c = 5)
#     #         )
#     #     ) 
#     
#     
#     
# }

## 

# guess_min <- 5
# guess_max <- 45
# guess_b <- .035
# guess_c <- .5
# 
# 
# initial_guess <- list( 
#     Xmin = guess_min, 
#     Xmax = guess_max, 
#     b = guess_b, 
#     c = guess_c)
# 
# 
# n_points <- 15  # Number of experiments
# 
# tol_design <- 1  # We cannot get closer to Xmin or Xmax than this (numerical singularities)
# 
# n_exp <- 500
# sigma_sq_mu <- .2
# 
# do_everything(initial_guss, n_points, tol_design,
#               n_exp, sigma_sq_mu)






