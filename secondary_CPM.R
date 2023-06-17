
library(tidyverse)
library(FME)
library(MEIGOR)
library(GGally)
library(lamW)
library(cowplot)
library(RColorBrewer)
library(ggthemes)
library(wesanderson)

source("biodesigner.R")

## 

set.seed(12412)

## Model parameters

model <- "CPM"

my_p <- list(  # Silva et al. (https://doi.org/10.1016/j.foodres.2020.109476)
    Xmin = -1.425, 
    Xopt = 38.17,
    Xmax = 44.36, 
    n = 2,
    mu_opt = 0.976)

sigma_mu <- 0.08
n_experiments <- 8000

min_T <- 10
max_T <- 30

maxeval <- 10000

psi <- 10

## OEDs

OED_D <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "D"
            )
        )

OED_E <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "E"
        )
    )

OED_penalty <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "penalty",
            psi = psi
        )
    )

## OEDs fixing n

OED_D_noN <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "D",
            include = c("Xmin", "Xopt", "Xmax", "mu_opt")
        )
    )

OED_E_noN <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "E",
            include = c("Xmin", "Xopt", "Xmax", "mu_opt")
        )
    )

OED_penalty_noN <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "penalty",
            psi = psi,
            include = c("Xmin", "Xopt", "Xmax", "mu_opt")
        )
    )

## OEDs fixing n and Xmax

OED_D_noN_Xmax <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "D",
            include = c("Xmin", "Xopt", "mu_opt")
        )
    )

OED_E_noN_Xmax <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "E",
            include = c("Xmin", "Xopt", "mu_opt")
        )
    )

OED_penalty_noN_Xmax <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "penalty",
            psi = psi,
            include = c("Xmin", "Xopt", "mu_opt")
        )
    )

## Visualize the designs - Figure 2

aa <- unlist(my_p)[c("Xmin", "Xmax")]
bb <- unlist(my_p)[c("Xopt")]

p <- list(`D-none` = OED_D,
     `D-n` = OED_D_noN,
     `D-n_Xmax` = OED_D_noN_Xmax,
     `E-none` = OED_E,
     `E-n` = OED_E_noN,
     `E-n_Xmax` = OED_E_noN_Xmax,
     `penalty-none` = OED_penalty,
     `penalty-n` = OED_penalty_noN,
     `penalty-n_Xmax` = OED_penalty_noN_Xmax
     ) %>%
    map(.,
        ~ map(., ~ tibble(x = sort(.$xbest)))
        ) %>%
    map(.,
        ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
    ) %>%
    imap_dfr(., ~ mutate(.x, foo = .y)) %>%
    separate(foo, into = c("design", "fixed"), sep = "-") %>%
    mutate(x = round(x, 0)) %>%
    group_by(x, n_points, fixed, design) %>%
    summarize(n = n()) %>% 
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  "OED+Penalty"))) %>%
    mutate(fixed = ifelse(fixed == "none", "None fixed",
                          ifelse(fixed == "n_Xmax", "n & Xmax fixed", "n fixed"))) %>%
    mutate(fixed = factor(fixed, levels = c("None fixed", "n fixed", "n & Xmax fixed"))) %>%
    ggplot() +
    geom_vline(xintercept = aa,
               linetype = 1, colour = "black", linewidth = 1) +
    geom_vline(xintercept = bb,
               linetype = 2, colour = "black", linewidth = 1) +
    geom_text(aes(x = x, y = n_points, 
                  label = paste0("(", n, ")"), 
                  colour = n)) +
    facet_grid(design ~ fixed) +
    # geom_vline(xintercept = aa,
    #            linetype = 1, colour = "red") +
    # geom_vline(xintercept = bb,
    #            linetype = 2, colour = "red") +

    xlab("Storage temperature (ºC)") + ylab("Number of growth experiments in the design") +
    theme_hc(base_size = 14) +
    scale_colour_gradientn(colours = brewer.pal(8, "YlGnBu")[-c(1,2,3)]) +
    theme(legend.position = "none")

ggsave(p, filename = "Figure2.png",
       width = 9, height = 9)

## Compare with sensitivities - Sup. Figure 3

p1 <- get_sensitivities("CPM", my_p) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p2 <- get_sensitivities("CPM", my_p) %>%
    select(-n) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D_noN$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E_noN$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty_noN$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p3 <- get_sensitivities("CPM", my_p) %>%
    select(-n, -Xmax) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D_noN_Xmax$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E_noN_Xmax$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty_noN_Xmax$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO")

ggsave(p, filename = "suppFigure3.png",
       width = 9, height = 6,
       bg = "white")

## Aggregated time - Figure 4

all_times <- list(`D-none` = OED_D,
     `D-n` = OED_D_noN,
     `D-n_Xmax` = OED_D_noN_Xmax,
     `E-none` = OED_E,
     `E-n` = OED_E_noN,
     `E-n_Xmax` = OED_E_noN_Xmax,
     `penalty-none` = OED_penalty,
     `penalty-n` = OED_penalty_noN,
     `penalty-n_Xmax` = OED_penalty_noN_Xmax
) %>%
    map(.,
        ~ map(., ~ .$xbest)
    ) 

aa <- c(min_T:max_T) %>%
    set_names(., .) %>%
    imap(., ~ seq(my_p$Xmin + 1, my_p$Xmax - 1, length = as.numeric(.y)))

all_times <- c(all_times, list(`uniform-none` = aa,
                               `uniform-n` = aa,
                               `uniform-n_Xmax` = aa))

my_cols <- wes_palette("Cavalcanti1", 4)


p <- all_times %>%
    map(.,
        ~ map_dfc(., ~ get_total_time(model, my_p, .))
        ) %>%
    map(.,
        ~ pivot_longer(., everything(), names_to = "n_points", values_to = "aggreg_time")
        ) %>%
    imap_dfr(., ~ mutate(.x, foo = .y, n_points = as.numeric(n_points))) %>%
    separate(foo, into = c("design", "fixed"), sep = "-") %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform"))))%>%
    mutate(fixed = ifelse(fixed == "none", "None fixed",
                          ifelse(fixed == "n_Xmax", "n & Xmax fixed", "n fixed"))) %>%
    mutate(fixed = factor(fixed, levels = c("None fixed", "n fixed", "n & Xmax fixed"))) %>%
    ggplot(aes(x = n_points, y = aggreg_time/24, colour = design)) +
    geom_point() +
    geom_line() +
    facet_wrap("fixed") +
    ylab("Aggregated time (days)") +
    scale_y_log10() +
    scale_color_manual(values = my_cols) +
    # theme_bw(base_size = 14) +
    theme_hc(base_size = 14) +
    xlab("Number of growth experiments") +
    theme(legend.title = element_blank(),
          legend.position = "top")

ggsave(p, filename = "Figure4.png",
       width = 12, height = 6,
       bg = "white")

## Simulate experiments

models_all <- list(D = all_times$`D-none`,
                   E = all_times$`E-none`,
                   penalty = all_times$`penalty-none`,
                   uniform = all_times$`uniform-none`)  %>%
    map(.,
        ~ map(.,
              ~ simulate_experiments(model = model, 
                                     pars_fit = my_p, 
                                     pars_known = list(),  
                                     sampled_x = .,  
                                     n_repetitions = 2,  
                                     sigma_sq_mu = sigma_mu, 
                                     n_exp = n_experiments
              )
        )
    ) 

p0 <- my_p[c("Xmin", "Xopt", "Xmax", "mu_opt")]
pfix <- my_p[c("n")]

models_fixN <- list(D = all_times$`D-n`,
                    E = all_times$`E-n`,
                    penalty = all_times$`penalty-n`,
                    uniform = all_times$`uniform-none`) %>%
    map(.,
        ~ map(.,
              ~ simulate_experiments(model = model, 
                                     pars_fit = p0, 
                                     pars_known = pfix,
                                     sampled_x = .,  
                                     n_repetitions = 2,  
                                     sigma_sq_mu = sigma_mu, 
                                     n_exp = n_experiments
              )
        )
    ) 

p0 <- my_p[c("Xmin", "Xopt", "mu_opt")]
pfix <- my_p[c("n", "Xmax")]

models_fixN_Xmax <- list(D = all_times$`D-n_Xmax`,
                         E = all_times$`E-n_Xmax`,
                         penalty = all_times$`penalty-n_Xmax`,
                         uniform = all_times$`uniform-none`)  %>%
    map(.,
        ~ map(.,
              ~ simulate_experiments(model = model, 
                                     pars_fit = p0, 
                                     pars_known = pfix,
                                     sampled_x = .,  
                                     n_repetitions = 2,  
                                     sigma_sq_mu = sigma_mu, 
                                     n_exp = n_experiments
              )
        )
    ) 

## Distribution of the parameters

# p1 <- models_all %>%  # Distribution of parameter estimates
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ coef(.))
#         )
#     ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#     ) %>%
#     imap_dfr(., ~ mutate(.x, foo = .y)) %>%
#     separate(foo, into = c("design")) %>%
#     pivot_longer(-c(n_points, design)) %>%
#     ggplot() +
#     geom_boxplot(aes(y = value, x = factor(n_points), colour = design)) +
#     facet_wrap("name", scales = "free_y", nrow = 1) +
#     geom_hline(aes(yintercept = x), linetype = 2,
#                data = tibble(x = unlist(my_p), name = names(my_p))
#     )
# 
# p2 <- models_fixN %>%  # Distribution of parameter estimates
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ coef(.))
#         )
#     ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#     ) %>%
#     imap_dfr(., ~ mutate(.x, foo = .y)) %>%
#     separate(foo, into = c("design")) %>%
#     pivot_longer(-c(n_points, design)) %>%
#     ggplot() +
#     geom_boxplot(aes(y = value, x = factor(n_points), colour = design)) +
#     facet_wrap("name", scales = "free_y", nrow = 1) +
#     geom_hline(aes(yintercept = x), linetype = 2,
#                data = tibble(x = unlist(my_p), name = names(my_p))
#     )
# 
# 
# p3 <- models_fixN_Xmax %>%
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ coef(.))
#               )
#         ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#         ) %>%
#     imap_dfr(., ~ mutate(.x, foo = .y)) %>%
#     separate(foo, into = c("design")) %>%
#     pivot_longer(-c(n_points, design)) %>%
#     ggplot() +
#     geom_boxplot(aes(y = value, x = factor(n_points), colour = design)) +
#     facet_wrap("name", scales = "free_y", nrow = 1) +
#     geom_hline(aes(yintercept = x), linetype = 2,
#                data = tibble(x = unlist(my_p), name = names(my_p))
#                )
# 
# plot_grid(p1, p2, p3, ncol = 1)

# ## Distribution of the CV
# 
# p1 <- models_all %>%
#     map(.,
#         ~ map(.,
#               ~ map(., ~ summary(.)$par)
#         )
#     ) %>%
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ as_tibble(., rownames = "par"))
#         )
#     ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#     )  %>%
#     imap_dfr(.,
#              ~ mutate(.x, foo = .y)
#     ) %>%
#     separate(foo, into = c("design")) %>%
#     mutate(cv = abs(`Std. Error`/Estimate)) %>%
#     # select(par, cv) %>%
#     filter(cv < 1) %>%  # Too large, assumed it did not converge
#     ggplot() +
#     geom_boxplot(aes(y = cv, x = factor(n_points), colour = design)) +
#     facet_wrap("par")
# 
# p2 <- models_fixN %>%
#     map(.,
#         ~ map(.,
#               ~ map(., ~ summary(.)$par)
#         )
#     ) %>%
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ as_tibble(., rownames = "par"))
#         )
#     ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#     )  %>%
#     imap_dfr(.,
#              ~ mutate(.x, foo = .y)
#     ) %>%
#     separate(foo, into = c("design")) %>%
#     mutate(cv = abs(`Std. Error`/Estimate)) %>%
#     # select(par, cv) %>%
#     filter(cv < 1) %>%  # Too large, assumed it did not converge
#     ggplot() +
#     geom_boxplot(aes(y = cv, x = factor(n_points), colour = design)) +
#     facet_wrap("par")
# 
# p3 <- models_fixN_Xmax %>%
#     map(.,
#         ~ map(.,
#               ~ map(., ~ summary(.)$par)
#         )
#     ) %>%
#     map(.,
#         ~ map(.,
#               ~ map_dfr(., ~ as_tibble(., rownames = "par"))
#         )
#     ) %>%
#     map(.,
#         ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
#     )  %>%
#     imap_dfr(.,
#              ~ mutate(.x, foo = .y)
#     ) %>%
#     separate(foo, into = c("design")) %>%
#     mutate(cv = abs(`Std. Error`/Estimate)) %>%
#     # select(par, cv) %>%
#     filter(cv < 1) %>%  # Too large, assumed it did not converge
#     ggplot() +
#     geom_boxplot(aes(y = cv, x = factor(n_points), colour = design)) +
#     facet_wrap("par")
# 
# plot_grid(p1, p2, p3, ncol = 1)

## Expected CV - Figure 8
    
errors_all <- models_all %>%  # Distribution of parameter estimates
    map(.,
        ~ map(., 
              ~ map(., ~ summary(.)$par)
        )
    ) %>%
    map(.,
        ~ map(.,
              ~ map_dfr(., ~ as_tibble(., rownames = "par"))
        )
    ) %>%
    map(.,
        ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
    )  %>%
    imap_dfr(.,
             ~ mutate(.x, foo = .y)
    ) %>%
    separate(foo, into = c("design", "fix"), sep = "-") %>%
    # filter(fix == "none") %>%
    mutate(cv = abs(`Std. Error`/Estimate)) %>%
    # select(par, cv) %>%
    filter(cv < 3) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9),
              n = n())

my_cols <- wes_palette("Cavalcanti1", 4)

p1 <- errors_all %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform")))) %>%
    ggplot(aes(x = n_points, y = m_cv, colour = design)) +
    geom_point() +
    geom_line() +
    # geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("par", nrow = 1, scales = "free") +
    # theme_bw(base_size = 14) +
    xlab("Number of growth experiments") +
    ylab("Expected CV") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = my_cols)

errors_fixN <- models_fixN %>%  # Distribution of parameter estimates
    map(.,
        ~ map(., 
              ~ map(., ~ summary(.)$par)
        )
    ) %>%
    map(.,
        ~ map(.,
              ~ map_dfr(., ~ as_tibble(., rownames = "par"))
        )
    ) %>%
    map(.,
        ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
    )  %>%
    imap_dfr(.,
             ~ mutate(.x, foo = .y)
    ) %>%
    separate(foo, into = c("design", "fix"), sep = "-") %>%
    # filter(fix == "n") %>%
    mutate(cv = abs(`Std. Error`/Estimate)) %>%
    # select(par, cv) %>%
    filter(cv < 3) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9))

p2 <- errors_fixN %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform")))) %>%
    ggplot(aes(x = n_points, y = m_cv, colour = design)) +
    geom_point() +
    geom_line() +
    # geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("par", nrow = 1, scales = "free") +
    # theme_bw(base_size = 14) +
    xlab("Number of growth experiments") +
    ylab("Expected CV") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = my_cols)

errors_fixN_Xmax <- models_fixN_Xmax %>%  # Distribution of parameter estimates
    map(.,
        ~ map(., 
              ~ map(., ~ summary(.)$par)
        )
    ) %>%
    map(.,
        ~ map(.,
              ~ map_dfr(., ~ as_tibble(., rownames = "par"))
        )
    ) %>%
    map(.,
        ~ imap_dfr(., ~ mutate(.x, n_points = as.numeric(.y)))
    )  %>%
    imap_dfr(.,
             ~ mutate(.x, foo = .y)
    ) %>%
    separate(foo, into = c("design", "fix"), sep = "-") %>%
    mutate(cv = abs(`Std. Error`/Estimate)) %>%
    # select(par, cv) %>%
    filter(cv < 3) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9))
    
p3 <- errors_fixN_Xmax %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform")))) %>%
    ggplot(aes(x = n_points, y = m_cv, colour = design)) +
    geom_point() +
    geom_line() +
    # geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("par", nrow = 1, scales = "free") +
    # theme_bw(base_size = 14) +
    xlab("Number of growth experiments") +
    ylab("Expected CV") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "right") +
    scale_color_manual(values = my_cols)

p <- plot_grid(p1, p2, p3, ncol = 1, labels = "AUTO")


ggsave(p, filename = "Figure8.png",
       width = 9, height = 12)

## Correlation between parameters - Figure 6

models_fixN$uniform$`20` %>%
    map_dfr(coef) %>%
    # filter(Xopt < 42) %>%
    cor(method = "spearman")

p1 <- models_fixN$uniform$`20` %>%
    map_dfr(coef) %>%
    filter(Xmax < 54) %>%
    # filter(Xopt < 42) %>%
    ggplot(aes(x = Xopt, y = Xmax)) +
    geom_point(shape = 1) +
    geom_smooth(color = my_cols[1]) +
    xlab("Topt (ºC)") + ylab("Xmax (ºC") +
    theme_hc(base_size = 14)

p2 <- models_fixN$uniform$`20` %>%
    map_dfr(coef) %>%
    filter(Xmax < 54) %>%
    # filter(Xopt < 42) %>%
    ggplot(aes(x = mu_opt, y = Xmax)) +
    geom_point(shape = 1) +
    geom_smooth(color = my_cols[1]) +
    xlab("mu_opt (log CFU/h)") + ylab("Xmax (ºC)") + 
    theme_hc(base_size = 14)

p <- cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(p, filename = "Figure6.png",
       width = 10, height = 5)

## Correlations when fitting n - Supp Figure 4

models_all$uniform$`20` %>%
    map_dfr(coef) %>%
    # filter(n<5) %>%
    cor(method = "spearman")


p1 <- models_all$uniform$`20` %>%
    map_dfr(coef) %>%
    # filter(n<5) %>%
    ggplot(aes(x = n, y = Xmin)) +
    geom_point(shape=1) +
    geom_smooth(color = my_cols[1]) +
    theme_hc(base_size = 14) +
    xlab("n") + ylab("Xmin (ºC)") +
    coord_cartesian(ylim = c(-10, 5), xlim = c(1, 5))

p3 <- models_all$uniform$`20` %>%
    map_dfr(coef) %>%
    filter(Xmax < 51) %>%
    # filter(Xmin > -50, n<5) %>%
    ggplot(aes(x = Xopt, y = Xmax)) +
    geom_point(shape=1) +
    geom_smooth(color = my_cols[1]) +
    theme_hc(base_size = 14) + 
    xlab("Xopt (ºC)") + ylab("Xmax (ºC)")

p2 <- models_all$uniform$`20` %>%
    map_dfr(coef) %>%
    filter(n < 5) %>%
    # filter(Xmin > -50, n<5) %>%
    # filter(Xmin > -50) %>%
    ggplot(aes(x = n, y = Xopt)) +
    geom_point(shape=1) +
    geom_smooth(color = my_cols[1]) +
    theme_hc(base_size = 14) +
    xlab("n") + ylab("Xopt (ºC)")

p <- plot_grid(p1, p2, p3,
          nrow = 1, labels = "AUTO")

ggsave(p, filename = "suppFigure4.png",
       width = 10, height = 5)

# ## Distribution of parameter estimates
# 
# exp_vals <- models_fixN %>%
#     map(., ~.$`10`) %>%
#     map(.,
#         ~ map_dfr(., coef)
#     ) %>%
#     imap_dfr(., ~ mutate(.x, design = .y)) %>%
#     # filter(c<1.5, b<.5, Xmax<48, Xmin<5) %>%
#     pivot_longer(-design) %>%
#     group_by(name, design) %>%
#     summarize(expected = median(value)) %>%
#     left_join(., tibble(name = c("Xmin", "Xmax", "Xopt", "mu_opt"),
#                         par = c("Tmin (ºC)", "Xmax (ºC)", "Topt (ºC)", "mu_opt (log CFU/h)"))) %>%
#     ungroup() %>%
#     select(-name)
# 
# exp_vals %>%
#     pivot_wider(names_from = par, values_from = expected) %>% View()
# 
# models_fixN %>%
#     map(., ~.$`20`) %>%
#     map(.,
#         ~ map_dfr(., coef)
#     ) %>%
#     imap_dfr(., ~ mutate(.x, design = .y)) %>%
#     filter(mu_opt<10, Xmax < 48, Xopt < 50, Xmin<5) %>%
#     pivot_longer(-design) %>%
#     left_join(., tibble(name = c("Xmin", "Xmax", "Xopt", "mu_opt"),
#                         par = c("Tmin (ºC)", "Xmax (ºC)", "Topt (ºC)", "mu_opt (log CFU/h)"))) %>%
#     ggplot() +
#     geom_density(aes(value, colour = design)) +
#     facet_wrap("par", scales = "free") +
#     geom_vline(aes(xintercept = x),
#                linetype = 1,
#                size = 1,
#                data = tibble(par = c("Tmin (ºC)", "Xmax (ºC)", "Topt (ºC)", "mu_opt (log CFU/h)"),
#                              x = c(my_p$Xmin, my_p$Xmax, my_p$Xopt, my_p$mu_opt)
#                )
#     ) +
#     geom_vline(aes(xintercept = expected, colour = design),
#                data = exp_vals,
#                linetype = 2,
#                size = 1) +
#     theme_bw(base_size = 14) +
#     xlab("") + ylab("Density") +
#     theme(legend.title = element_blank(),
#           legend.position = "top")



