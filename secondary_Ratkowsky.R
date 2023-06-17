
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

## Model definition

model <- "Ratkowsky" 

my_p <- list( 
    Xmin = -1.425, 
    Xmax = 44.36, 
    b = .027, 
    c = .4)

# xopt <- (lambertW0(exp(-my_p$Xmin*my_p$c + my_p$Xmax*my_p$c + 1)) + my_p$c*my_p$Xmin - 1)/my_p$c
# xopt
# 
# my_p$b*(xopt - my_p$Xmin)*(1 - exp(my_p$c*(xopt - my_p$Xmax)))

sigma_mu <- 0.08
n_experiments <- 8000

min_T <- 10
max_T <- 30

maxeval <- 10000

psi <- 50

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

## OEDs fixing Xmax

OED_D_noXmax <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "D",
            include = c("Xmin", "b", "c")
        )
    )

OED_E_noXmax <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "E",
            include = c("Xmin", "b", "c")
        )
    )

OED_penalty_noXmax <- c(min_T:max_T) %>%
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
            include = c("Xmin", "b", "c")
        )
    )

## OEDs fixing Xmax and c

OED_D_noXmax_c <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "D",
            include = c("Xmin", "b")
        )
    )

OED_E_noXmax_c <- c(min_T:max_T) %>%
    set_names(., .) %>%
    map(.,
        ~ calculate_OED(
            model = model, 
            pars = my_p,
            n_points = .,
            tol_design = 1,
            maxeval = maxeval,
            criterion = "E",
            include = c("Xmin", "b")
        )
    )

OED_penalty_noXmax_c <- c(min_T:max_T) %>%
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
            include = c("Xmin", "b")
        )
    )

## Visualize the designs - Figure 1

aa <- unlist(my_p)[c("Xmin", "Xmax")]

p <- list(`D-none` = OED_D,
     `D-Xmax` = OED_D_noXmax,
     `D-Xmax_c` = OED_D_noXmax_c,
     `E-none` = OED_E,
     `E-Xmax` = OED_E_noXmax,
     `E-Xmax_c` = OED_E_noXmax_c,
     `penalty-none` = OED_penalty,
     `penalty-Xmax` = OED_penalty_noXmax,
     `penalty-Xmax_c` = OED_penalty_noXmax_c
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
                          ifelse(fixed == "Xmax", "Xmax fixed", "Xmax & c fixed"))) %>%
    mutate(fixed = factor(fixed, levels = c("None fixed", "Xmax fixed", "Xmax & c fixed"))) %>%
    ggplot() +
    geom_vline(xintercept = aa,
               linetype = 1, colour = "black", linewidth = 1) +
    geom_text(aes(x = x, y = n_points, 
                  label = paste0("(", n, ")"), 
                  colour = n)) +
    facet_grid(design ~ fixed) +
    xlab("Storage temperature (ºC)") + ylab("Number of growth experiments in the design") +
    theme_hc(base_size = 14) +
    scale_colour_gradientn(colours = brewer.pal(8, "YlGnBu")[-c(1,2,3)]) +
    theme(legend.position = "none")

ggsave(p, filename = "Figure1.png",
       width = 9, height = 9)

## Compare with sensitivities  - Sup. figure 2

p1 <- get_sensitivities("Ratkowsky", my_p) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p2 <- get_sensitivities("Ratkowsky", my_p) %>%
    select(-Xmax) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D_noXmax$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E_noXmax$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty_noXmax$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p3 <- get_sensitivities("Ratkowsky", my_p) %>%
    select(-c, -Xmax) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x, value, colour = name)) +
    geom_vline(xintercept = OED_D_noXmax_c$`10`$xbest, linetype = 2, colour = "red") +
    geom_vline(xintercept = OED_E_noXmax_c$`10`$xbest, linetype = 3, colour = "green") +
    geom_vline(xintercept = OED_penalty_noXmax_c$`10`$xbest, linetype = 4, colour = "blue") +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    xlab("Temperature (ºC)") + ylab("Local sensitivity")

p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO")

ggsave(p, filename = "suppFigure2.png",
       width = 9, height = 6,
       bg = "white")

## Aggregated time - Figure 3

all_times <- list(`D-none` = OED_D,
                  `D-Xmax` = OED_D_noXmax,
                  `D-Xmax_c` = OED_D_noXmax_c,
                  `E-none` = OED_E,
                  `E-Xmax` = OED_E_noXmax,
                  `E-Xmax_c` = OED_E_noXmax_c,
                  `penalty-none` = OED_penalty,
                  `penalty-Xmax` = OED_penalty_noXmax,
                  `penalty-Xmax_c` = OED_penalty_noXmax_c
                  ) %>%
    map(.,
        ~ map(., ~ .$xbest)
    ) 

aa <- c(min_T:max_T) %>%
    set_names(., .) %>%
    imap(., ~ seq(my_p$Xmin + 1, my_p$Xmax - 1, length = as.numeric(.y)))

all_times <- c(all_times, list(`uniform-none` = aa,
                               `uniform-Xmax` = aa,
                               `uniform-Xmax_c` = aa))

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
                          ifelse(fixed == "Xmax", "Xmax fixed", "Xmax & c fixed"))) %>%
    ggplot(aes(x = n_points, y = aggreg_time/24, colour = design)) +
    geom_point() +
    geom_line() +
    facet_wrap("fixed") +
    ylab("Aggregated time (days)") +
    scale_y_log10(limits = c(10, 10000)) +
    scale_color_manual(values = my_cols) +
    # theme_bw(base_size = 14) +
    theme_hc(base_size = 14) +
    xlab("Number of growth experiments") +
    theme(legend.title = element_blank(),
          legend.position = "top")


ggsave(p, filename = "Figure3.png",
       width = 12, height = 6,
       bg = "white")

## Simulate experiments 

set.seed(679678)

models_all <- list(D = all_times$`D-none`,
     E = all_times$`E-none`,
     penalty = all_times$`penalty-none`,
     uniform = all_times$`uniform-none`) %>%
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

p0 <- my_p[c("Xmin", "b", "c")]
pfix <- my_p[c("Xmax")]

models_fixXmax <- list(D = all_times$`D-Xmax`,
                       E = all_times$`E-Xmax`,
                       penalty = all_times$`penalty-Xmax`,
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

p0 <- my_p[c("Xmin", "b")]
pfix <- my_p[c("c", "Xmax")]

models_fixXmax_c <- list(D = all_times$`D-Xmax_c`,
                         E = all_times$`E-Xmax_c`,
                         penalty = all_times$`penalty-Xmax_c`,
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

# ## Distribution of the parameters
# 
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
# p2 <- models_fixXmax %>%  # Distribution of parameter estimates
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
# p3 <- models_fixXmax_c %>%  
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
# plot_grid(p1, p2, p3, ncol = 1)
# 
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
# p2 <- models_fixXmax %>% 
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
# p3 <- models_fixXmax_c %>% 
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

## Expected CV - Figure 7

my_cols <- wes_palette("Cavalcanti1", 4)

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
    mutate(cv = abs(`Std. Error`/Estimate)) %>%
    # select(par, cv) %>%
    filter(cv < 1) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9))

p1 <- errors_all %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform"))))%>%
    mutate(fixed = ifelse(par == "none", "None fixed",
                          ifelse(par == "Xmax", "Xmax", "Xmax & c"))) %>%
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

errors_fixXmax <- models_fixXmax %>%  # Distribution of parameter estimates
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
    filter(cv < 1) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9))

p2 <- errors_fixXmax %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform"))))%>%
    mutate(fixed = ifelse(par == "none", "None fixed",
                          ifelse(par == "Xmax", "Xmax", "Xmax & c"))) %>%
    ggplot(aes(x = n_points, y = m_cv, colour = design)) +
    geom_point() +
    geom_line() +
    # geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("par", nrow = 1, scales = "free") +
    xlab("Number of growth experiments") +
    ylab("Expected CV") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = my_cols)

errors_fixXmax_c <- models_fixXmax_c %>%  # Distribution of parameter estimates
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
    # filter(fix == "n_Xmax") %>%
    mutate(cv = abs(`Std. Error`/Estimate)) %>%
    # select(par, cv) %>%
    filter(cv < 1) %>%  # Too large, assumed it did not converge
    group_by(design, n_points, par) %>%
    summarize(m_cv = median(cv, na.rm = TRUE), 
              q10 = quantile(cv, .1),
              q90 = quantile(cv, .9))

p3 <- errors_fixXmax_c %>%
    mutate(design = ifelse(design == "D", "D-optimal",
                           ifelse(design == "E", "E-optimal",
                                  ifelse(design == "penalty", "OED+Penalty", "Uniform"))))%>%
    mutate(fixed = ifelse(par == "none", "None fixed",
                          ifelse(par == "Xmax", "Xmax", "Xmax & c"))) %>%
    ggplot(aes(x = n_points, y = m_cv, colour = design)) +
    geom_point() +
    geom_line() +
    # geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("par", nrow = 1, scales = "free") +
    xlab("Number of growth experiments") +
    ylab("Expected CV") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "right") +
    scale_color_manual(values = my_cols)

p <- plot_grid(p1, p2, p3, ncol = 1, labels = "AUTO")

ggsave(p, filename = "Figure7.png",
       width = 9, height = 12)

## Correlation between parameters - Figure 5

models_all$uniform$`14` %>%
    map_dfr(coef) %>%
    cor(method = "spearman")

p1 <- models_all$uniform$`14` %>%
    map_dfr(coef) %>%
    ggplot(aes(x = Xmin, y = b)) +
    geom_point(shape=1) +
    geom_smooth(color = my_cols[1]) +
    # theme_bw(base_size = 14) +
    xlab("Xmin (ºC)") + ylab("b (log CFU/h/ºC)") +
    theme_hc(base_size = 14)

p2 <- models_all$uniform$`14` %>%
    map_dfr(coef) %>%
    filter(Xmax > 42) %>%
    ggplot(aes(x = Xmax, y = c)) +
    geom_point(shape=1) +
    geom_smooth(color = my_cols[1]) +
    # theme_bw(base_size = 14) +
    xlab("Xmax (ºC)") + ylab("c (1/ºC)") +
    scale_y_log10() +
    theme_hc(base_size = 14)

p <- plot_grid(p1, p2, labels = "AUTO")

ggsave(p, filename = "Figure5.png",
       width = 10, height = 5)


## Distribution of parameter estimates

exp_vals <- models_all %>%
    map(., ~.$`10`) %>%
    map(.,
        ~ map_dfr(., coef)
        ) %>%
    imap_dfr(., ~ mutate(.x, design = .y)) %>%
    filter(c<1.5, b<.5, Xmax<48, Xmin<5) %>%
    pivot_longer(-design) %>%
    group_by(name, design) %>%
    summarize(expected = median(value)) %>%
    left_join(., tibble(name = c("Xmin", "Xmax", "b", "c"),
                        par = c("Xmin (ºC)", "Xmax (ºC)", "b (log CFU/h/ºC)", "c (1/ºC)"))) %>%
    ungroup() %>%
    select(-name)

exp_vals %>%
    pivot_wider(names_from = par, values_from = expected) %>% View()

models_all %>%
    map(., ~.$`20`) %>%
    map(.,
        ~ map_dfr(., coef)
    ) %>%
    imap_dfr(., ~ mutate(.x, design = .y)) %>%
    filter(c<1.5, b<.5, Xmax<48, Xmin<5) %>%
    pivot_longer(-design) %>%
    left_join(., tibble(name = c("Xmin", "Xmax", "b", "c"),
                        par = c("Xmin (ºC)", "Xmax (ºC)", "b (log CFU/h/ºC)", "c (1/ºC)"))) %>%
    ggplot() +
    geom_density(aes(value, colour = design)) +
    facet_wrap("par", scales = "free") +
    geom_vline(aes(xintercept = x),
               linetype = 1,
               size = 1,
               data = tibble(name = names(my_p),
                             x = unlist(my_p)
                             ) %>% left_join(., 
                                             tibble(name = c("Xmin", "Xmax", "b", "c"),
                                                    par = c("Xmin (ºC)", "Xmax (ºC)", "b (log CFU/h/ºC)", "c (1/ºC)"))
                                             )
               ) +
    geom_vline(aes(xintercept = expected, colour = design),
               data = exp_vals,
               linetype = 2,
               size = 1) +
    theme_bw(base_size = 14) +
    xlab("") + ylab("Density") +
    theme(legend.title = element_blank(),
          legend.position = "top")















