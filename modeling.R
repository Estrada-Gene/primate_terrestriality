####### statistical modeling for Estrada and Marshall (2024) #######

library(tidyverse)
library(phytools)
library(brms)
library(brmstools)
library(performance)
library(cowplot)
library(modelsummary)
library(kableExtra)

#read in primate data 
primates <- read.csv(file = "data/primates_final.csv", header = TRUE)

#saving ordinal terrestriality variable as ordered factor
primates$Terrestrial.ordinal <- factor(primates$Terrestrial.ordinal, 
                                       levels = c("arboreal", "semiterrestrial", "terrestrial"),
                                       ordered = TRUE)

#read in '10k Trees' primate phylogenetic tree
full_tree <- read.nexus(file = "data/10kTrees_Primates_full.nex")

#remove underscores from Latin binomials
full_tree$tip.label <- sub("_"," ", full_tree$tip.label)
full_tree$tip.label <- sub("_"," ", full_tree$tip.label)

#create variance-covariance matrix
A <- vcv.phylo(full_tree)

#two-part function to create multi-model coefficient plots for brms objects
compare_posteriors_data <- function(..., dodge_width = 0.5) {
  dots <- rlang::dots_list(..., .named = TRUE)
  draws <- lapply(dots, function(x) {
    if (class(x)[1] == "stanreg") {
      posterior::subset_draws(posterior::as_draws(x$stanfit),
                              variable = names(fixef(x))
      )
    } else if (class(x)[1] == "brmsfit") {
      brm_draws <- posterior::subset_draws(posterior::as_draws(x$fit),
                                           variable = paste0("b_", rownames(fixef(x)))
      )
      posterior::variables(brm_draws) <- stringr::str_split(posterior::variables(brm_draws), "_", simplify = T)[, 2]
      posterior::rename_variables(brm_draws, `(Intercept)` = Intercept)
    } else {
      stop(paste0(class(x)[1], " objects not supported."))
    }
  })
  intervals <- lapply(draws, bayesplot::mcmc_intervals_data)
  combined <- dplyr::bind_rows(intervals, .id = "model")
  return(combined) 
}

compare_posteriors_plots <- function(posterior_data, parameters_keep, dodge_width = -0.5){
  new_posterior_data <- posterior_data %>%
    filter(parameter %in% parameters_keep)
  ggplot(new_posterior_data, aes(x = m, y = parameter, color = model, group = model)) +
    geom_linerange(aes(xmin = l, xmax = h), size = 2, position = position_dodge(dodge_width)) +
    geom_linerange(aes(xmin = ll, xmax = hh), position = position_dodge(dodge_width)) +
    geom_point(color = "black", position = position_dodge(dodge_width)) +
    geom_vline(xintercept = 0, linetype = "dashed")
}

#set seed
set.seed(1120)

#########################################################################################################################
####### MODELS 

##### binary models

#sex-averaged model
brms1 <- brm(Terrestrial.binary ~
               Estrada.Continent +                           #continent primate located in
               AvBody.mass.comb.sc +                         #sex-averaged body mass
               n.carn.sc  +                                  #number of carnivores in range & >25% primate mass
               Perc.fruit.comb.sc +                          #percent fruit in diet
               WC1_area_mean.sc +                            #max temperature
               WC12_area_mean.sc +                           #annual rainfall   
               FH.Median.sc +                                #forest height
               NDVI.MED.sc +                                 #median NDVI
               TC.median.sc +                                #canopy cover
               (1|gr(TenK.tax, cov = A)),                    #control for phylogeny
             family = bernoulli(link = "logit"),
             data = primates,
             data2 = list(A = A),
             iter = 6000, warmup = 2000, thin = 1,
             prior = c(prior(normal(0,10), class = "Intercept"),
                       prior(normal(0,1), class = "b")),
             control = list(adapt_delta = 0.9999),
             cores = 12)

#model summary
summary(brms1)
#checking for multicollinearity 
check_collinearity(brms1)
#posterior predictive check 
pp_check(brms1, ndraws = 1000)


#female-only model
brms1.1 <- brm(Terrestrial.binary ~
                 Estrada.Continent +                           
                 FBody.mass.comb.sc +     #female body mass only
                 n.carn.sc.f  +           #number of carnivores in range & >25% female mass
                 Perc.fruit.comb.sc +                          
                 WC1_area_mean.sc +                            
                 WC12_area_mean.sc +                             
                 FH.Median.sc +                                
                 NDVI.MED.sc +
                 TC.median.sc +
                 (1|gr(TenK.tax, cov = A)),
               family = bernoulli(link = "logit"),
               data = primates,
               data2 = list(A = A),
               iter = 6000, warmup = 2000, thin = 1,
               prior = c(prior(normal(0,10), class = "Intercept"),
                         prior(normal(0,1), class = "b")),
               control = list(adapt_delta = 0.9999),
               cores = 12)

summary(brms1.1)
check_collinearity(brms1.1)
pp_check(brms1.1, ndraws = 1000)

#compiling brms model output parameters into one object 
cpd1 <- compare_posteriors_data(brms1, brms1.1)

#changing parameter names for two female-only pars so they line up with the associated pars in sex-avg model 
cpd1 <- cpd1 %>%
  mutate(parameter = gsub("n.carn.sc.f", "n.carn.sc", parameter)) %>%
  mutate(parameter = gsub("FBody.mass.comb.sc", "AvBody.mass.comb.sc", parameter)) %>%
  mutate(parameter = fct_relevel(parameter, 
                                 "AvBody.mass.comb.sc", "Perc.fruit.comb.sc", "FH.Median.sc", 
                                 "TC.median.sc", "NDVI.MED.sc", "WC1", "WC12", "n.carn.sc"))

#plotting both models together
p1 <- compare_posteriors_plots(cpd1,
                               c("AvBody.mass.comb.sc", "n.carn.sc", "Perc.fruit.comb.sc",
                                 "WC1", "WC12", "NDVI.MED.sc", "FH.Median.sc", "TC.median.sc"))

p1 +
  labs(title = "Binary Models", x = "coefficient estimate", y = "", fill = "") + 
  scale_y_discrete(labels = c("AvBody.mass.comb.sc" = "body mass",
                              "Perc.fruit.comb.sc" = "% fruit in diet",
                              "FH.Median.sc" = "forest height",
                              "TC.median.sc" = "canopy cover",
                              "NDVI.MED.sc" = "NDVI",
                              "WC1" = "max. temp.",
                              "WC12" = "annual rainfall",
                              "n.carn.sc" = "# carnivores in range"), 
                   limits = rev) + 
  scale_colour_discrete(name = "", labels = c("m/f avg.","female only")) +
  theme_classic() + 
  theme(legend.position = c(.85, 0.5), 
        legend.background = element_blank())


##### ordinal models

#sex-averaged model
brms2 <- brm(Terrestrial.ordinal ~
               Estrada.Continent +                           
               AvBody.mass.comb.sc +
               n.carn.sc  +                         
               Perc.fruit.comb.sc +                          
               WC1_area_mean.sc +                            
               WC12_area_mean.sc +                              
               FH.Median.sc +                                
               NDVI.MED.sc +                                 
               TC.median.sc +
               (1|gr(TenK.tax, cov = A)),
             family = cumulative(link = "probit"),
             data = primates,
             data2 = list(A = A),
             iter = 6000, warmup = 2000, thin = 1,
             prior = c(prior(normal(0,10), class = "Intercept"),
                       prior(normal(0,1), class = "b")),
             control = list(adapt_delta = 0.9999),
             cores = 12)

summary(brms2)
check_collinearity(brms2)
pp_check(brms2, ndraws = 1000)

#female-only model
brms2.1 <- brm(Terrestrial.ordinal ~
                 Estrada.Continent +                           
                 FBody.mass.comb.sc +
                 n.carn.sc.f  +                         
                 Perc.fruit.comb.sc +                          
                 WC1_area_mean.sc +                            
                 WC12_area_mean.sc +                             
                 FH.Median.sc +                                
                 NDVI.MED.sc +                                 
                 TC.median.sc +
                 (1|gr(TenK.tax, cov = A)),
               family = cumulative(link = "probit"),
               data = primates,
               data2 = list(A = A),
               iter = 6000, warmup = 2000, thin = 1,
               prior = c(prior(normal(0,10), class = "Intercept"),
                         prior(normal(0,1), class = "b")),
               control = list(adapt_delta = 0.9999),
               cores = 12)

summary(brms2.1)
check_collinearity(brms2.1)
pp_check(brms2.1, ndraws = 1000)

#compiling brms model output parameters into one object 
cpd2 <- compare_posteriors_data(brms2, brms2.1)

#changing parameter names
cpd2 <- cpd2 %>%
  mutate(parameter = gsub("n.carn.sc.f", "n.carn.sc", parameter)) %>%
  mutate(parameter = gsub("FBody.mass.comb.sc", "AvBody.mass.comb.sc", parameter)) %>%
  mutate(parameter = fct_relevel(parameter, 
                                 "AvBody.mass.comb.sc", "Perc.fruit.comb.sc", "FH.Median.sc", "TC.median.sc",
                                 "NDVI.MED.sc", "WC1", "WC12", "n.carn.sc"))

#plotting both models
p2 <- compare_posteriors_plots(cpd2,
                                  c("AvBody.mass.comb.sc", "n.carn.sc", "Perc.fruit.comb.sc",
                                    "WC1", "WC12", "NDVI.MED.sc", "FH.Median.sc", "TC.median.sc"))

p2 +
  labs(title = "Ordinal Models", x = "coefficient estimate", y = "", fill = "") + 
  scale_y_discrete(labels = c("AvBody.mass.comb.sc" = "body mass",
                              "Perc.fruit.comb.sc" = "% fruit in diet",
                              "FH.Median.sc" = "forest height",
                              "TC.median.sc" = "canopy cover",
                              "NDVI.MED.sc" = "NDVI",
                              "WC1" = "max. temp.",
                              "WC12" = "annual rainfall",
                              "n.carn.sc" = "# carnivores in range"), 
                   limits = rev) + 
  scale_colour_discrete(name = "", labels = c("m/f avg.","female only")) +
  theme_classic() + 
  theme(legend.position = c(.85, 0.5), 
        legend.background = element_blank())


##### continuous models

#sex-averaged model
brms3 <- brm(ground.use.prop ~ 
               Estrada.Continent +                           
               AvBody.mass.comb.sc +
               n.carn.sc  +                         
               Perc.fruit.comb.sc +                          
               WC1_area_mean.sc +                            
               WC12_area_mean.sc +                              
               FH.Median.sc +                                
               NDVI.MED.sc +                                 
               TC.median.sc +
               (1|gr(TenK.tax, cov = A)),
             family = zero_inflated_beta(link = "logit"),
             data = primates,
             data2 = list(A = A),
             iter = 8000, warmup = 3000, thin = 1,
             prior = c(prior(normal(0,10), class = "Intercept"),
                       prior(normal(0,1), class = "b")),
             control = list(adapt_delta = 0.9999),
             cores = 12)

summary(brms3)
check_collinearity(brms3)
pp_check(brms3, ndraws = 1000)

#female-only model
brms3.1 <- brm(ground.use.prop ~ 
                 Estrada.Continent +                           
                 FBody.mass.comb.sc +                          
                 Perc.fruit.comb.sc +                          
                 WC1_area_mean.sc +                            
                 WC12_area_mean.sc +                             
                 FH.Median.sc +                                
                 NDVI.MED.sc +                                 
                 n.carn.sc.f +                                 
                 TC.median.sc +
                 (1|gr(TenK.tax, cov = A)),
               family = zero_inflated_beta(link = "logit"),
               data = primates,
               data2 = list(A = A),
               iter = 8000, warmup = 3000, thin = 1,
               prior = c(prior(normal(0,10), class = "Intercept"),
                         prior(normal(0,1), class = "b")),
               control = list(adapt_delta = 0.9999),
               cores = 12)

summary(brms3.1)
check_collinearity(brms3.1)
pp_check(brms3.1, ndraws = 1000)


#compiling brms model output parameters
cpd3 <- compare_posteriors_data(brms3, brms3.1)

#changing parameter names
cpd3 <- cpd3 %>%
  mutate(parameter = gsub("n.carn.sc.f", "n.carn.sc", parameter)) %>%
  mutate(parameter = gsub("FBody.mass.comb.sc", "AvBody.mass.comb.sc", parameter)) %>%
  mutate(parameter = fct_relevel(parameter, 
                                 "AvBody.mass.comb.sc", "Perc.fruit.comb.sc", "FH.Median.sc", "TC.median.sc",
                                 "NDVI.MED.sc", "WC1", "WC12", "n.carn.sc"))

#plotting both models
p3 <- compare_posteriors_plots(cpd3, 
                                  c("AvBody.mass.comb.sc", "n.carn.sc", "Perc.fruit.comb.sc",
                                    "WC1", "WC12", "NDVI.MED.sc", "FH.Median.sc", "TC.median.sc"))

p3 +
  labs(title = "Continuous Models", x = "coefficient estimate", y = "", fill = "") + 
  scale_y_discrete(labels = c("AvBody.mass.comb.sc" = "body mass",
                              "Perc.fruit.comb.sc" = "% fruit in diet",
                              "FH.Median.sc" = "forest height",
                              "TC.median.sc" = "canopy cover",
                              "NDVI.MED.sc" = "NDVI",
                              "WC1" = "max. temp.",
                              "WC12" = "annual rainfall",
                              "n.carn.sc" = "# carnivores in range"), 
                   limits = rev) + 
  scale_colour_discrete(name = "", labels = c("m/f avg.","female only")) +
  theme_classic() + 
  theme(legend.position = c(.85, 0.5), 
        legend.background = element_blank())


##### sub-model one

#sex-averaged model
brms_s1 <- brm(ground.use.prop ~ 
                 AvBody.mass.comb.sc +                         
                 Perc.fruit.comb.sc +                          
                 WC1_area_mean.sc +                            
                 WC12_area_mean.sc +                              
                 FH.Median.sc +                                
                 NDVI.MED.sc +                                 
                 n.carn.sc +                                   
                 TC.median.sc +
                 (1|gr(TenK.tax, cov = A)),
               family = zero_inflated_beta(link = "logit"),
               data = primates[primates$radiation == "Cercopithecoidea" | primates$radiation == "Hominoidea",], 
               data2 = list(A = A),
               iter = 8000, warmup = 3000, thin = 1,
               prior = c(prior(normal(0,10), class = "Intercept"),
                         prior(normal(0,1), class = "b")),
               control = list(adapt_delta = 0.9999),
               cores = 12)

summary(brms_s1)
check_collinearity(brms_s1)
pp_check(brms_s1, ndraws = 1000)

#female
brms_s1f <- brm(ground.use.prop ~ 
                  FBody.mass.comb.sc +                          
                  Perc.fruit.comb.sc +                          
                  WC1_area_mean.sc +                            
                  WC12_area_mean.sc +                             
                  FH.Median.sc +                                
                  NDVI.MED.sc +                                 
                  n.carn.sc.f +                                 
                  TC.median.sc +
                  (1|gr(TenK.tax, cov = A)),
                family = zero_inflated_beta(link = "logit"),
                data = primates[primates$radiation == "Cercopithecoidea" | primates$radiation == "Hominoidea",], 
                data2 = list(A = A),
                iter = 8000, warmup = 3000, thin = 1,
                prior = c(prior(normal(0,10), class = "Intercept"),
                          prior(normal(0,1), class = "b")),
                control = list(adapt_delta = 0.9999),
                cores = 12)

summary(brms_s1f)
check_collinearity(brms_s1f)
pp_check(brms_s1f, ndraws = 1000)


#compiling brms model output parameters
cpd_s1 <- compare_posteriors_data(brms_s1, brms_s1f)

#changing parameter names for two female-only models pars so that they line up with the associated pars in full model 
cpd_s1 <- cpd_s1 %>%
  mutate(parameter = gsub("n.carn.sc.f", "n.carn.sc", parameter)) %>%
  mutate(parameter = gsub("FBody.mass.comb.sc", "AvBody.mass.comb.sc", parameter)) %>%
  mutate(parameter = fct_relevel(parameter, 
                                 "AvBody.mass.comb.sc", "Perc.fruit.comb.sc", "FH.Median.sc", "TC.median.sc",
                                 "NDVI.MED.sc", "WC1", "WC12", "n.carn.sc"))

#plotting 
p4 <- compare_posteriors_plots(cpd_s1,
                                  c("AvBody.mass.comb.sc", "n.carn.sc", "Perc.fruit.comb.sc",
                                    "WC1", "WC12", "NDVI.MED.sc", "FH.Median.sc", "TC.median.sc"))

p4 +
  labs(title = "Sub-model 1", subtitle = "Cercopithecoidea & Hominoidea",
       x = "coefficient estimate", y = "") + 
  scale_y_discrete(labels = c("AvBody.mass.comb.sc" = "body mass",
                              "Perc.fruit.comb.sc" = "% fruit in diet",
                              "FH.Median.sc" = "forest height",
                              "TC.median.sc" = "canopy cover",
                              "NDVI.MED.sc" = "NDVI",
                              "WC1" = "max. temp.",
                              "WC12" = "annual rainfall",
                              "n.carn.sc" = "# carnivores in range"), 
                   limits = rev) + 
  xlim(-2,2.3) +
  scale_colour_discrete(name = "", labels = c("m/f avg.","female only")) +
  theme_classic() + 
  theme(legend.position = "none")


##### sub-model two

#sex-averaged model
brms_s2 <- brm(ground.use.prop ~ 
                 AvBody.mass.comb.sc +                         
                 Perc.fruit.comb.sc +                          
                 WC1_area_mean.sc +                            
                 WC12_area_mean.sc +                             
                 FH.Median.sc +                                
                 NDVI.MED.sc +                                 
                 n.carn.sc +                                   
                 TC.median.sc +
                 (1|gr(TenK.tax, cov = A)),
               family = zero_inflated_beta(link = "logit"),
               data = primates[primates$radiation == "Platyrrhines" | primates$radiation == "Lemuriformes",], 
               data2 = list(A = A),
               iter = 8000, warmup = 3000, thin = 1,
               prior = c(prior(normal(0,10), class = "Intercept"),
                         prior(normal(0,1), class = "b")),
               control = list(adapt_delta = 0.9999),
               cores = 12)

summary(brms_s2)
check_collinearity(brms_s2)
pp_check(brms_s2, ndraws = 1000)

#female
brms_s2f <- brm(ground.use.prop ~ 
                  FBody.mass.comb.sc +                          
                  Perc.fruit.comb.sc +                          
                  WC1_area_mean.sc +                            
                  WC12_area_mean.sc +                             
                  FH.Median.sc +                                
                  NDVI.MED.sc +                                 
                  n.carn.sc.f +                                
                  TC.median.sc +
                  (1|gr(TenK.tax, cov = A)),
                family = zero_inflated_beta(link = "logit"),
                data = primates[primates$radiation == "Platyrrhines" | primates$radiation == "Lemuriformes",], 
                data2 = list(A = A),
                iter = 8000, warmup = 3000, thin = 1,
                prior = c(prior(normal(0,10), class = "Intercept"),
                          prior(normal(0,1), class = "b")),
                control = list(adapt_delta = 0.9999),
                cores = 12)

summary(brms_s2f)
check_collinearity(brms_s2f)
pp_check(brms_s2f, ndraws = 1000)

#compiling brms model output parameters
cpd_s2 <- compare_posteriors_data(brms_s2, brms_s2f)

#changing parameter names
cpd_s2 <- cpd_s2 %>%
  mutate(parameter = gsub("n.carn.sc.f", "n.carn.sc", parameter)) %>%
  mutate(parameter = gsub("FBody.mass.comb.sc", "AvBody.mass.comb.sc", parameter)) %>%
  mutate(parameter = fct_relevel(parameter, 
                                 "AvBody.mass.comb.sc", "Perc.fruit.comb.sc", "FH.Median.sc", "TC.median.sc",
                                 "NDVI.MED.sc", "WC1", "WC12", "n.carn.sc"))

#plotting 
p5 <- compare_posteriors_plots(cpd_s2,
                                  c("AvBody.mass.comb.sc", "n.carn.sc", "Perc.fruit.comb.sc",
                                    "WC1", "WC12", "NDVI.MED.sc", "FH.Median.sc", "TC.median.sc"))
p5 +
  labs(title = "Sub-model 2", subtitle = "Platyrrhines & Lemuriformes",
       x = "coefficient estimate", y = "") + 
  scale_y_discrete(limits = rev) + 
  scale_colour_discrete(name = "", labels = c("m/f avg.","female only")) +
  #xlim(-2, 2) +
  theme_classic() + 
  theme(legend.position = c(.85, 0.5), 
        legend.background = element_blank(),
        axis.text.y=element_blank())


#comparing sub-models
plot_grid(p4, p5, 
          align = "h", 
          nrow = 1, 
          rel_widths = c(3/5, 2/5))


##### model comparison table

table1 <- modelsummary(list("(m/f)" =  brms1,     #binary m/f
                            "(fem)" =  brms1.1,   #binary fem
                            "(m/f)" =  brms2,     #ordinal m/f
                            "(fem)" =  brms2.1,   #binary fem
                            "(m/f)" =  brms3,     #continuous m/f
                            "(fem)" =  brms3.1,   #continuous fem
                            "(m/f)" =  brms_s1,   #sub1 m/f
                            "(fem)" =  brms_s1f,  #sub1 fem
                            "(m/f)" =  brms_s2,   #sub2 m/f
                            "(fem)" =  brms_s2f), #sub2 fem
                       fmt = 3,
                       centrality = "mean",
                       estimate = "{estimate} ({std.dev})",
                       gof_map = c("nobs", "r.squared"),
                       coef_map = c("b_AvBody.mass.comb.sc" = "Body Mass (sex-averaged)",
                                    "b_FBody.mass.comb.sc" = "Body Mass (female only)",
                                    "b_Perc.fruit.comb.sc" = "Percent Fruit in Diet",
                                    "b_FH.Median.sc" = "Forest Height",
                                    "b_TC.median.sc" = "Canopy Cover",
                                    "b_NDVI.MED.sc" = "NDVI",
                                    "b_WC1_area_mean.sc" = "Maximum Temperature",
                                    "b_WC12_area_mean.sc" = "Annual Rainfall",
                                    "b_n.carn.sc" = "Carnivores in Range (sex-averaged)",
                                    "b_n.carn.sc.f" = "Carnivores in Range (female only)"))

#add headers
kable1 <- add_header_above(table1, c(" ", "Binary" = 2, "Ordinal" = 2, "Continuous" = 2, 
                                     "Sub-model 1" = 2, "Sub-model 2" = 2))

#save
#save_kable(kable1, file = "kable1.png")
