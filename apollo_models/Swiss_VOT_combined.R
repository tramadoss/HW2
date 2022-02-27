# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "Swiss_VOT_combined",
  modelDescr      = "Model from Axhausen",
  nCores          = 4,
  indivID         = "P_NR", 
  outputDirectory = "output"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package

database = readRDS("swiss_vot_rdata.RDS")

database <- database %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(T_DIST_C = if_else(T_DIST_C == 0, 1 , T_DIST_C)) %>%
  inner_join(
    tibble(
      N_O_S = 1:6,
      PT_A = c(0,0,1,0,1,1),
      PT_B = c(1,1,1,0,1,1)
    )
  ) %>%
  mutate(CAR_A = 1-PT_A) %>%
  mutate(CAR_B = 1-PT_B)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(asc_a = 0,
              asc_b = 0,
              
              b_tt_pt_business = -0.1086,
              b_tt_pt_commuters = -0.1353,
              b_tt_pt_leisure = -0.0571,
              b_tt_pt_shopping = -0.1066,
              
              b_tt_car_business = -0.1100,
              b_tt_car_commuters = -0.1491,
              b_tt_car_leisure = -0.0764,
              b_tt_car_shopping = -0.1462,
              
              b_dc = -0.0572,
              
              b_tc_business = -0.1314,
              b_tc_commuters = -0.2920,
              b_tc_leisure = -0.1570,
              b_tc_shopping = -0.3607,
              
              b_ic_business = -1.0318,
              b_ic_commuters = -1.4274,
              b_ic_leisure = -1.1492,
              b_ic_shopping = -1.2692,
              
              b_hw_business = -0.0326,
              b_hw_commuters = -0.0544,
              b_hw_leisure = -0.0350,
              b_hw_shopping = -0.0510,
              
              b_bus_commuters = 2.9428,
              b_bus_leisure = -1.0980,
              b_bus_shopping = 2.4257,
              
              b_rail_business = 1.4948,
              b_rail_commuters = 1.3332,
              b_rail_leisure = 0.2886,
              b_rail_shopping = 1.2117,
              
              b_car_inertia   = 1.9076,
              b_car_available   = 0.5880,
              b_car_male   = -0.3295,
              
              b_bus_disc = 1.1958,
              b_rail_disc = 1.6995,
              b_bus_ga = 3.6511,
              b_rail_ga = 1.7218,
              
              lam_dist_tt_pt_business = -0.3135,
              lam_dist_tt_pt_commuters = -0.2368,
              lam_dist_tt_pt_leisure = -0.2837,
              lam_dist_tt_pt_shopping = -0.2009,
              
              lam_dist_tt_car_business = -0.3573,
              lam_dist_tt_car_commuters = -0.1321,
              lam_dist_tt_car_leisure = -0.3744,
              lam_dist_tt_car_shopping = -0.1905,
              
              lam_dist_tc = -0.5949,
              lam_inc_tc_business = -0.8922,
              lam_inc_tc_commuters = -0.1697,
              
              mu_mc_car_bus = 0.4082,
              mu_mc_car_rail = 0.5051,
              mu_rc_bus = 0.9878,
              mu_rc_car = 1.2254,
              mu_rc_rail_by_car = 0.7688,
              mu_rc_rail = 1
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_a","mu_rc_rail")

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### Variables
  
  dist = T_DIST_C
  tt_pt_elas = # need to multiply by TT_PT
    b_tt_pt_business * PUR_N * (dist/mean(dist))^lam_dist_tt_pt_business +
    b_tt_pt_commuters * PUR_P * (dist/mean(dist))^lam_dist_tt_pt_commuters +
    b_tt_pt_leisure * PUR_T * (dist/mean(dist))^lam_dist_tt_pt_leisure +
    b_tt_pt_shopping * PUR_E * (dist/mean(dist))^lam_dist_tt_pt_shopping
  
  tt_car_elas = # need to multiply by TT_CAR
    b_tt_car_business * PUR_N * (dist/mean(dist))^lam_dist_tt_car_business +
    b_tt_car_commuters * PUR_P * (dist/mean(dist))^lam_dist_tt_car_commuters +
    b_tt_car_leisure * PUR_T * (dist/mean(dist))^lam_dist_tt_car_leisure +
    b_tt_car_shopping * PUR_E * (dist/mean(dist))^lam_dist_tt_car_shopping
  
  inc = HH_INC_A
  tc_elas = 
    b_tc_business * PUR_N * (inc/mean(inc))^lam_inc_tc_business +
    b_tc_commuters * PUR_P * (inc/mean(inc))^lam_inc_tc_commuters +
    b_tc_leisure * PUR_T +
    b_tc_shopping * PUR_E
  
  tc_elas = tc_elas *(dist/mean(dist))^lam_dist_tc 
  
  ic_p =               
    b_ic_business * PUR_N +
    b_ic_commuters * PUR_P +
    b_ic_leisure * PUR_T +
    b_ic_shopping * PUR_E
  
  hw_p = 
    b_hw_business * PUR_N +
    b_hw_commuters * PUR_P +
    b_hw_leisure * PUR_T +
    b_hw_shopping * PUR_E
  
  bus_p = b_bus_commuters * PUR_P + 
    b_bus_leisure * PUR_T + 
    b_bus_shopping * PUR_E
  
  rail_p = 
    b_rail_business * PUR_N + 
    b_rail_commuters * PUR_P + 
    b_rail_leisure * PUR_T + 
    b_rail_shopping * PUR_E
  
  car_mc  = b_car_inertia * CAR   + b_car_available * VEH_AVA1 + b_car_male * MALE
  bus_mc  = b_bus_disc  * DISCOUNT + b_bus_ga  * GA + bus_p
  rail_mc = b_rail_disc * DISCOUNT + b_rail_ga  * GA + rail_p
  
  car_mc = car_mc  * as.integer(N_O_S < 3)
  bus_mc = bus_mc  * as.integer(N_O_S == 1)
  rail_mc= rail_mc * as.integer(N_O_S == 2)
  
  mu = 
    mu_mc_car_bus * as.integer(N_O_S == 1) +
    mu_mc_car_rail * as.integer(N_O_S == 2) +
    mu_rc_bus * as.integer(N_O_S == 3) +
    mu_rc_car * as.integer(N_O_S == 4) +
    mu_rc_rail_by_car * as.integer(N_O_S == 5) +
    mu_rc_rail * as.integer(N_O_S == 6)
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  
  V[["a"]]    = asc_a + b_dc * S_O_C_A * CAR_A  +  tt_pt_elas * TT_A * PT_A  +  
    tt_car_elas * TT_A * CAR_A + tc_elas * TC_A + (ic_p * CHANGE_A + hw_p* HW_B) * PT_A +
    car_mc
  
  V[["b"]]    = asc_b + b_dc * S_O_C_B * CAR_B  +  tt_pt_elas * TT_B * PT_B  +  
    tt_car_elas * TT_B * CAR_B + tc_elas * TC_B + (ic_p * CHANGE_B + hw_p * HW_B) * PT_B +
    bus_mc + rail_mc
  
  V[["a"]]  = mu * V[["a"]]
  V[["b"]]  = mu * V[["b"]]  
  
  mnl_settings = list(
    alternatives  = c(a = 1, b = 2), 
    avail         = list(a = 1, b = 1), 
    choiceVar     = CHOICE,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
  
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model)


