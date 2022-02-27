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
  modelName       = "Swiss_VOT",
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

              b_tt_pt_business = 0,
              b_tt_pt_commuters = 0,
              b_tt_pt_leisure = 0,
              b_tt_pt_shopping = 0,

              b_tt_car_business = 0,
              b_tt_car_commuters = 0,
              b_tt_car_leisure = 0,
              b_tt_car_shopping = 0,

              b_dc = 0,
              
              b_tc_business = 0,
              b_tc_commuters = 0,
              b_tc_leisure = 0,
              b_tc_shopping = 0,

              b_ic_business = 0,
              b_ic_commuters = 0,
              b_ic_leisure = 0,
              b_ic_shopping = 0,

              b_hw_business = 0,
              b_hw_commuters = 0,
              b_hw_leisure = 0,
              b_hw_shopping = 0,

              b_bus_commuters = 0,
              b_bus_leisure = 0,
              b_bus_shopping = 0,

              b_rail_business = 0,
              b_rail_commuters = 0,
              b_rail_leisure = 0,
              b_rail_shopping = 0,

              b_car_inertia   = 0,
              b_car_available   = 0,
              b_car_male   = 0,

              b_bus_disc = 0,
              b_rail_disc = 0,
              b_bus_ga = 0,
              b_rail_ga = 0,

              lam_dist_tt_pt_business = 1,
              lam_dist_tt_pt_commuters = 1,
              lam_dist_tt_pt_leisure = 1,
              lam_dist_tt_pt_shopping = 1,

              lam_dist_tt_car_business = 1,
              lam_dist_tt_car_commuters = 1,
              lam_dist_tt_car_leisure = 1,
              lam_dist_tt_car_shopping = 1,

              lam_dist_tc = 1,
              lam_inc_tc_business = 1,
              lam_inc_tc_commuters = 1,

              mu_mc_car_bus = 1,
              mu_mc_car_rail = 1,
              mu_rc_bus = 1,
              mu_rc_car = 1,
              mu_rc_rail_by_car = 1,
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
  bus_p = 
    b_bus_commuters * PUR_P + 
    b_bus_leisure * PUR_T + 
    b_bus_shopping * PUR_E
  
  rail_p = 
    b_rail_business * PUR_N + 
    b_rail_commuters * PUR_P + 
    b_rail_leisure * PUR_T + 
    b_rail_shopping * PUR_E
  
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
     b_tc_business * PUR_N * (inc/mean(inc))^lam_inc_tc_business * (dist/mean(dist))^lam_dist_tc +
     b_tc_commuters * PUR_P * (inc/mean(inc))^lam_inc_tc_commuters * (dist/mean(dist))^lam_dist_tc +
     b_tc_leisure * PUR_T * (dist/mean(dist))^lam_dist_tc +
     b_tc_shopping * PUR_E * (dist/mean(dist))^lam_dist_tc
   
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
   
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["car"]]  = b_car_inertia * CAR   + b_car_available * VEH_AVA1 + b_car_male * MALE
  V[["bus"]]  = b_bus_disc  * DISCOUNT + b_bus_ga  * GA + bus_p
  V[["rail"]] = b_rail_disc * DISCOUNT + b_rail_ga  * GA + rail_p
  
  V[["a"]]    = asc_a + b_dc * S_O_C_A * CAR_A  +  tt_pt_elas * TT_A * PT_A  +  
    tt_car_elas * TT_A * CAR_A + tc_elas * TC_A + (ic_p + hw_p) * PT_A
  
  V[["b"]]    = asc_b + b_dc * S_O_C_B * CAR_B  +  tt_pt_elas * TT_B * PT_B  +  
    tt_car_elas * TT_B * CAR_B + tc_elas * TC_B + (ic_p + hw_p) * PT_B
  
  ### Compute probabilities for the RP part of the data using MNL model
  mnl_settings_1 = list(
    alternatives  = c(car=1, bus = 2), 
    avail         = list(car= 1, bus=1), 
    choiceVar     = CHOICE, 
    utilities     = list(car  = mu_mc_car_bus*(V[["car"]] + V[["a"]]),
                         bus  = mu_mc_car_bus*(V[["bus"]] + V[["b"]])),
    rows          = (N_O_S==1)
  )
  
  P[["1"]] = apollo_mnl(mnl_settings_1, functionality)
  
  mnl_settings_2 = list(
    alternatives  = c(car=1, rail = 2), 
    avail         = list(car= 1, rail=1), 
    choiceVar     = CHOICE, 
    utilities     = list(car  = mu_mc_car_rail*(V[["car"]] + V[["a"]]),
                         rail = mu_mc_car_rail*(V[["rail"]]+ V[["b"]])),
    rows          = (N_O_S==2)
  )
  
  P[["2"]] = apollo_mnl(mnl_settings_2, functionality)
  
  mnl_settings_3 = list(
    alternatives  = c(a=1, b = 2), 
    avail         = list(a = 1, b = 1), 
    choiceVar     = CHOICE, 
    utilities     = list(a  = mu_rc_bus *V[["a"]],
                         b  = mu_rc_bus *V[["b"]]),
    rows          = (N_O_S==3)
  )
  
  P[["3"]] = apollo_mnl(mnl_settings_3, functionality)
  
  mnl_settings_4 = list(
    alternatives  = c(a=1, b = 2), 
    avail         = list(a = 1, b = 1), 
    choiceVar     = CHOICE, 
    utilities     = list(a  = mu_rc_car *V[["a"]],
                         b  = mu_rc_car *V[["b"]]),
    rows          = (N_O_S==4)
  )
  
  P[["4"]] = apollo_mnl(mnl_settings_4, functionality)
  
  mnl_settings_5 = list(
    alternatives  = c(a=1, b = 2), 
    avail         = list(a = 1, b = 1), 
    choiceVar     = CHOICE, 
    utilities     = list(a  = mu_rc_rail_by_car *V[["a"]],
                         b  = mu_rc_rail_by_car *V[["b"]]),
    rows          = (N_O_S==5)
  )
  
  P[["5"]] = apollo_mnl(mnl_settings_5, functionality)
  
  mnl_settings_6 = list(
    alternatives  = c(a=1, b = 2), 
    avail         = list(a = 1, b = 1), 
    choiceVar     = CHOICE, 
    utilities     = list(a  = mu_rc_rail *V[["a"]],
                         b  = mu_rc_rail *V[["b"]]),
    rows          = (N_O_S==6)
  )
  
  P[["6"]] = apollo_mnl(mnl_settings_6, functionality)
  
  ### Combined model
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
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


