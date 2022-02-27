library(tidyverse)
data_test <- readRDS("swiss_vot_rdata.RDS")

t1_general <- data_test %>% 
  with_groups(P_NR, mutate, rank = row_number())%>%
  filter(rank == 1) %>%
  select(DISCOUNT, GA, VEH_AVAI, EDUC, WORKING, HH_INC_A) %>%
  mutate(across(everything(), as.numeric)) %>%
  na.omit() 

t1_general %>%
  count(eval(parse(text = "EDUC"))) %>%
  mutate(perc = 100*`n`/sum(`n`))

sample_prop <- function(var, df) {
  new_df <- df %>%
    count(eval(parse(text = var))) %>%
    mutate(Percent = 100*`n`/sum(`n`)) %>%
    mutate(Variable = var, .before = 1) %>%
    select(-n)
    
    colnames(new_df) <- c("Variable", "Categories", "Percent")
    return(new_df)
}

socio <- c("DISCOUNT", "GA", "VEH_AVAI", "EDUC", "WORKING", "HH_INC_A")

t1 <- map_df(socio, function(x) sample_prop(x, t1_general)) %>%
  filter(!(Categories == 0 & Variable != "WORKING"))

## Table 2
t2 <- data_test %>% select(P_NR,contains("DIST"), contains("PUR")) %>% 
  mutate(across(everything(), as.numeric)) %>%
  with_groups(P_NR, mutate, rank = row_number())%>%
  filter(rank == 1) %>% 
  mutate(Distance = cut(T_DIST_C, 
                        c(-Inf, 5, 10, 20, 30, 50, 75, 100, Inf), 
                        labels = c("<5","5-10","10-20","20-30",
                                   "30-50","50-75","75-100",">100"))) %>%
  count(Distance, T_PUR) %>%
  inner_join(
    tibble(
      T_PUR = 1:4,
      Purpose = c("Commuting", "Shopping", "Business", "Leisure")
    )
  ) %>%
  pivot_wider(
    id_cols = Distance,
    names_from = Purpose,
    values_from = n
  ) %>%
  mutate(Total = rowSums(select_if(., is.numeric))) %>%
  mutate(across(-Distance, ~ 100* .x / sum(.x))) 


## Table 3
t3_counts <- data_test %>% 
  mutate(across(everything(), as.numeric)) %>%
  inner_join(
    tibble(
      T_PUR = 1:4,
      Purpose = c("Commuting", "Shopping", "Business", "Leisure")
    )
  ) %>%
  count(Purpose, N_O_S) %>%
  pivot_wider(
    id_cols = N_O_S,
    names_from = Purpose,
    values_from = n
  ) %>%
  select(-N_O_S) %>%
  rbind(colSums(.))

rownames(t3_counts) <- c("Mode choice: car vs. bus", "Mode choice: car vs. rail", 
                         "Route choice: bus for bus users", "Route choice: car for car users",
                         "Route choice: rail for car users", "Route choice: rail for rail users",
                         "Total")

t3_avg <- data_test %>% 
  select(T_PUR, HH_INC_A, T_DIST_C) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(T_DIST_C = if_else(T_DIST_C == 0, 1 , T_DIST_C)) %>%
  inner_join(
    tibble(
      T_PUR = 1:4,
      Purpose = c("Commuting", "Shopping", "Business", "Leisure")
    )
  ) %>%
  select(-T_PUR) %>%
  rbind(cbind(summarise(.,across(-Purpose, mean)), tibble(Purpose = "Total"))) %>%
  with_groups(Purpose, summarise, across(everything(),mean, na.rm = TRUE)) 
  
t3_totals <- data_test %>% 
  with_groups(P_NR, mutate, rank = row_number()) %>%
  filter(rank == 1) %>%
  select(T_PUR) %>%
  mutate(across(everything(), as.numeric)) %>%
  inner_join(
    tibble(
      T_PUR = 1:4,
      Purpose = c("Commuting", "Shopping", "Business", "Leisure")
    )
  ) %>%
  count(Purpose) %>% 
  rbind(cbind(summarise(.,across(-Purpose, sum)), tibble(Purpose = "Total"))) %>%
  rename(Respondents = n)

t3_addl <- inner_join(t3_totals, t3_avg) #%>%

colnames(t3_addl) <- c("Trip Purpose", "Total", "Avg Income", "Avg Trip Distance")

