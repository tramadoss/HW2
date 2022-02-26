#  Set up Swiss VOT Data
library(haven)
library(expss)
swiss_vot_rdata = haven::read_spss("vot_swiss_F1_dsb0.sav")
# add missing 'labelled' class
swiss_vot_rdata = add_labelled_class(swiss_vot_rdata) 
#
saveRDS(swiss_vot_rdata, file = "swiss_vot_rdata.RDS")
rm(swiss_vot_rdata)
data_test <- readRDS("swiss_vot_rdata.RDS")
