# @Title: Preview Culture Sensitivities Data
# @Author: Steph Reynolds (Stephanie.Reynolds@ucsf.edu)
# @Project: Mindscape
# @DateCreated: 01-11-2022
# @DateModified: 

# Install packages 
# install.packages("AMR")
# install.packages("cleaner")

# Load required packages
library(tidyverse)
library(lubridate)
library(AMR)     # To analyze AMR data
library(cleaner)     # To create frequency tables for AMR data

# Read in culture sensitivities data and assign to `cs`
cs <- read_delim("/Users/sreynolds2/Downloads/ucsf_data_pull_11.08.21/Culture sensitivities 11.08.21.csv")

# Transform variables 
cs <- cs %>% 
  select(deid_enc_id, org_name, abx_name, resistance, sensitivity_value, analyzed_time) %>% 
  mutate(mo = as.mo(org_name),    # get microbial ID based on given organism 
         ab = as.ab(abx_name),    # get antibiotic name from  given antibiotic
         resistance = as.rsi(resistance),   # transform to 'rsi' class (ordered factor) 
         analyzed_time = ymd_hms(analyzed_time)) %>%    
  select(-c(org_name, abx_name)) %>% 
  distinct(deid_enc_id, mo, analyzed_time, ab, .keep_all = T)  # remove duplicates

# Save data as csv 
write_csv(cs, "/Users/sreynolds2/Documents/GitHub/MS-Antibiogram/data/clean/cs_01.11.22.csv")

# Frequency table of organism names (`mo`)
cs %>% freq(mo_name(mo), nmax = 20)

# Frequency table of antibiotic names (`ab`)
cs %>% freq(ab_name(ab), nmax= 20)

# Bar plot for AMR data anlalysis 
cs %>% 
  ggplot_rsi(combine_SI = F, datalabels = F)   # do not combine 'S' and 'I'

