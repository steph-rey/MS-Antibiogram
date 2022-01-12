# @Title: Preview Culture Sensitivities Data
# @Author: Steph Reynolds (Stephanie.Reynolds@ucsf.edu)
# @Project: Mindscape
# @DateCreated: 01-11-2022
# @DateModified: 01-11-2022 at 12:55PM PT

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

# Save data as csv file
# write_csv(cs, "/Users/sreynolds2/Documents/GitHub/MS-Antibiogram/data/clean/cs_01.11.22.csv")

# Summary stats
summary(cs)

sum(is.na(cs$resistance))
# missing resistance (S/I/R) data for n=7398 rows 

# Frequency table of organism names (`mo`)
freq_mo <- cs %>% freq(mo_name(mo), nmax = 20, title = "Frequency Counts of Organisms\nUCSF Hospital Data from 11/2/19 to 11/2/21")

# Create list of top 10 organism names 
mo_names_10 <- as.mo(freq_mo[1:10,1])

# Filter dataset to include include top 10 organism, assign to `sample`
sample <- cs %>% filter(mo %in% mo_names_10)
sample

# Frequency table of antibiotic names (`ab`)
freq_ab <- cs %>% freq(ab_name(ab), nmax= 20, title = "Frequency Counts of Antibiotics\nUCSF Hospital Data from 11/2/19 to 11/2/21")

# Create list of antibiotic names 
ab_names <- as.list(freq_ab[[1]])

# Bar plot showing proportion of microbes that are S/I/R (susceptible, intermediate, resistant)
sample %>% ggplot_rsi(combine_SI = F, datalabels = F, facet = 'mo')  # for top 10 organisms
cs %>% ggplot_rsi(combine_SI = F, datalabels = F, facet = 'mo')  # for all organisms

cs %>% group_by(mo, ab) %>% summarize(count = n()) %>% arrange(desc(count)) 

library(tableone)

CreateCatTable(data=cs, vars=c("mo", "resistance"), strata="ab")

longtbl <- cs %>%
  group_by(mo, ab) %>% 
  #filter(resistance == "S") %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = ab, values_from = n, values_fill = 0)



