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

# Read in culture sensitivities data and assign to `cs_raw`
cs_raw <- read_delim("/Users/sreynolds2/Downloads/ucsf_data_pull_11.08.21/Culture sensitivities 11.08.21.csv")

# Remove unnecessary variables
cs <- cs_raw %>% select(deid_enc_id, order_name, org_name, abx_name, resistance, analyzed_time)
cs_raw <- NULL

# Summarize counts of all abx 
all_abx <- cs %>%
  group_by(abx_name) %>% 
  count() %>% 
  arrange(desc(n), abx_name)
all_abx

# With count data, create list of most common abx (those with n>2500)
top_abx <- all_abx %>% 
  filter(n>2500)
top_abx

# Summarize counts of all order names
all_order_names <- cs %>% 
  group_by(order_name) %>% 
  count() %>% 
  arrange(desc(n), order_name)

# Summarize counts of all organisms/isolates
all_isolates <- cs %>% 
  group_by(org_name) %>% 
  count() %>% 
  arrange(desc(n), org_name)

# Read in `microorganisms` dataset to obtain full organism names; assign to `microorg`
microorg <- microorganisms %>% 
  select(mo, fullname) %>%              # only keep `mo` and `fullname` variables
  mutate(mo = as.factor(mo))            # change `mo` to factor class

# Tranform `org_name` and `abx_name` to factor class 
cs$org_name <- as.factor(cs$org_name)
cs$abx_name <- as.factor(cs$abx_name)

# Transform `resistance` to rsi class (ordinal factor)
cs$resistance <- as.rsi(cs$resistance)

# Count number of missing `resistance` values 
sum(is.na(cs[, "resistance"]))

# Create table of frequency counts of unique organisms 
num_orgs <- cs %>% 
  filter(!is.na(resistance)) %>%                        # remove rows where `resistance` is NA 
  distinct(deid_enc_id, org_name, .keep_all = TRUE) %>% # table of unique isolates per patient ID 
  group_by(org_name) %>%                                # group by organism 
  summarize(n_isolates = n()) %>%                       # count num of organisms/isolates
  filter(n_isolates>=30) %>%                            # remove rows where `n_isolates`<30 
  arrange(desc(n_isolates)) %>%                         # sort by `n_iolates` in descending order 
  inner_join(microorg, by = c('org_name' = 'fullname')) # join to `microorg` table to link with `fullname`

n_iso

sum(n_iso$n_isolates)                               # 7,278 total organisms/isolates 

# Create list of names included in n_iso$mo
mo_names <- list(n_iso[1:13, 1])

abx_tbl_grp_by_mo <- cs %>% 
  filter(!is.na(resistance)) %>% 
  filter(mo %in% c(mo_names)) %>% 
  select(-c(deid_enc_id, analyzed_time, sensitivity_value)) %>%
  group_by(mo, ab) %>% 
  summarize(n(resistance=='S'))

freq_mo
freq_mo_unique <- n_iso %>% freq(mo_name(mo))
freq_mo_unique

# Frequency table of organism names (`mo`)
freq_mo <- cs %>% freq(mo_name(mo), nmax = 20, title = "Frequency Counts of Organisms\nUCSF Hospital Data from 11/2/19 to 11/2/21")

# Create list of top 10 organism names 
mo_names_10 <- as.mo(freq_mo[1:10,1])







cs <- cs_raw %>% 
  select(deid_enc_id, org_name, abx_name, analyzed_time, resistance) %>% 
  mutate(mo = as.mo(org_name),    # get microbial ID based on given organism 
         ab = as.ab(abx_name),    # get antibiotic name from  given antibiotic
         resistance = as.rsi(resistance),   # transform to 'rsi' class (ordered factor) 
         analyzed_time = ymd_hms(analyzed_time)) %>%    
  select(-c(org_name, abx_name, analyzed_time)) %>% 
  filter(!is.na()) %>% 
  distinct(deid_enc_id, mo, analyzed_time, ab, .keep_all = T)  # remove duplicates


sum(is.na(cs[,4]))
cs <- cs %>% mutate(mo = as.mo(org_name),    # get microbial ID based on given organism 
                ab = as.ab(abx_name),    # get antibiotic name from  given antibiotic
                resistance = as.rsi(resistance),   # transform to 'rsi' class (ordered factor) 
                analyzed_time = ymd_hms(analyzed_time)) %>% 
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

# Create table showing frequency counts for each organism/isolate 

microorganisms <- microorganisms %>% select(mo, fullname) %>% mutate(mo = as.factor(mo))

cs$mo <- as.factor(cs$mo)
levels(cs$mo)

n_iso <- cs %>% 
  filter(!is.na(resistance)) %>%                    # remove rows where `resistance` is NA 
  distinct(deid_enc_id,mo, .keep_all = TRUE) %>%    # table of unique isolates per patient ID 
  group_by(mo) %>%                                  # group by organism 
  summarize(n_isolates = n()) %>%                   # count num of organisms/isolates
  filter(n_isolates>=30) %>%                        # remove rows where `n_isolates`<30 
  arrange(desc(n_isolates)) %>%                     # sort by `n_iolates` in descending order 
  inner_join(microorganisms, by = 'mo')             # join to `microorganisms` table to link with `fullname``

n_iso

sum(n_iso$n_isolates)                               # 7,278 total organisms/isolates 

# Create list of names included in n_iso$mo
mo_names <- list(n_iso[1:13, 1])

abx_tbl_grp_by_mo <- cs %>% 
  filter(!is.na(resistance)) %>% 
  filter(mo %in% c(mo_names)) %>% 
  select(-c(deid_enc_id, analyzed_time, sensitivity_value)) %>%
  group_by(mo, ab) %>% 
  summarize(n(resistance=='S'))

freq_mo
freq_mo_unique <- n_iso %>% freq(mo_name(mo))
freq_mo_unique

# Frequency table of organism names (`mo`)
freq_mo <- cs %>% freq(mo_name(mo), nmax = 20, title = "Frequency Counts of Organisms\nUCSF Hospital Data from 11/2/19 to 11/2/21")

# Create list of top 10 organism names 
mo_names_10 <- as.mo(freq_mo[1:10,1])




library(tableone)

#CreateCatTable(data=cs, vars=c("mo", "resistance"), strata="ab")

longtbl <- cs %>%
  group_by(mo, ab) %>% 
  #filter(resistance == "S") %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = ab, values_from = n, values_fill = 0) 

all <- cs %>% 
  filter(resistance=="S") %>%
  group_by(mo, ab) %>% 
  summarize(Tota_Num_Isolates = n()) %>% 
  arrange(desc(Tota_Num_Isolates))


s <- cs %>% 
  filter(resistance=="S") %>%
  group_by(mo, ab) %>% 
  summarize(Num_Susc_Isolates = n()) %>% 
  arrange(desc(Num_Susc_Isolates))

longtbl <- cs %>%
  group_by(mo, ab) %>% 
  summarise(sum = sum()) %>% 
  pivot_wider(names_from = ab, values_from = sum, values_fill = 0)

