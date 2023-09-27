library(MCMCglmm)
library(tidyverse)
library(naniar)


#  expect a single argument from Rscript which is the array_idx
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 1) {
  stop("Must be called from Rscript with array_idx as the sole argument after the name of the script.")
} else {
  array_idx = as.integer(args[1])
}


Ped_length_age <- read_csv(
  "RR_2007-2020_pedigree_ages_lengths.csv",
  col_names = T,
  cols(
    kid = col_character(),
    ma = col_character(),
    pa = col_character(),
    kid_dayofyear = col_character(),
    parent_dayofyear = col_character(),
    parent_sg = col_character(),
    pa_sg = col_character(),
    ma_sg = col_character(),
    kid_sex = col_character(),
    new_kid_sex = col_character(),
    pa_sex = col_character(),
    new_pa_sex = col_character(),
    ma_sex = col_character(),
    new_ma_sex = col_character(),
    SpawnYear = col_double(),
    kid_year = col_double(),
    kid_year2 = col_double(),
    pa_year = col_double(),
    ma_year = col_double(),
    kid_sg = col_character(),
    kid_hatchery = col_character(),
    pa_hatchery = col_character(),
    ma_hatchery = col_character(),
    OffspCollection = col_character(),
    PopName = col_character(),
    FDR = col_double(),
    Pvalue = col_double(),
    LOD = col_double(),
    P.Pr.C_Se_Se = col_double(),
    P.Pr.Max = col_double(),
    MaxP.Pr.Relat = col_character(),
    TotPaNonExc = col_double(),
    TotMaNonExc = col_double(),
    TotUnkNonExc = col_double(),
    TotPairsMendCompat = col_double(),
    TotPairsMendAndLogL = col_double(),
    TotParsMendLoglAndRank = col_double(),
    TotPairsNonExc = col_double(),
    KidMiss = col_double(),
    PaMiss = col_double(),
    MaMiss = col_double(),
    MI.Kid.Pa = col_double(),
    MI.Kid.Ma = col_double(),
    MI.Trio = col_double(),
    MendIncLoci = col_character(),
    kid_age = col_double(),
    new_parent_sg = col_date(format = ""),
    new_kid_sg = col_date(format = ""),
    pa_age = col_double(),
    ma_age = col_double(),
    kid_length = col_double(),
    ma_length = col_double(),
    pa_length = col_double()
  )
)

Known_Sexes_fish_list <- read_csv(
  "RR_steelhead_07-20_known_sexes.csv",
  col_names = T,
  cols(
    indiv = col_character(),
    year = col_character(),
    spawner_group = col_character(),
    hatchery = col_character(),
    length = col_double(),
    watershed = col_character(),
    Wild = col_character(),
    sex = col_character(),
    new_sex = col_character()
  )
)
Known_Sexes_fish_list <- Known_Sexes_fish_list %>%
  replace_with_na(replace = list(length = 203.2))


#remove the parents from 2018 that could only have 2 yr old offspring in the dataset
Pedigree.prep <- Ped_length_age %>%
  filter(new_ma_sex == "Female") %>%
  filter(new_pa_sex == "Male") %>%
  select(kid_year, kid, pa, ma, SpawnYear) %>%
  filter(SpawnYear < 2018) %>%
  select(-SpawnYear) %>%
  arrange(kid_year)

ped.stuff <- Known_Sexes_fish_list %>%
  select(indiv) %>%
  rename(animal = indiv)
ped.stuff2 <- Pedigree.prep %>%
  select(kid, ma, pa) %>%
  rename(animal = kid)
ped.stuff2$pa <- recode(
  ped.stuff2$pa,
  M035631_Coyote_Valley = "M035631",
  M047366_Warm_Springs = "M047366",
  M065002_Coyote_Valley = "M065002"
)

pedigree <- full_join(ped.stuff, ped.stuff2, by = "animal") %>%
  arrange(!is.na(pa), animal) %>%
  data.frame()


RR_w_Omy <- read_csv("RR_2007-2020_Omy5_meta_pedigree.csv")

OmyModel <- RR_w_Omy %>%
  filter(kid_SH114448_types != "NA_NA") %>%
  filter(!is.na(new_kid_sex)) %>%
  filter(!is.na(new_ma_sex)) %>%
  filter(!is.na(new_pa_sex)) %>%
  filter(kid_age < 5) %>%
  mutate(age_binary = ifelse(kid_age == 2, yes = 0, no = 1))

OmyModel$pa <- recode(
  OmyModel$pa,
  M035631_Coyote_Valley = "M035631",
  M047366_Warm_Springs = "M047366",
  M065002_Coyote_Valley = "M065002"
)

Pheno <- OmyModel %>%
  select(kid, ma, new_kid_sex, age_binary, ma_hatchery, SpawnYear, kid_SH114448_types) %>%
  rename(animal = kid) %>%
  data.frame()
Pheno$SpawnYear <- as_factor(Pheno$SpawnYear)

Ped <- prunePed(pedigree, Pheno$animal, make.base = T)

prior <- list(R=list(V = diag(2), nu = 1.002),
              G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                     G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                     G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

outdir <- "outputs/genetic_corr"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile_path <- paste0(outdir, "/RR_genetic_corr_", sprintf("%02d", array_idx), ".rds")

# make sure to set the seed differently for each run
set.seed(array_idx)

#changed the hatchery to ma_hatchery so the fish are grouped by hatchery of origin
Model.binary.age.correlation <- MCMCglmm(
  age_binary ~ new_kid_sex + ma_hatchery + kid_SH114448_types,
  random = ~ us(new_kid_sex):animal + ma + SpawnYear,
  rcov = ~ idh(new_kid_sex):units,
  data = Pheno,
  pedigree = Ped,
  family = "threshold",
  prior = prior,
  slice = TRUE,
  verbose = TRUE,
  nitt = 775000,
  burnin = 25000,
  thin = 500
)


write_rds(Model.binary.age.correlation, path = outfile_path, compress = "xz")

