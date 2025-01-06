### Distribution extraction
library(sf)
sf_use_s2(FALSE)
library(stars)
library(abind)
library(tidyverse)
library(terra)



{
  sdm_wd   <- "/Volumes/projects/bioing/user/slisovsk/ArcticSDM/SDM_Results/"
  unres_wd <- "/Volumes/projects/bioing/user/slisovsk/ArcticSDM/finalRun/Arrays_Unrestricted/"
  restr_wd <- "/Volumes/projects/bioing/user/slisovsk/ArcticSDM/finalRun/Arrays_Restricted/"
  out_wd   <- "/Volumes/projects/bioing/user/slisovsk/ArcticSDM/finalRun/"
}

spTable <- tibble(species = list.files(unres_wd)) %>%
              bind_rows(tibble(species = list.files(restr_wd))) %>%
              filter(!duplicated(species))

save(spTable, file = glue::glue("{out_wd}/spTable.rda"))

### Grid
load("/Volumes/projects/bioing/user/slisovsk/ArcticSDM/data/grid_25km.rda")
########

predArray <- array(dim = c(nrow(spTable), nrow(grid), 4, 3, 2))

for(sp in spTable$species) {
  
  ### unrestricted
  tryCatch(load(glue::glue("{unres_wd}/{sp}/binArray_Unconstraint.rda")), error = function(e) NULL)
  if(!is.null(binarArray)) unc <- binarArray else unc <- array(dim = c(1, nrow(grid), 4, 3))
  
  ### restricted
  tryCatch(load(glue::glue("{restr_wd}/{sp}/binArray_1.rda")), error = function(e) NULL)
  if(!is.null(binarArray)) con <- binarArray else con <- array(dim = c(1, nrow(grid), 4, 3))
  
  out <- abind::abind(unc, con, along = 5)
  
  predArray[which(spTable$species==sp),,,,] <- out
  
}

save(predArray, file = glue::glue("{out_wd}/predArray.rda"))
