library(spatstat)
library(terra)
library(geodata)
library(maps)
library(rnaturalearth)
library(sf)
library(mgcv)
library(INLA)
library(fmesher)
library(geosphere)
library(here)   # should start from the directory that contains your subdirectories called "code", "data", and "output"

# load all functions 
source(here("code", "functions.R"), verbose=FALSE)

# run this script to generate all figures in note 2
source(here("code", "code_for_note2.R"), verbose=FALSE)

# run this script to generate all figures in discussion 1 of note 3 scale 1
source(here("code", "scale1_nonstructure_dd_spatial.R"), verbose=FALSE)

# run this script to generate all figures in discussion 2 of note 3 scale 1
source(here("code", "scale1_structure_dd_spatial.R"), verbose=FALSE)

# run this script to generate all figures in note 3 scale 2
source(here("code", "scale2.R"), verbose=FALSE)
