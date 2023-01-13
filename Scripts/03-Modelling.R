
### load library

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(modleR))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(progressr))


### define path

project_data_ <- '/home/shpc_100785/seaweed_area/'

env_ <- paste0(project_data_, 'Data/ENV_data/current.grd')
occ_all_ <- paste0(project_data_, 'Data/Occurance_data')
output_folder <- paste0(project_data_, 'outputs_083_400_20221122')

### load raw occurance data

load_raw_occ_csv <- function(x){
  df <- read_csv(x, col_types = cols(), show_col_types = F)
  df <- df |> 
    dplyr::select(dis_name, longitude, latitude) |> 
    rename(sp = dis_name,
           lon = longitude,
           lat = latitude)
  return(df |> as.data.frame())
}

get_names <- function(x){
  df <- read_csv(x, show_col_types = F)
  dis_name <- unique(df$dis_name)
  
  return(dis_name)
}

occ_csv_ <- list.files(occ_all_, '.csv', full.names = T)

occs <- purrr::map(occ_csv_, load_raw_occ_csv)
species <- unlist(purrr::map(occ_csv_, get_names))

names(occs) <- species

### Loading Environmental data

env <- raster::stack(env_)

### Create buffer / Pre-processing data for input

plan(multisession, workers = 24) # on server

y <- as.list(names(occs))

with_progress({
  p <- progressor(steps = length(occs))
  
  nouse <- future_map2(.x = occs,
                       .y = y,
                       ~ {
                         p()
                         setup_sdmdata(
                           species_name = .y,
                           occurrences = .x,
                           partition_type = 'crossvalidation',
                           cv_n = 10,
                           cv_partitions = 10,
                           clean_nas = T,
                           clean_dupl = T,
                           clean_uni = T,
                           buffer_type = 'distance',
                           dist_buf = 4,
                           predictors = env,
                           models_dir = output_folder,
                           n_back = 10000
                         )
                       },
                       .options = furrr_options(seed = 20220724))
})

### Modelling

plan(multisession, workers = 24) # on server, longtime run > 120h

with_progress({
  p <- progressor(length(y))
  
  furrr::future_map(y, 
                    ~{
                      options( java.parameters = c("-Xss2560k", "-Xmx100g") )
                      do_many(species_name = .,
                              predictors = env,
                              models_dir = output_folder,
                              bioclim = T,
                              maxent = T,
                              maxnet = T,
                              # rf = T,
                              # svmk = T,
                              #svme = T,
                              # domain = F,
                              glm = T,
                              # mahal = T,
                              equalize = T,
                              png_partitions = T,
                              write_bin_cut = F,
                              write_rda = T)
                      
                      p()
                    })
})


### Combine sub-models for each species to a final model

plan(multisession, .cleanup = T, workers= 24)
with_progress({
  p <- progressor(steps = length(species))
  
  species %>%
    as.list(.) %>%
    furrr::future_map(
      ~ {final_model(
        species_name = .,
        consensus_level = 0.5,
        models_dir =  output_folder,
        which_models = c("raw_mean",
                         "bin_mean",
                         "bin_consensus"),
        overwrite = TRUE,
        uncertainty = T
      )
        p()
        }
    )
})
  
### Ensemble all species of all models

#### consensus-all 4

plan(multisession, .cleanup = T, workers= 24)
with_progress({
  p <- progressor(steps = length(species))
  
  occs %>% purrr::map2(
    .x = .,
    .y = as.list(names(.)),
    ~ {ensemble_model(
      species_name = .y,
      occurrences = .x,
      which_ensemble = "consensus",
      png_ensemble = TRUE,
      models_dir = output_folder,
      overwrite = F,
      consensus_level = .5,
      ensemble_dir = 'all_ensemble_05'
    )
      p()
      }
  )
})


