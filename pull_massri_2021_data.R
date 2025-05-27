#script to download, unpack, and read in data from Massri 2021

#make new directory for these data
dir.create(here::here("objects", "Massri_data"), showWarnings = FALSE) 

#extend timeout for large download
options(timeout = max(1200, getOption("timeout"))) #20 minute timeout (default is 1 minute)
#download from SRA (~1.8Gb)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184538&format=file&file=GSE184538%5Flv%5Fumap%2Erds%2Egz",
              destfile = here::here("objects", "Massri_data", "GSE184538_lv_umap.rds.gz"),
              mode = "wb")
#download from SRA (~790 Kb)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184538&format=file&file=GSE184538%5Flv%5Fumap%5Fmarkers%2Ecsv%2Egz",
              destfile = here::here("objects", "Massri_data", "GSE184538_lv_umap_markers.csv.gz"),
              mode = "wb")

#extract
R.utils::gunzip(here::here("objects", "Massri_data", "GSE184538_lv_umap.rds.gz"), remove = TRUE)
R.utils::gunzip(here::here("objects", "Massri_data", "GSE184538_lv_umap_markers.csv.gz"), remove = TRUE)

