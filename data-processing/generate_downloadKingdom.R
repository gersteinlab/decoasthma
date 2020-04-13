#
# generate_downloadKingdom.R
# 
# Script to generate a file with a series of wget commands to download the most 
# recent release of ensemble genomes.
# 
# Created 04 Jun 2017 by DS

library(RCurl)
library(tidyverse)

kingdomsOfInterest <- c("fungi", "bacteria", "plants", "metazoa", "protists")

# Get the directories at the ensemble ftp
url <- list()
for (i in kingdomsOfInterest) {
  url[[i]] <- paste("ftp://ftp.ensemblgenomes.org/pub/", i, "/", sep = "")
}

# Retrieve all urls
dirs <- lapply(url, getURL)

# Split to individual lines
dirnames <- dirs %>%
  lapply(function(x) strsplit(x, split = "\n")) %>%
  lapply(unlist)

# Retrieve the current release
currenturl <- list()
for (i in names(dirnames)) {
  currentRelease <- dirnames[[i]][grep("current", dirnames[[i]])] %>%
    strsplit(split = "-> ") %>%
    unlist
  currenturl[[i]] <- paste(url[[i]], currentRelease[2], "/fasta/", sep = "")
}

# Get the directories within the current release
dirscurr <- currenturl %>%
  lapply(getURL) %>%
  lapply(function(x) strsplit(x, split = "\n")) %>%
  lapply(unlist)

# Get vector of organisms
orgs <- dirscurr %>%
  lapply(function(x) 
    lapply(x, function(y) strsplit(y, split = " "))) %>%
  lapply(function(x) 
    lapply(x, function(y) unlist(y))) %>%
  lapply(function(x) 
    sapply(x, function(y) y[length(y)]))


# Special case: collections
# Find the collections within the list of orgs
collections <- orgs %>%
  lapply(function(x) grep("collection", x))

# Subset to just collections and non-collections
collec_dirs <- list()
for (i in names(orgs)) {
  collec_dirs[[i]] <- orgs[[i]][collections[[i]]]
}

org_dirs <- list()
for (i in names(orgs)) {
  if (any(collections[[i]])){
    org_dirs[[i]] <- orgs[[i]][-collections[[i]]]
  } else {
    org_dirs[[i]] <- orgs[[i]]
  }
  if (length(orgs[[i]] == length(collections[[i]]))) {
    orgs[[i]] <- NULL
  }
}

# urls of each collection
col_url <- list()
for (i in names(collec_dirs)) {
  if (length(collec_dirs[[i]] > 0)) {
    for (j in 1:length(collec_dirs[[i]])) {
      col_url[[i]][j] <- paste(currenturl[[i]], collec_dirs[[i]][j], "/", sep = "")
    }
  } else {
    col_url[[i]] <- NULL
  }
}

# Get organism names for each collection
# This was erroring out on my machine so I moved it to grace and then brought it back in
save(col_url, file = "col_url")
load("~/Documents/projects_local/chas/candida_alignments/col_url_str.RData")

# col_url_str <- list()
# for (i in names(col_url)){
#   for (j in 1:length(col_url[[i]])){
#     col_url_str[[i]][j] <- getURL(col_url[[i]][j], ssl.verifypeer = FALSE,
#                                   dirlistonly = TRUE, verbose = TRUE)
#   }
# }

# Take the websites and parse to just names of orgs
# Keep within lists to follow the collection directory structure
col_url_names <- col_url_str %>%
  lapply(function(x)
    lapply(x, function(y) strsplit(y, split = "\n"))) %>%
  lapply(function(x) 
    lapply(x, unlist))

# Write a wget statement for each organism NOT in a collection
notcol_wget <- vector(mode = "list", length = length(kingdomsOfInterest))
names(notcol_wget) <- kingdomsOfInterest
for (n in names(org_dirs)) {
  if (length(org_dirs[[n]]) > 0) {
    for (i in 1:length(org_dirs[[n]])) {
      notcol_wget[[n]][i] <- paste("wget -qO- ", currenturl[[n]], 
                                   org_dirs[[n]][i], 
                                   '/dna/*.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>', 
                                   capitalize(n), ':', org_dirs[[n]][i],
                                   ':/" >> ', capitalize(n), '.fa', sep = "")
    } 
  } 
}

# Write a wget statement for each organism in the collection
# Preallocate memory
collections_wget <- c(col_url_names, 
                      kingdomsOfInterest[!(kingdomsOfInterest %in% names(col_url_names))])
for (n in names(col_url_names)) {
  for (i in 1:length(col_url_names[[n]])) {
    for (j in 1:length(col_url_names[[n]][[i]])) {
      collections_wget[[n]][[i]][j] <- paste("wget -qO- ", col_url[[n]][i], 
                                             col_url_names[[n]][[i]][j], 
                                             '/dna/*.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>', 
                                             capitalize(n), ':', col_url_names[[n]][[i]][j],
                                             ':/" >> ', capitalize(n), '.fa', sep = "")
    }
  }
}

# Combine all wget lines
all_wget <- vector(mode = "list", length = length(kingdomsOfInterest))
names(all_wget) <- kingdomsOfInterest
for (n in names(notcol_wget)) {
  all_wget[[n]] <- c(notcol_wget[[n]], collections_wget[[n]])  
}

# Unlist for writing
all_wget <- lapply(all_wget, unlist)

# Write output to file
for (i in names(all_wget)) {
  fileOut <- file(paste("download", capitalize(i), ".sh", sep = ""))
  writeLines(c("#!/bin/bash", 
               "#SBATCH --ntasks=1 --nodes=1",
               "#SBATCH --time=12:00:00",
               paste("#SBATCH --job-name=download", capitalize(i), sep = ""),
               "#SBATCH --mail-user=daniel.spakowicz@yale.edu",
               "#SBATCH --mail-type=ALL",
               "# ",
               "# Generated by generate_downloadKingdom.R",
               "#  https://github.com/dspak/chas/generate_downloadKingdom.R",
               "# ",
               paste("# ", Sys.time(), sep = ""),
               "",
               paste("cd /project/fas/gerstein/djs88/references/", i, "/", sep = ""),
               "",
               all_wget[[i]]),
             fileOut)
  close(fileOut)
}
