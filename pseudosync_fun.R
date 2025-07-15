# (C) Irene Sophia Plank
# 
# Pseudosynchrony is calculated using the CCF function from rMEA. 
# 
# Two approaches can be used by setting the method input variable, both are 
# described in in Moulder et al. (2018, Psychol Methods) 
# 
# Choose "data" to test the following hypothesis: "There is no time 
# dependency at all between these two time series." (data shuffling). 
# This data shuffling, from Moulder et al., "For each pair of time series (X, Y) 
# each data point, Xi, is randomly shuffled to create a new time series, Xs 
# until no single data point existed at its original time point Xi!=Xjs.". 
# This new time series is then paired with an unchanged time series Y. 
#
# Choose "seg" to to test the following hypothesis: "Synchrony does not 
# exist between sections of size m in these two time series" (section sliding). 
# This segment shuffling, from Moulder et al., "requires researchers to cut a 
# time series X into shorter sections of size m which are randomly appended to 
# one another to create a new time series Xs until no section is in its original 
# position". This shuffled time series is then paired with an unchanged time 
# series Y. 
#
# Input: 
#     * method: string, either "data" or "seg"
#     * mea: a list of MEA or fakeMEA objects
#     * sampRate: sampling rate per second
#     * lagSec, winSec, incSec: length of lag, window and increment in seconds
#     * n: how many permutations to perform per MEA object (default = 100)
#     * peak: whether to apply peak picking or take grandaverage (default = F)
#     * r2Z: whether to apply z transformation (BOOLEAN, default T)
#     * ABS: whether to report absolute values (BOOLEAN, default T)
# Output:
#     * df.psync: data frame in long format containing pseudosync values
#

# Packages ----------------------------------------------------------------

# if packman is not installed yet, install it
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
pacman::p_load(tidyverse, rMEA) 

# Data shuffling ----------------------------------------------------------

pseudosync = function(method, mea, sampRate, lagSec, winSec, incSec, 
                      n = 100, peak = F, r2Z = T, ABS = T) {
  
  # create new data frame for pseudosync values
  cols = c("name", "psync")
  df.psync = data.frame(matrix(nrow = 0, ncol = length(cols)))
  colnames(df.psync) = cols
  
  # check if there are elements in the list
  if (length(mea) < 1) stop('There are no elements in the list provided!')

  # go through all list elements
  for (i in 1:length(mea)) {
    
    print(sprintf("%s: %i of %i", Sys.time(), i, length(mea)))
    
    if (method == "seg") {
      # divide left and right into segments
      mea_1 = split(mea[[i]][["MEA"]][[1]], floor(seq_along(mea[[i]][["MEA"]][[1]])/(winSec*sampRate)))
      mea_2 = split(mea[[i]][["MEA"]][[2]], floor(seq_along(mea[[i]][["MEA"]][[2]])/(winSec*sampRate)))
    }
    
    # shuffle each side n/2 times
    for (j in 1:(n/2)) {
      
      # shuffling first person
      if (method == "data") {
        mea[[i]][["MEA"]][,1] = sample(mea[[i]][["MEA"]][,1])
      } else if (method == "seg") {
        mea[[i]][["MEA"]][,1] = unlist(sample(mea_1))
        mea[[i]][["MEA"]][,2] = unlist(mea_2)
      } else {
        stop("Only data or seg can be used as method.")
      }
      # compute synchrony
      mea.ccf = MEAccf(mea[[i]], lagSec, winSec, incSec, r2Z, ABS)
      # extract peak or grandaverage
      if (peak) {
        psync = mean(apply(mea.ccf[["ccf"]], 1, max, na.rm = T ), na.rm = T)
      } else {
        psync = mea.ccf[["ccfRes"]][["grandAver"]]
      }
      # create new row in the data frame
      df.psync[nrow(df.psync) + 1,] = c(names(mea)[i], psync)      
      
      # shuffling second person
      if (method == "data") {
        mea[[i]][["MEA"]][,2] = sample(mea[[i]][["MEA"]][,2])
      } else if (method == "seg") {
        mea[[i]][["MEA"]][,1] = unlist(mea_1)
        mea[[i]][["MEA"]][,2] = unlist(sample(mea_2))
      } else {
        stop("Only data or seg can be used as method.")
      }
      # compute synchrony
      mea.ccf = MEAccf(mea[[i]], lagSec, winSec, incSec, r2Z, ABS)
      # compute pseudosync
      if (peak) {
        psync = mean(apply(mea.ccf[["ccf"]], 1, max, na.rm = T ), na.rm = T)
      } else {
        psync = mea.ccf[["ccfRes"]][["grandAver"]]
      }
      # create new row in the data frame
      df.psync[nrow(df.psync) + 1,] = c(names(mea)[i], psync)
    }
  }
  
  return(df.psync)
  
}