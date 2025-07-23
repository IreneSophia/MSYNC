ls.packages = c("knitr",# kable
                "tidyverse",        # tibble stuff, ggplot
                "ggpubr",           # ggarrange
                "BayesFactor",      # anovaBF
                "effectsize",       # interpretBF
                "ggrain",           # geom_rain
                "rMEA"
)

lapply(ls.packages, library, character.only=TRUE)

# set number of permutations
n = 2

# load pseudosync function
source("pseudosync_fun.R")

# Initialize function to create a fake MEA object out of two vectors
# Input: 
#     * s1, s2: numeric vectors containing the values to be correlated
#     * sampRate: sampling rate per second
#     * s1Name, s2Name: name for the values to be correlated, default is "s1Name" and "s2Name"
# Output:
#     * fake MEA object that pretends to be a MEA object
#
fakeMEA = function(s1, s2, sampRate, s1Name = "s1Name", s2Name = "s2Name") {
  mea = structure(list(all_01_01 = structure(list(MEA = structure(list(
    s1Name = s1, s2Name = s2), row.names = c(NA, -length(s1)), class = "data.frame"), 
    ccf = NULL, ccfRes = NULL), id = "01", session = "01", group = "all", sampRate = sampRate, 
    filter = "raw", ccf = "", s1Name = s1Name, s2Name = s2Name, uid = "all_01_01", 
    class = c("MEA","list"))), class = "MEAlist", nId = 1L, n = 1L, groups = "all", sampRate = sampRate, 
    filter = "raw", s1Name = s1Name, s2Name = s2Name, ccf = "")
  return(mea)
}

# set the path to the data folders
dt.path  = "/media/emba/emba-2/MSYNC/data"

# find out where we are to set a new seed
fls = list.files(path = file.path(dt.path, "pseudosyncOP"), 
                 pattern = "df_psync_MG-seg")
i = length(fls) + 1
write_csv(as.data.frame(0), file.path(dt.path, "pseudosyncOP", 
                            sprintf("df_psync_MG-seg_%02d.csv", i)))


# start with the conversation task 
attach(file.path(dt.path, "preprocessedOP", "MSYNC_OP_MG.Rdata"), name = 'MG')
meth = "seg"  
df.mea = pseudosync(meth, 
                    mea.ccf,
                    sampRate = 120, 
                    lagSec = 2,
                    winSec = 30, 
                    incSec = 1, 
                    n = n,
                    log = '/home/emba/pCloudDrive/data/NEVIA/output/LOGFILES/log_pseudo.txt',
                    seed = i)
# detach the CT data
detach('MG')

# save the pseudosynchrony values
write_csv(df.mea, file.path(dt.path, "pseudosyncOP", 
                            sprintf("df_psync_MG-seg_%02d.csv", i)))
