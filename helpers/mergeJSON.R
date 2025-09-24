# Script to combine OpenPose output (JSON files) to one CSV per video. 
# Created losely based on 1Merge_JSON.R, a file shared by Ken Fujiwara on OSF: https://osf.io/ueg7j

library(jsonlite)
library(zoo)
library(tidyverse)

# path to the JSON folders
task = 'MG'
dir.path = paste0("C:/Users/NEVIA/Desktop/ISP_MSYNC/", task, "_SC_out")
folders  = list.dirs(path = dir.path,
                     recursive = F)

# some settings about the videos
fps      = 120
skip     = 10
frame_width = 1920

StartFrame = fps*skip
StartFrame = StartFrame - fps/2 # for rollmean

if (task == "CT") {
  duration = 10*60
  EndFrame   = StartFrame + duration*fps
  EndFrame   = EndFrame   + fps/2 # for rollmean
}

# create column names
ls.cols = c(paste0(rep(paste0("L", 1:25), each = 3), c("x", "y", "c")), 
            paste0(rep(paste0("R", 1:25), each = 3), c("x", "y", "c")))

# number of dyads
nod = length(folders)

for (f in 1:nod){
  
  files = list.files(pattern = "\\.json$", path = folders[f])
  
  if (file.exists(file.path(dir.path, 
                            paste0(gsub('.*out/', "\\1", folders[f]),
                                   '.csv')))) {
    next
  }
  
  if (task == 'MG') {
    EndFrame = length(files)
  }
  
  print(folders[f])
  
  merged.matrix = c()
  ls.frame = c()
  perc = -1
  l = (EndFrame - StartFrame)/10
  for (i in StartFrame:EndFrame){
    # print if new 10 percent
    if (floor(perc) != floor(i/l)) {
      print(sprintf('%s : %.0f%%', Sys.time(), i*10/l))
    }
    # update perc
    perc = i/l
    filepath = file.path(folders[f], files[i])
    if (!file.exists(filepath)) {
      break
    }
    tmp    = fromJSON(filepath) # read JSON file
    tmp.P1 = unlist(tmp$people$pose_keypoints_2d[1]) # read person 1
    tmp.P2 = unlist(tmp$people$pose_keypoints_2d[2]) # read person 2 
    if (is.null(tmp.P1)) tmp.P1 = c(rep(-1,75))
    if (is.null(tmp.P2)) tmp.P2 = c(rep(-1,75))
    # who is left and who is right?
    if ((tmp.P1[1] + tmp.P1[4])/2 > (frame_width/2)) { # half of the screen width
      tmp.right = tmp.P1
      tmp.left  = tmp.P2
    } else {
      tmp.right = tmp.P2
      tmp.left  = tmp.P1
    }
    merged.matrix = rbind(merged.matrix, c(tmp.left, tmp.right))
    ls.frame = c(ls.frame, gsub(".*cut_(.+)_key.*", "\\1", files[i]))
  }
  df = as.data.frame(merged.matrix) %>%
    mutate(frame = ls.frame)
  colnames(df) = ls.cols

  write_csv(df, 
            file = file.path(dir.path, 
                             paste0(gsub('.*out/', "\\1", folders[f]),
                                    '.csv'))
            )
}