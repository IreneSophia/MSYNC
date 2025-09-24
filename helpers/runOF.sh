#!/bin/bash

# Bash script to run OpenFace on Ubuntu

# set our directories
dirout="preprocessedOF"

# list of files
files=(vids4OF/*)

# set up parallisation
cores=12
starts=($(seq 0 $cores $((nos-1))))

# start with first batch
for i in "${starts[@]}"
do
   : 

   # assess the last index in this batch
   if [[ $((i+cores)) > $nos ]]; then
      end=$((nos-1))
   else
      end=$((i+cores-1))
   fi
   
   # loop through all subjects in this batch
   sequ=($(seq $i 1 $end))
   echo "$(date) starting: $i to $end" >> "$dirout/logfile.txt"
   for j in "${sequ[@]}"
   do
      :
      video=${files[$j]}
      base="$(basename -- $video .mp4)"
      file="$dirout/${base}_of_details.txt"
      if [ ! -f "$file" ]; then
	      echo "$(date) starting: $base" >> "$dirout/logfile.txt"

	      # run OpenFace
	      /home/nevia-admin/Downloads/opencv-4.1.0/build/OpenFace/build/bin/FeatureExtraction -f "$video" -aus -pose -tracked -out_dir "$dirout" &
      fi
      
   done
   
   # wait for background to finish
   wait
   
   # print that this batch is finished
   echo "$(date) finished: $i to $end" >> "$dirout/logfile.txt"
   
done

