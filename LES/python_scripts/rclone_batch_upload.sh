#!/bin/bash -l
#PBS -N rclone
#PBS -l walltime=0:23:59:00,nodes=1:ppn=1,pmem=1000mb

echo "This script is used for automatically uploading folders to Google Drive"
echo "RCLONE is used in this script. You should configure it before use."

## User definitions
folderPrefix="your_prefix"
localParentFolder="local_folder"
remoteParentFolder="remote:parent_folder"
iStart=1 
iEnd=10
iStep=1
# mode=1: upload. mode=2: download
mode=1 

## Main program
if [ $mode -eq 1 ]
then
  for idx in $(seq $iStart $iStep $iEnd)
    do 
      echo "Processing ${folderPrefix}_${idx}"
      rclone mkdir $remoteParentFolder/${folderPrefix}_${idx}
      rclone -v --log-file=$HOME/rclone_${folderPrefix}_${idx}.log --tpslimit 5 copy ${localParentFolder}/${folderPrefix}_${idx}/ $remoteParentFolder/${folderPrefix}_${idx}/
      ## Set a 10-minutes pause between the two uploading task
      sleep 10m
    done
fi

echo "Mission completed."
