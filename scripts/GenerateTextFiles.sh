#!/bin/bash

#
# Generate the input files for the carving algorithm
#  ./GenerateTextFiles.sh  ImagesFolder CamerasFolder
#

if [[ ! "$1" =~ ^/ ]];then
  ls -1 $PWD/$1/*.png > Images.txt
  ls -1 $PWD/$2/*.txt > Cameras.txt
else
  ls -1 $1/*.png > Images.txt
  ls -1 $2/*.txt > Cameras.txt
fi

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#$DIR/../../build/bin/Carving Cameras.txt Images.txt $3 $4
