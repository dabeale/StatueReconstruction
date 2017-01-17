#!/bin/bash

echo "Checking pre-requisites"
if [ ! -e  /usr/local/lib/libopencv_core.dylib ]
then 
  echo "OpenCV is not installed"
  echo "Install from homebrew using the command:"
  echo "  brew install opencv"
  exit 1
fi

echo "Building C files"
qm=$(mdfind -name qmake | grep bin/qmake | grep 5 | tail -n 1)
qml=$(echo $qm | wc -c)
nargs=$#
if((qml == 0 && nargs == 0))
then
  echo "qmake must be version 5 or higher"
  echo "If you are sure you have qmake and know where it is installed"
  echo "try again using ./Build.sh pathtoqmake" 
  exit 1
else
  if(( nargs == 0 ))
  then
    mkdir build
    cd build
    #$qm -spec macx-clang-omp ../src/StatueReconstruction.pro CONFIG+=release CONFIG+=x86_64 #-d
    $qm ../src/StatueReconstruction.pro CONFIG+=release CONFIG+=x86_64 #-d 
  else
    qm=$1
    if [ -e $qm ]
    then
      mkdir build
      cd build
      $qm ../src/StatueReconstruction.pro CONFIG+=release CONFIG+=x86_64 #-d 
    else
      echo $qm "is not a valid location for qmake"
      exit 1;
    fi
  fi 
  make -f Makefile -j4
  cd ..
fi

echo "Build complete."
