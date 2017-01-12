#!/bin/bash

echo "Building C files"
qm=$(mdfind -name qmake | grep bin/qmake | grep 5 | tail -n 1)
qml=$(echo $qm | wc -c)
if((qml == 0))
then
  echo "qmake must be version 5 or higher"
  exit 1
else
  mkdir build
  cd build
  #$qm -spec macx-clang-omp ../src/StatueReconstruction.pro CONFIG+=release CONFIG+=x86_64 #-d
  $qm ../src/StatueReconstruction.pro CONFIG+=release CONFIG+=x86_64 #-d  
  make -f Makefile -j4
  cd ..
fi

echo "C Build Complete"

echo "Build complete."
