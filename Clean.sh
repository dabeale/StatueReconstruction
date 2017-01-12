#!/bin/bash

if [[ $1 == "oct" ]]
then
  rm -rfv oct/*
  exit 0
fi

find . -name .qmake.stash -delete
find . -regex .*~ -print -delete
find . -regex .*.pro.user -print -delete
find . -name Makefile -print -delete
rm -rfv build*
