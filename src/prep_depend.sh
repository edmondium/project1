#!/bin/sh

echo 'Build dependencies.' 
{
  MOD=$(grep -l "end module" $*)
  module=$( grep -h 'end module' $MOD | awk '{print $3}')
  module=$(sed 's/ /\n/g' <<< $module | sort -u)
  i=0
  for mod in $module 
  do
    modfile=$( grep -l 'end module '$mod' *$' $MOD )
    depfile=$( grep -l 'use '$mod $*)
    if [ $(echo "$modfile" | wc -w) -gt 1 ]; then
      if [ $mod != "readwrite" ]; then
        echo "ERROR: Multiple module definition is not for module readwrite" > /dev/stderr
        exit
      fi  
      for file in $modfile      
      do
        if [ $file != "readwrite.f" -a ${file:0:10} != "readwrite_" ]; then	  
          echo "ERROR: File for module readwrite named wrongly: $file" > /dev/stderr
          exit
	fi
      done
      modfile='readwrite$(DFT).f'
    fi
    echo MOD$i = $modfile
    echo -n DEP$i =
    for file in $depfile
    do
      [ $file != $modfile ] && echo -n ' '$file
    done
    echo
    let i=i+1
  done
  n=$i
  i=0
  while [ $i -lt $n ]
  do
    echo '$('DEP$i:%.f=%.o'): $(findstring $('MOD$i':%.f=%.o),$(SRC:%.f=%.o)) $@'
    let i=i+1
  done 
} > depend.mk

