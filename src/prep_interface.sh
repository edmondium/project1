#!/bin/sh
	 
if [ "$#" == 0 ]; then
  echo 'Missing arguments.'
  exit 1
fi

[ -d interface ] || mkdir interface

files=$( grep -l 'c begin interface' $* )

for file in $files
do
  { 
  while IFS= read -r a
  do
    if [ "$a" == "c begin interface" ]; then
      IFS= read -r a || exit 1
      name=$(expr "$a" : ' *subroutine *\([^( ]*\)')
      if [ ! "$name" ]; then
        echo 'Subroutine not found.'
	exit 1
      fi
      echo "Prepare interface "$name" in "$file
      echo "      interface" > interface/$name.inc
    fi
    if [ "$name" ]; then      
      use=$(expr "$a" : '.* use[ ,]')
      inc=$(expr "$a" : '.* !inc')
      if [ $use == 0 -o $inc != 0 ]; then
        if [ "$a" == "c end interface" ]; then
          echo "      end subroutine "$name >> interface/$name.inc
          echo "      end interface" >> interface/$name.inc
          name=''
        else
          echo "$a" >> interface/$name.inc
        fi
      fi
    fi
  done
  } < $file
done  
