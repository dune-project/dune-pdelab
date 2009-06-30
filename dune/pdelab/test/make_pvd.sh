#!/bin/bash

prefix="$1"

# .pvd (output) file name
pvdFile="$prefix.pvd"

# File prefix, ending and number of files in your animation series
# Filename generation:  filename=$filePrefix$fileNumber$ending
# with filenumber=`printf $formating $i`
# e.g. filename=hairy-0045.vtu

suffix=.vtu
formatting="%d"

# pvd xml file header
echo  "<?xml version=\"1.0\"?>" > "$pvdFile"
echo  "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" >> "$pvdFile"
echo "<Collection>" >> "$pvdFile"

# Files list
i=0
while filename=`printf "%s.$formatting%s" "$prefix" "$i" "$suffix"`
      [ -e "$filename" ]
do echo "<DataSet timestep=\"$i \" part=\"001\" file=\"$filename\"/>" >> "$pvdFile"
   i=`expr $i + 1`
done

# pvd end
echo "</Collection>" >> "$pvdFile"
echo "</VTKFile>" >> "$pvdFile"
