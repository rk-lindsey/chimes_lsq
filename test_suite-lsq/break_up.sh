#!/bin/bash

for i in `find . -name A.txt`;
do 
	dir=`dirname $i`
	
	cd $dir
	
	echo "Breaking up files in $dir"
	
	# Breaks up into 95M chunks, i.e., < Github's 100M file size limit
	# Then removes the big single A.txt file
	
	split -b95M A.txt A.txt.
	rm -f A.txt

	cd - > /dev/null 2>& 1
done
