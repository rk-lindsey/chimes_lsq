#!/bin/bash

# NOTE: If this command fails, edit utilities/local_pydoc.py to call the same
# python path as you use to execute this driver

# NOTE: local_pydoc.py is simply the code located at: 
# https://github.com/python/cpython/blob/2.7/Lib/pydoc.py

rm -rf doc

cp utilities/local_pydoc.py .

for i in `ls *py`; do ./local_pydoc.py -w ${i%*.py}; done

mkdir doc
mv *html doc

rm -f local_pydoc.py

browser=`which firefox`

if [[ $browser == *"no firefox in"* ]] ; then

	echo "Cannot find browser firefox. Exiting"

	exit 0
else
	echo "Opening documentation"

	firefox doc/main.html &
fi
