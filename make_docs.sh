#!/bin/bash


# Check Sphinx version

ver=`sphinx-build --version 2>&1 | awk '/Sphinx/{split($2,s,"."); split(s[1],v,"v");print v[2]}'`

if [ $ver -lt 3 ] ; then 
	echo "Error: Please upgrade Sphinx to version 3 or greater"
	exit 0
fi

# Remove any existing documentation .pdf

rm -f chimes_lsq.pdf


# Remove any build files, compile 

cd doc

make clean
make html
make latex

# Convert built latex files to a .pdf, move to the docs folder

cd build/latex

TEXFILE="chimes_lsq.tex"
pdflatex  $TEXFILE


cd ../..

cp build/latex/chimes_lsq.pdf .


echo ""
echo "*****************************************************************"
echo "                       DOCUMENTATION BUILT!                      "
echo "*****************************************************************"
echo ""
echo " To view documentation, either:"
echo " (1) Open ./docs/build/html/index.html in a browser (recommended)"
echo " (2) View ./docs/chimes_lsq.pdf"
echo ""
echo "Type \"make clean\" in ./doc to remove all built files." 
echo ""
echo "*****************************************************************"
echo "" 
