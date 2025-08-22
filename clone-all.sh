#!/bin/bash

if [ ! -d "imports" ] ; then
	 mkdir imports
fi

cd imports

if [ ! -d chimes_calculator ] ; then

    # If chimes_lsq compilation isn't working, try cloning via: 
    # git clone git@github.com:RKLindsey/chimes_calculator-myLLfork.git chimes_calculator
    
    git clone https://github.com/rk-lindsey/chimes_calculator.git chimes_calculator
    cd chimes_calculator
    git checkout main
    ./install.sh 
    cd - 2>&1> /dev/null
    
    if [ $? -ne 0 ] ; then
        echo 'git clone failed'
        exit
    fi
else
    echo "chimes_calculator already found in imports. Not cloning."
fi
