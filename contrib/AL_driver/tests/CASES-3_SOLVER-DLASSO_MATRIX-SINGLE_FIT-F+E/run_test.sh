#!/bin/bash

# Usage:
#
# ./run_test.sh run 		# Runs the test job
# ./run_test.sh compare 	# Compares resulting parameter files against expected results
# ./run_test.sh clean   	# Removes newly-generated files/directories
#
# Notes: I recommend running this script from screen


if [ "$#" -ne 1 ]; then

	echo "Error: Exactly one argument is expected."
	echo "Received ${#}."
	echo "Allowed arguments are run, compare, or clean"
	echo "Exiting."
	exit 0
fi

if [[ "$1" == "run" ]] ; then

	# Set up the test directory...

	if [ -d "run_test" ]; then  rm -rf run_test; fi

	mkdir run_test
	cp -r ALL_BASE_FILES run_test

	mkdir -p run_test/ALC_src
	cp ../../*.* run_test/ALC_src
	cp -r ../../utilities run_test/ALC_src

	# Run the test

	cp config.py run_test
	cd run_test
	unbuffer python ../../../../PYTHON_DRIVER/main.py 0 1 2 | tee driver.log 
	cd ..

elif [[ "$1" == "compare" ]] ; then

	tot_diff_lines=0

	for i in 0 1 2
	do
		echo -e "\nComparing ALC-${i} parameter files"
		
		python compare.py expected_output/ALC-${i}.params.txt run_test/ALC-${i}/GEN_FF/params.txt > ALC-${i}.diff

		lines=`wc -l ALC-${i}.diff | awk '{print $1}'`
		
		if [ $lines -gt 0 ] ; then
		
			echo -e "\tfiles differ by $lines lines."
		fi
		
		let tot_diff_lines=tot_diff_lines+$lines
	done
	
	if [ $tot_diff_lines -lt 1 ] ; then
	
		echo -e "\nNo differences found!"
	fi


elif [[ "$1" == "clean" ]] ; then

	rm -rf run_test
	rm -f  ALC-{0..2}.diff

else
	echo "ERROR: Unknown option $1."
	echo "Allowed arguments are run, compare, or clean"
	echo "Exiting."
	exit 0
fi
