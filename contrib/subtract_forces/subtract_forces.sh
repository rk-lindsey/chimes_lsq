#!/bin/bash

# Run with: ./this_file.sh <frames> <input.xyzf file> <input md file (i.e. run_md-compare_force.in)>

# Read input from arguments

IN_XYZF=$2	# ="input-orig.xyzf"	 	# Input xyzf file.. MUST end in .xyzf
IN_FRMS=$1	# =47 			  	# Frames in the above file
IN_MDFL=$3	# ="run_md-compare_force.in"	# Input file for a compare force md run

NEW_TRAJ="${IN_XYZF%.*}_#ALL-subtracted.xyzf"

echo ""
echo "Will process: "
echo "" 
echo "	frames:        $IN_FRMS"
echo "	xyzf file:     $IN_XYZF"
echo "	md input file: $IN_MDFL"
echo ""
echo "Will produce:  "
echo "	output file:   $NEW_TRAJ"
echo ""

# Set the source code base directory

CONTRB_BASE=`dirname $0`
SOURCE_BASE="${CONTRB_BASE}/../../src/"

BREAK_APART="${CONTRB_BASE}/break_apart_xyz.py"
SUBTR_FORCE="${CONTRB_BASE}/subtract_forces.py"
HOUSE_MD_EX="${SOURCE_BASE}/house_md"

if   [ ! -e $BREAK_APART ] ; then
	echo "ERROR: Cannot find script:"
	echo $BREAK_APART
	exit 0
elif [ ! -e $SUBTR_FORCE ] ; then
	echo "ERROR: Cannot find script:"
	echo $SUBTR_FORCE
	exit 0
elif [ ! -e $HOUSE_MD_EX ] ; then
	echo "ERROR: Cannot find executable:"
	echo $HOUSE_MD_EX
	exit 0	
fi	

############ START THE SCRIPT

echo "Searching for break apart files named like: ${IN_XYZF%.*}_#<number>.xyzf"
echo "Beware, these files will be deleted after creation/processing"
echo "as will files of the type *_#ALL-subtracted*"
echo ""

rm -f $NEW_TRAJ
rm -f *_#ALL-subtracted*

# 1. Break the original .xyzf file into its individual components

python $BREAK_APART $IN_FRMS $IN_XYZF

# 2. For each file, run the compare md code, and update xyzf with the resulting forces 

for i in `ls ${IN_XYZF%.*}_*.xyzf`
do
	# 0. Get the number of frames associated with the file
	
	FRAMES=${i%.*}
	FRAMES=${FRAMES#*#}
	
	# 1. Replace xyzf file name in the md input file 
	
	awk -v newstr="${i}" '/CRDFILE/{print;getline;sub($1,newstr)}{print}' run_md-compare_force.in > tmp.in

	# 2. Run the compare force calculation
	
	$HOUSE_MD_EX < tmp.in > /dev/null # We don't need the default output of this md call.
	
	# 3. Update the .xyzf file by subtracting off the computed forces 
	
	python $SUBTR_FORCE 1 $i forceout.txt
	
	NEW="${i%.*}_subtracted.xyzf"
	
	# 4. Re-assemble the .xyzf trajectory file
	
	NEW="${i%.*}_subtracted.xyzf"
	
	cat $NEW >> $NEW_TRAJ
	
	rm -f $NEW $i
done

rm -f forceout.txt forceout-labeled.txt



