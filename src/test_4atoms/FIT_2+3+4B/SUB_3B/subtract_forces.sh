


###############################################################################
###############################################################################
########                                                               ########
######## DID YOU SET A_pen TO ZERO IN THE PARAMETER FILE FIRST?!?!?!?! ########
########                                                               ########
###############################################################################
###############################################################################

FRAMES=251
TRFILE="dump_2b_sub.xyzf"	# Expected to end in ".xyzf," forces in Hartree/Bohr
BASEFL="subtract_forces.in"

cp $BASEFL subtr_forces.in

# Break the input traj apart

python /Users/lindsey11/Desktop/FORCE_MATCHING_VERSIONS/CURR_TRACKED-GEN/contrib/subtract_forces/break_apart_xyz.py $FRAMES $TRFILE

# Subtract forces on each frame, cat together

rm -f  ${TRFILE%.*}_forces_subtracted.xyzf
rm -f  ${TRFILE%.*}_forces_subtracted-FORCES.dat

LIST=`ls ${TRFILE%.*}_\#*f`

for i in $LIST
do
	IDX=${i%.*}; IDX=${IDX#*\#}
	FRCFILE=`ls ${i%\#*}FORCES_\#${IDX}*`

	# Replace frame filename
	
	awk -v CURR="${i}"       '/CRDFILE/{print; getline; sub($1,CURR)}{print}' subtr_forces.in > tmp; mv tmp subtr_forces.in
	awk -v CURR="${FRCFILE}" '/SUBTFRC/{print; getline; sub($2,CURR)}{print}' subtr_forces.in > tmp; mv tmp subtr_forces.in

	# Do the subtraction
	
	/Users/lindsey11/Desktop/FORCE_MATCHING_VERSIONS/TRACKED_4BODY/src/house_md < subtr_forces.in | tee subtr_forces.log | grep RMS
	
	cat ${i}_forces_subtracted.xyz >>  ${TRFILE%.*}_forces_subtracted.xyzf
done

# Cleanup

rm -f *_\#*
rm -f traj*
rm -f subtr_forces.*
rm -f md_statistics.out 

	
	
