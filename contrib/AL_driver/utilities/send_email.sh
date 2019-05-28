#!/bin/bash

address=$1
shift
status=""

for i 
do
	status="${status} $i"
done

echo $status | mailx -s "ALC Driver Update" $address
