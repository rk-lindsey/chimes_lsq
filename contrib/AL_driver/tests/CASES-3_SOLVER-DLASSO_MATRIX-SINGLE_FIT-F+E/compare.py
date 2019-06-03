import sys

# Read in the two parameter files, and only save non-comment and non-empty lines

if len(sys.argv) != 3:

	print "ERROR: Usage is: python <this script> paramter_file_1.dat parameter_file_2.dat"
	print "Received:",sys.argv
	print "Exiting"
	
	exit()

file_1 = open(sys.argv[1],'r')
file_2 = open(sys.argv[2],'r')

tmp_1 = file_1.readlines()
tmp_2 = file_2.readlines()

contents_1 = []
contents_2 = []

for i in xrange(len(tmp_1)):

	tmp_1[i] = tmp_1[i].split()
	
	if len(tmp_1[i]) > 0 and tmp_1[i][0] != '!':	
		contents_1.append(tmp_1[i])
		
for i in xrange(len(tmp_2)):

	tmp_2[i] = tmp_2[i].split()
	
	if len(tmp_2[i]) > 0 and tmp_2[i][0] != '!':	
		contents_2.append(tmp_2[i])	
		

# Compare the two files. If the last line entry is a float with >=8 digits(contains a decimal), only compare the first 6 digits

if len(contents_1) != len(contents_2):

	print "ERROR: File contents are of different lengths"
	exit()
	
for i in xrange(len(contents_1)):

	if len(contents_1[i]) != len(contents_2[i]):

		print "ERROR: File contents lines are of different lengths"
		print contents_1[i]
		print contents_2[i]
		
		exit()	

	# Check the first n-1 columns of two lines

	mismatch = False

	for j in xrange(len(contents_1[i])-1):
	
		if contents_1[i][j] != contents_2[i][j]:
		
			mismatch = True	
		
	# Check the last column of the two lines, treating special if float
	
	lastcol_1 = contents_1[i][len(contents_1[i])-1]
	lastcol_2 = contents_2[i][len(contents_2[i])-1]
	
	keep = 6

	
	if '.' in lastcol_1 or '.' in lastcol_2:
	
		lastcol_1.strip('-')
		lastcol_2.strip('-')
		
		lastcol_1.strip('+')
		lastcol_2.strip('+')
		
		lastcol_1.split('.')
		lastcol_2.split('.')
		
		pre_1 = len(lastcol_1[0])
		pre_2 = len(lastcol_2[0])
		
		post_1 = len(lastcol_1[1])
		post_2 = len(lastcol_2[2])


		if pre_1 > keep:
			tmp_1 = lastcol_1[0][0:keep]
		else:
			tmp_1 = lastcol_1[0] + '.' + lastcol_1[1][0:(keep-pre_1)]
			
		if pre_2 > keep:
		
			tmp_2 = lastcol_2[0][0:keep]	
		else:
			tmp_2 = lastcol_2[0] + '.' + lastcol_2[1][0:(keep-pre_2)]						

		
		if tmp_1 != tmp_2:

			print lastcol_1
			print lastcol_2
			print tmp_1
			print tmp_2
			
			mismatch = True
	else:
	
		if lastcol_1 != lastcol_2:
		
			mismatch = True
			
	if mismatch:
	
		print "Found mismatch"
		print '\t',contents_1[i]
		print '\t',contents_1[i]
		
	
		
		
	
