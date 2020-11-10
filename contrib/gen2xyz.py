import sys


traj = '.'.join(sys.argv[1].split('.')[0:-1])

print "Converting", traj+".gen", "to", traj+".xyz"

new_frame = 1

natoms  = None
atmtyps = None
latvec  = None
coords  = None

ofstream = open(traj+".xyz",'w')

with open (traj+".gen",'r') as ifstream:

	for line in ifstream:

		if new_frame >= 1:

			if new_frame == 1:
				natoms     = int(line.split()[0])
				new_frame += 1
				coords  = []
				latvec  = ""
				continue
				
			if new_frame == 2:
				atmtyps = line.split()
				new_frame += 1
				continue
				
			if (new_frame > 2) and (new_frame <= natoms + 2):
				line = line.split()
				coords.append(line[1:])
				coords[-1][0] = atmtyps[int(coords[-1][0])-1]
				coords[-1] = ' '.join(coords[-1])
				
				new_frame += 1
				continue
				
			if new_frame == natoms + 3:
				new_frame += 1
				continue
			
			if (new_frame > natoms + 3) and (new_frame <= natoms + 6):
			
				latvec += ' '.join(line.split()) + " "
				new_frame += 1
								
			if new_frame == natoms + 7:
				ofstream.write(str(natoms)+'\n')
				ofstream.write(str(latvec)+'\n')
				
				for i in xrange(natoms):
					ofstream.write(str(coords[i])+'\n')
				
				new_frame = 1
				
			
		
	
	
