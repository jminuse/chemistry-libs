import os, string, sys, re

def job(run_name, basis, in_file, out_file, queue):
	processors = {'batch':'4','huge':'8','long':'8'}[queue]
	memory = {'batch':'1GB','huge':'8GB','long':'4GB'}[queue]
	head = '''#!/bin/csh
##NBS-queue: '''+queue+'''
##NBS-nproc: '''+processors+'''
##NBS-speed: 3000
##NBS-name: "'''+run_name+'''"

setenv g09root /usr/local/g09
source $g09root/g09/bsd/g09.login
g09 <<END > '''+out_file+'''
%NProcShared='''+processors+'''
%Mem='''+memory+'''
#t '''+basis+''' Opt !Pop=CHelpG !SP !Counterpoise=2 Opt Freq !SCRF=(Solvent=n-Pentane)

opt before counterpoise

0,1
'''

	tail = '''

eof
cd /tmp
rm *.rwf
'''
	atoms = [ [line.split()[0]]+[float(s) for s in line.split()[1:]] for line in open(in_file)]

	xyz = '\n'.join( ['\t'.join([str(s) for s in a]) for a in atoms] ) 

	os.chdir('gaussian')
	open(run_name+'.gjf', 'w').write(head+xyz+tail)
	os.system('jsub '+run_name+'.gjf')
	os.chdir('..')

def parse(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		print 'Job did not finish'

	a = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
	print a

	last_coordinates = contents.rindex('Coordinates (Angstroms)')

	start = contents.index('---\n', last_coordinates)+4
	end = contents.index('\n ---', start)
	lines = contents[start:end].splitlines()

	for line in lines:
		columns = line.split()
		element = columns[1]
		x,y,z = columns[3:6]
		f.write( '\t'.join([element, x, y, z]) + '\n' )

def run(levels_of_theory, queue='batch'): #blocks until done
	for theory in levels_of_theory:
		job(run_name, theory, in_file, out_file, queue)
		while true:
			check_if_done
			os.sleep(60)

def bond_energies(xyz, bonds, theory):
	

def angle_energies(xyz, angles, theory):
	

def dihedral_energies(xyz, dihedrals, theory):
	

