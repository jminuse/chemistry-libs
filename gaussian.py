import os, string, sys, re

def job(atoms, basis, queue, run_name, job_type, extra_section='', chkfile=None):
	processors = {'batch':'4','huge':'8','long':'8'}[queue]
	memory = {'batch':'1GB','huge':'8GB','long':'4GB'}[queue]
	if not chkfile:
		chkfile = +run_name+'.chk'
	head = '''#!/bin/csh
##NBS-queue: '''+queue+'''
##NBS-nproc: '''+processors+'''
##NBS-speed: 3000
##NBS-name: "'''+run_name+'''"

setenv g09root /usr/local/g09
source $g09root/g09/bsd/g09.login
g09 <<END > '''+run_name+'''.log
%NProcShared='''+processors+'''
%Mem='''+memory+'''
%Chk='''+chkfile+'''
#t '''+basis+''' '''+job_type+'''

gaussian.py job

0,1
'''

	tail = '''

eof
cd /tmp
rm *.rwf
'''
	xyz = '\n'.join( ['\t'.join([str(s) for s in [a.element, a.x, a.y, a.z]]) for a in atoms] ) 

	os.chdir('gaussian')
	open(run_name+'.gjf', 'w').write(head+xyz+extra_section+tail)
	os.system('jsub '+run_name+'.gjf')
	os.chdir('..')

def parse_coords(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None

	a = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
	print a

	last_coordinates = contents.rindex('Coordinates (Angstroms)')

	start = contents.index('---\n', last_coordinates)+4
	end = contents.index('\n ---', start)
	coords = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		element = columns[1]
		x,y,z = [float(s) for s in columns[3:6]]
		coords.append( (x, y, z) )
	return coords
	
def parse_chelpg(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None
	
	start = contents.rindex('Fitting point charges to electrostatic potential')
	end = contents.index('-----------------', start)
	charges = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		if len(columns)==3:
			charges.append( float(columns[2]) )
	return charges

def minimize(atoms, levels_of_theory, queue='batch'): #blocks until done
	for theory in levels_of_theory:
		run_number = 0
		while True:
			run_name = 'minimize_'+theory.translate( string.maketrans('/(),', '----') )+'_'+run_number
			run_number += 1
			if not os.path.exists(run_name+'.log'): break
		
		job(atoms, theory, queue, run_name, 'Opt')
		while(True):
			jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
			if run_name in jlist:
				os.sleep(60)
			else:
				break
		coords = parse_coords(run_name+'.log')
		if coords:
			for i,xyz in enumerate(coords):
				atoms[i].x, atoms[i].y, atoms[i].z = xyz
			return run_name

def chelpg(atoms, level_of_theory, queue='batch', chkfile=None):
	run_number = 0
	while True:
		run_name = 'chelpg_'+theory.translate( string.maketrans('/(),', '----') )+'_'+run_number
		run_number += 1
		if not os.path.exists(run_name+'.log'): break
	
	job(atoms, theory, queue, run_name, 'Pop=CHelpG', chkfile=chkfile)
	while(True):
		jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
		if run_name in jlist:
			os.sleep(60)
		else:
			break
	charges = parse_chelpg(run_name+'.log')
	if charges:
		for i,charge in enumerate(charges):
			atoms[i].charge = charge
		return run_name

def bond_energies(xyz, bonds, theory):
	pass

def angle_energies(xyz, angles, theory):
	pass

def dihedral_energies(xyz, dihedrals, theory):
	pass

