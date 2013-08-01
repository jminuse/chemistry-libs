import os, string, sys, re, shutil
import jsub, utils

def job(atoms, basis, queue, run_name, job_type, extra_section='', procs=None, alternate_coords=None):
	head = '''#t '''+basis+''' '''+job_type+'''

run by gaussian.py

0,1
'''
	if alternate_coords:
		xyz = '\n'.join( ["%s %f %f %f" % ((a.element,)+tuple(alternate_coords[i])) for i,a in enumerate(atoms)] ) + '\n\n'
	else:
		xyz = '\n'.join( [( "%s %f %f %f" % (a.element, a.x, a.y, a.z) ) for a in atoms] ) + '\n\n'

	os.chdir('gaussian')
	with open(run_name+'.inp', 'w') as inp:
		inp.write(head+xyz+extra_section)
	#os.system('g09sub '+run_name+' -chk -disk 8192 -nproc 4 -queue huge -xhost sys_eei sys_icse')
	os.system('g09sub '+run_name+' -chk -queue '+queue+(' -nproc '+str(procs)+' ' if procs else '')+' -xhost sys_eei sys_icse')
	#os.system('g09sub '+run_name)
	os.chdir('..')

def parse_coords(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None

	energy_line = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
	energy = float(re.search('SCF Done: +\S+ += +(\S+)', energy_line).group(1))

	last_coordinates = contents.rindex('Coordinates (Angstroms)')

	start = contents.index('---\n', last_coordinates)+4
	end = contents.index('\n ---', start)
	coords = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		element = columns[1]
		x,y,z = [float(s) for s in columns[3:6]]
		coords.append( (x, y, z) )
	return energy, coords
	
def parse_chelpg(input_file):
	#os.system('ls gaussian') #test
	#print ''
	print os.getcwd()
	
	with open(input_file) as inp:
		contents = inp.read()
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

def minimize(atoms, levels_of_theory, queue='batch', name=''): #blocks until done
	for theory in levels_of_theory:
		run_name = utils.unique_filename('gaussian/', 'min_'+name+'_'+theory[:8].translate( string.maketrans('/(),*', '_____') ), '.log')
		
		job(atoms, theory, queue, run_name, 'Opt')
		jsub.wait(run_name)
		energy, coords = parse_coords('gaussian/'+run_name+'.log')
		if coords:
			for i,xyz in enumerate(coords):
				atoms[i].x, atoms[i].y, atoms[i].z = xyz
			return run_name, energy

def chelpg(atoms, theory, queue='batch', chkfile_run_name=None, name=''):
	run_name = utils.unique_filename('gaussian/', 'chelpg_'+name+'_'+theory[:8].translate( string.maketrans('/(),*', '_____') ), '.log')
	if chkfile_run_name:
		shutil.copyfile('gaussian/'+chkfile_run_name+'.chk', 'gaussian/'+run_name+'.chk')
	job(atoms, theory, queue, run_name, 'Pop=CHelpG')
	jsub.wait(run_name)
	charges = parse_chelpg('gaussian/'+run_name+'.log')
	if charges:
		for i,charge in enumerate(charges):
			atoms[i].charge = charge
		return charges

def energy(atoms, bond, theory, queue, chkfile, async=False, name='', alternate_coords=None):
	run_name = utils.unique_filename('gaussian/', 'energy_'+name, '.inp')
	#if chkfile:
	#	shutil.copyfile('gaussian/'+chkfile+'.chk', 'gaussian/'+run_name+'.chk')
	job(atoms, theory, queue, run_name, 'SP', procs=1, alternate_coords=alternate_coords)
	if async:
		return run_name
	else:
		jsub.wait(run_name)
		energy, coords = parse_coords('gaussian/'+run_name)
		return energy, coords

def bond_energy(atoms, bond, theory, queue, chkfile, async=False, inc=None, name=''):
	run_name = utils.unique_filename('gaussian/', 'bond_'+name, '.inp')
	#if chkfile:
	#	shutil.copyfile('gaussian/'+chkfile+'.chk', 'gaussian/'+run_name+'.chk')
	job(atoms, theory, queue, run_name, 'SP Geom=ModRedundant', extra_section='B %d %d +=%f'%(bond.atoms[0].index, bond.atoms[1].index, inc or bond.d*0.05), procs=1 )
	if async:
		return run_name
	else:
		jsub.wait(run_name)
		energy, coords = parse_coords('gaussian/'+run_name)
		return energy, coords

def angle_energy(atoms, angle, theory, queue, chkfile, async=False, inc=3.0, name=''):
	run_name = utils.unique_filename('gaussian/', 'angle_'+name, '.inp')
	#if chkfile:
	#	shutil.copyfile('gaussian/'+chkfile+'.chk', 'gaussian/'+run_name+'.chk')
	job(atoms, theory, queue, run_name, 'SP Geom=ModRedundant', extra_section='A %d %d %d +=%f'%(angle.atoms[0].index, angle.atoms[1].index, angle.atoms[2].index, inc), procs=1 )
	if async:
		return run_name
	else:
		jsub.wait(run_name)
		energy, coords = parse_coords('gaussian/'+run_name)
		return energy, coords

def dihedral_energy(atoms, dihedral, theory, queue, chkfile, async=False, inc=None, name=''):
	run_name = utils.unique_filename('gaussian/', 'dihedral_'+name, '.inp')
	#if chkfile:
	#	shutil.copyfile('gaussian/'+chkfile+'.chk', 'gaussian/'+run_name+'.chk')
	job(atoms, theory, queue, run_name, 'SP Geom=ModRedundant', extra_section='D %d %d %d %d +=%f'%(dihedral.atoms[0].index, dihedral.atoms[1].index, dihedral.atoms[2].index, dihedral.atoms[2].index, inc), procs=1 )
	if async:
		return run_name
	else:
		jsub.wait(run_name)
		energy, coords = parse_coords('gaussian/'+run_name)
		return energy, coords

