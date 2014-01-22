import os, string, sys, re, shutil
import jsub, utils

def job(atoms, basis, queue, run_name, job_type, extra_section='', procs=None, alternate_coords=None, charge_and_multiplicity='0,1'):
	counterpoise = 'counterpoise' in job_type.lower()
	head = '#t '+basis+' '+job_type+'\n\nrun by gaussian.py\n\n'+charge_and_multiplicity+'\n'
	if alternate_coords:
		xyz = '\n'.join( ["%s %f %f %f" % ((a.element,)+tuple(alternate_coords[i])) for i,a in enumerate(atoms)] ) + '\n\n'
	else:
		if atoms and type(atoms[0])==type([]): #multiple lists of atoms (e.g. transistion state calculation)
			xyz = 'run by gaussian.py\n\n0,1\n'.join([('\n'.join( [( "%s %f %f %f" % (a.element, a.x, a.y, a.z) ) for a in atom_list] ) + '\n\n') for atom_list in atoms])
		else: #single list of atoms
			if not counterpoise:
				xyz = '\n'.join( [( "%s %f %f %f" % (a.element, a.x, a.y, a.z) ) for a in atoms] ) + '\n\n'
			else:
				xyz = '\n'.join( [( "%s(Fragment=%d) %f %f %f" % (a.element, a.fragment, a.x, a.y, a.z) ) for a in atoms] ) + '\n\n'

	os.chdir('gaussian')
	with open(run_name+'.inp', 'w') as inp:
		inp.write(head+xyz+extra_section)
	os.system('g09sub '+run_name+' -chk -queue '+queue+((' -nproc '+str(procs)+' ') if procs else '')+' -xhost sys_eei sys_icse')
	os.chdir('..')

def restart_job(old_run_name, job_type='ChkBasis Opt=Restart', queue='batch', procs=None):
	run_name = old_run_name+'r'
	os.chdir('gaussian')
	shutil.copyfile(old_run_name+'.chk', run_name+'.chk')
	with open(run_name+'.inp', 'w') as inp:
		inp.write('#t '+job_type+'\n\nrun by gaussian.py\n\n')
	os.system('g09sub '+run_name+' -chk -queue '+queue+((' -nproc '+str(procs)+' ') if procs else '')+' -xhost sys_eei sys_icse')
	os.chdir('..')

def parse_coords(input_file, get_modredundant=False, get_energy=True):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None

	if get_energy:
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
	
	if get_modredundant:
		try:
			dihedral = [int(re.search('The following ModRedundant input section has been read:\s+D'+('\s+(\S+)'*4), contents).group(i)) for i in range(1,5)]
		except: return None
		if get_energy:
			return energy, coords, dihedral
		else:
			return coords, dihedral
	else:
		if get_energy:
			return energy, coords
		else:
			return coords

def scan_modredundant(starting_atoms, basis, run_name_base, modredundant_section, param_starts, param_ends, n_steps):
	ranges = [param_ends[i]-param_starts[i] for i in range(len(param_ends))]
	increments = [r/n_steps for r in ranges]
	params = [p for p in param_starts]
	for step in xrange(n_steps):
		run_name = run_name_base+str(step)
		if step==0:
			job(starting_atoms, basis, 'batch', run_name, 'Opt=(CalcFC,ModRedundant)', procs=1, extra_section=modredundant_section % tuple(params) )
		else:
			old_run_name = run_name_base+str(step-1)
			atoms = parse_atoms('gaussian/'+old_run_name+'.log', check_convergence=False)[1]
			shutil.copy('gaussian/'+old_run_name+'.chk', 'gaussian/'+run_name+'.chk')
			job(atoms, basis, 'batch', run_name, 'Opt=(CalcFC,ModRedundant) Guess=Read', procs=1, extra_section=modredundant_section % tuple(params) )
		jsub.wait(run_name)
		for i in range(len(params)):
			params[i] += increments[i]

def parse_atoms(input_file, get_energy=True, check_convergence=True, get_time=False, counterpoise=False):
	contents = open(input_file).read()
	if check_convergence and get_energy and 'Normal termination of Gaussian 09' not in contents:
		return None

	if 'Summary of Optimized Potential Surface Scan' in contents:
		end_section = contents[contents.rindex('Summary of Optimized Potential Surface Scan'):]
		energy_lines = re.findall('Eigenvalues -- ([^\\n]+)', end_section)
		energy = [float(s) for line in energy_lines for s in re.findall('-[\d]+\.[\d]+', line)]
		
		minima = re.split('Stationary point found', contents)
		atoms = []
		for m in minima[1:]:
			coordinates = m.index('Coordinates (Angstroms)')
			
			start = m.index('---\n', coordinates)+4
			end = m.index('\n ---', start)
			atoms.append([])
			for line in m[start:end].splitlines():
				columns = line.split()
				element = columns[1]
				x,y,z = [float(s) for s in columns[3:6]]
				atoms[-1].append( utils.Struct(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z) )
		
		if get_energy:
			return energy, atoms
		
	elif get_energy:
		if not counterpoise:
			energy_line = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
			energy = float(re.search('SCF Done: +\S+ += +(\S+)', energy_line).group(1))
			#print '\n'.join(re.findall('Counterpoise: corrected energy = +(\S+)', contents))
			#print '\n'.join(re.findall('SCF Done: +\S+ += +(\S+)', contents))
		else:
			energy = float(re.findall('Counterpoise: corrected energy = +(\S+)', contents)[-1])

	last_coordinates = contents.rindex('Coordinates (Angstroms)')

	start = contents.index('---\n', last_coordinates)+4
	end = contents.index('\n ---', start)
	atoms = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		element = columns[1]
		x,y,z = [float(s) for s in columns[3:6]]
		atoms.append( utils.Struct(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z) )
	
	if get_time:
		m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
		time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))
		return energy, atoms, time
	if get_energy:
		return energy, atoms
	else:
		return atoms
	
def parse_starting_coords(input_file):
	contents = open(input_file).read()
	first_coordinates = contents.index('Charge =  0 Multiplicity = 1')
	start = contents.index('\n', first_coordinates)+1
	end = contents.index('\n \n', start)
	coords = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		element = columns[0]
		x,y,z = [float(s) for s in columns[1:4]]
		coords.append( (x, y, z) )
	return coords
	
def parse_scan(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None
	scan_steps = contents.split('on scan point')
	energy_list = []
	atoms_list = []
	
	scan_steps = [ scan_steps[i] for i in range(1,len(scan_steps)-1) if scan_steps[i][:10].split()[0]!=scan_steps[i+1][:10].split()[0] ]
	#print [int(s[:10].split()[0]) for s in scan_steps]
	#print len(scan_steps)
	
	for scan_step in scan_steps:
		energy_line = scan_step[scan_step.rindex('SCF Done'):scan_step.index('\n', scan_step.rindex('SCF Done'))]
		energy = float(re.search('SCF Done: +\S+ += +(\S+)', energy_line).group(1))

		last_coordinates = scan_step.rindex('Coordinates (Angstroms)')

		start = scan_step.index('---\n', last_coordinates)+4
		end = scan_step.index('\n ---', start)
		atoms = []
		for line in scan_step[start:end].splitlines():
			columns = line.split()
			element = columns[1]
			x,y,z = [float(s) for s in columns[3:6]]
			atoms.append( utils.Struct(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z) )
		energy_list.append(energy)
		atoms_list.append(atoms)
	return energy_list, atoms_list
	
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

def minimize(atoms, theory, queue='batch', name='', restrained=None, async=False): #blocks until done
	run_name = utils.unique_filename('gaussian/', 'min_'+name+'_'+theory[:8].translate( string.maketrans('/(),*', '_____') ), '.inp')
	
	if not restrained:
		job(atoms, theory, queue, run_name, 'Opt=CalcFC') #SCRF(Solvent=n-Hexane)
	else:
		dihedral = restrained
		job(atoms, theory, queue, run_name, 'Opt=(ModRedundant,Loose,CalcFC)', extra_section='D %d %d %d %d F'%tuple([a.index for a in dihedral.atoms]))
	if async:
		return run_name
	else:
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
	try:
		charges = parse_chelpg('gaussian/'+run_name+'.log')
	except: #sometimes file is not written yet
		time.sleep(10)
		charges = parse_chelpg('gaussian/'+run_name+'.log')
	if charges:
		for i,charge in enumerate(charges):
			atoms[i].charge = charge
		return charges

def energy(atoms, theory, queue='batch', chkfile=None, async=False, name='', alternate_coords=None):
	run_name = utils.unique_filename('gaussian/', 'E_'+name+'_'+theory[:8].translate( string.maketrans('/(),*', '_____') ), '.inp')
	if chkfile:
		shutil.copyfile('gaussian/'+chkfile+'.chk', 'gaussian/'+run_name+'.chk')
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

