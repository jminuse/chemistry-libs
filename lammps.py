import re, random, numpy, math, os, subprocess
import utils, jsub

def parse_opls_parameters(parameter_file):
	elements = {}; atom_types = []; bond_types = []; angle_types = []; dihedral_types = []
	for line in open(parameter_file):
		columns = line.split()
		if not columns: continue
		if columns[0]=='atom':
			m = re.match('atom +(\d+) +(\d+) +(\S+) +"([^"]+)" +(\d+) +(\S+) +(\d+)', line)
			atom_type = utils.Struct(index=int(m.group(1)), index2=int(m.group(2)), element_name=m.group(3), notes=m.group(4), element=int(m.group(5)), mass=float(m.group(6)), bond_count=int(m.group(7) ) )
			if atom_type.element not in elements: elements[atom_type.element] = []
			elements[atom_type.element].append(atom_type)
			if '(UA)' in atom_type.notes:
				atom_type.element = 0 #reject united-atom parameters
			atom_types.append(atom_type)
		elif columns[0]=='vdw':
			atom_types[int(columns[1])-1].vdw_r = max( float(columns[2]), 1.0)
			atom_types[int(columns[1])-1].vdw_e = max( float(columns[3]), 0.01)
		elif columns[0]=='charge':
			atom_types[int(columns[1])-1].charge = float(columns[2])
		elif columns[0]=='bond':
			bond_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:3]]),e=float(columns[3]),r=float(columns[4])) )
		elif columns[0]=='angle':
			angle_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:4]]),e=float(columns[4]),angle=float(columns[5])) )
		elif columns[0]=='torsion':
			dihedral_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:5]]),e=tuple([float(s) for s in columns[5::3]])) )
			if len(dihedral_types[-1].e)==3:
				dihedral_types[-1].e = dihedral_types[-1].e + (0.,)
	return elements, atom_types, bond_types, angle_types, dihedral_types

def guess_opls_parameters(atoms, bonds, angles, dihedrals, opls_parameter_file):
	elements, atom_types, bond_types, angle_types, dihedral_types = parse_opls_parameters(opls_parameter_file)
	
	atom_types_by_element_and_bond_count = {}
	for a in atom_types:
		if (a.element, a.bond_count) not in atom_types_by_element_and_bond_count:
			atom_types_by_element_and_bond_count[(a.element, a.bond_count)] = []
		atom_types_by_element_and_bond_count[(a.element, a.bond_count)].append(a)

	atomic_number = {'H':1, 'C':6, 'N':7}
	for a in atoms:
		a.type = atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))][0]
	
	dihedral_types_by_index2 = dict([ (t.index2s,t) for t in dihedral_types] + [ (t.index2s[::-1],t) for t in dihedral_types])
	angle_types_by_index2 = dict([ (t.index2s,t) for t in angle_types] + [ (t.index2s[::-1],t) for t in angle_types])
	bond_types_by_index2 = dict([ (t.index2s,t) for t in bond_types] + [ (t.index2s[::-1],t) for t in bond_types])
	
	def get_bond_type(options):
		for o1 in options[0]:
			for o2 in options[1]:
				if (o1.index2, o2.index2) in bond_types_by_index2:
					return o1, o2
	
	for b in bonds:
		options = [atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))] for a in b.atoms]
		if not get_bond_type(options):
			raise Exception(b)

	for ii in range(1000):
		b = random.choice(bonds)
		if (b.atoms[0].type.index2, b.atoms[1].type.index2) not in bond_types_by_index2:
			options = [atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))] for a in b.atoms]
			if get_bond_type(options):
				b.atoms[0].type, b.atoms[1].type = get_bond_type(options)
	
	def missing_params(close_ok):
		missing_count = 0
		atom_badness = [0 for a in atoms]
		for b in bonds:
			if (b.atoms[0].type.index2, b.atoms[1].type.index2) not in bond_types_by_index2:
				missing_count += 100
				for atom in b.atoms:
					atom_badness[atom.index-1] += 1
		for a in angles:
			if tuple([atom.type.index2 for atom in a.atoms]) not in angle_types_by_index2:
				missing_count += 2
				for atom in a.atoms:
					atom_badness[atom.index-1] += 1
		for d in dihedrals:
			if tuple([a.type.index2 for a in d.atoms]) not in dihedral_types_by_index2 and not (close_ok and (0,d.atoms[1].type.index2,d.atoms[2].type.index2,0) in dihedral_types_by_index2):
				missing_count += 1
				for atom in d.atoms:
					atom_badness[atom.index-1] += 1
		return missing_count, atoms[atom_badness.index(max(atom_badness))]
	
	def anneal_params(T_max, steps, close_ok):
		best_error, worst_atom = missing_params(close_ok)
		for T in numpy.arange(T_max,0.,-T_max/steps):
			a = worst_atom if random.random()<0.8 else random.choice(atoms)
			old_type = a.type
			a.type = random.choice(atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))])
			error, worst_atom = missing_params(close_ok)
			#print (best_error-error)
			if error < best_error or random.random() < math.exp( (best_error-error)/T ):
				best_error = error
			else:
				a.type = old_type
	
	best_types = None
	best_error = missing_params(False)[0]
	for i in range(5):
		anneal_params(1.,1000,False)
		anneal_params(0.1,1000,True)
		error = missing_params(False)[0]
		if error < best_error:
			best_types = [a.type for a in atoms]
			error = best_error
	
	for i,a in enumerate(atoms):
		if best_types: a.type = best_types[i]
		print a.type.notes
	print missing_params(True)
	
	
	c1,c2,c3 = 0,0,0
	for b in bonds:
		if (b.atoms[0].type.index2, b.atoms[1].type.index2) not in bond_types_by_index2:
			c1 += 1
	for a in angles:
		if tuple([atom.type.index2 for atom in a.atoms]) not in angle_types_by_index2:
			c2 += 1
	for d in dihedrals:
		if tuple([a.type.index2 for a in d.atoms]) not in dihedral_types_by_index2:
			c3 += 1
			if (0,d.atoms[1].type.index2,d.atoms[2].type.index2,0) in dihedral_types_by_index2:
				c3 -= 1
	print c1, c2, c3
	
	
	bond_types_used = [ bond_types_by_index2[(b.atoms[0].type.index2,b.atoms[1].type.index2)] for b in bonds]
	for i,b in enumerate(bonds): b.e = bond_types_used[i].e
	
	angle_types_used = []
	for a in angles:
		index2s = tuple([atom.type.index2 for atom in a.atoms])
		if index2s in angle_types_by_index2:
			angle_types_used.append(angle_types_by_index2[index2s])
		else:
			angle_types_used.append( utils.Struct(angle=180., e=0.) )
	for i,a in enumerate(angles): a.e = angle_types_used[i].e
	
	dihedral_types_used = []
	for d in dihedrals:
		index2s = tuple([a.type.index2 for a in d.atoms])
		if index2s in dihedral_types_by_index2:
			dihedral_types_used.append(dihedral_types_by_index2[index2s])
		else:
			dihedral_types_used.append( utils.Struct(e=(0.,0.,0.,0.)) )
	for i,d in enumerate(dihedrals): d.e = dihedral_types_used[i].e	
	
	for a in atoms:
		a.charge = a.type.charge
	net_charge = sum([a.charge for a in atoms])
	overcharged_count = len( [a for a in atoms if a.charge*net_charge>0] )
	for a in atoms:
		if a.charge*net_charge>0:
			a.charge -= net_charge/overcharged_count
	
	return bond_types_used, angle_types_used, dihedral_types_used



def write_data_file(atoms, bonds, angles, dihedrals, starting_params, run_name):
	bond_types, angle_types, dihedral_types = starting_params
	atom_types = [a.type for a in atoms]
	atom_types_used = dict( [(t,True) for t in atom_types] )
	bond_types_used = dict( [(t,True) for t in bond_types] )
	angle_types_used = dict( [(t,True) for t in angle_types] )
	dihedral_types_used = dict( [(t,True) for t in dihedral_types] )
	
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types_used)] )
	bond_type_numbers = dict( [(t,i+1) for i,t in enumerate(bond_types_used)] )
	angle_type_numbers = dict( [(t,i+1) for i,t in enumerate(angle_types_used)] )
	dihedral_type_numbers = dict( [(t,i+1) for i,t in enumerate(dihedral_types_used)] )
	
	box_size = (100,100,100)
	
	os.chdir('lammps')
	f = open(run_name+'.data', 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals))+'''  dihedrals
0  impropers

'''+str(len(atom_types_used))+'''  atom types
'''+str(len(bond_types_used))+'''  bond types
'''+str(len(angle_types_used))+'''  angle types
'''+str(len(dihedral_types_used))+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%f" % (atom_type_numbers[t], t.mass) for t in atom_types_used]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%f\t%f" % (atom_type_numbers[t], t.vdw_e, t.vdw_r) for t in atom_types_used])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (bond_type_numbers[t], t.e, t.r) for t in bond_types_used]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (angle_type_numbers[t], t.e, t.angle) for t in angle_types_used]))
	try:
		if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((dihedral_type_numbers[t],)+t.e) for t in dihedral_types_used]))
	except:
		print t

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in [a.index, 1, atom_type_numbers[a.type], a.charge, a.x, a.y, a.z]] ) for a in atoms]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, bond_type_numbers[bond_types[i]], b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, angle_type_numbers[angle_types[i]]]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, dihedral_type_numbers[dihedral_types[i]]]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()
	os.chdir('..')



def write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, pair_coeffs_included=True, atom_types=None):
	if not atom_types:
		atom_types = dict( [(t.type,True) for t in atoms] ).keys() #unique set of atom types
	bond_types = dict( [(t.type,True) for t in bonds] ).keys()
	angle_types = dict( [(t.type,True) for t in angles] ).keys()
	dihedral_types = dict( [(t.type,True) for t in dihedrals] ).keys()
	
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] ) #assign unique numbers
	bond_type_numbers = dict( [(t,i+1) for i,t in enumerate(bond_types)] )
	angle_type_numbers = dict( [(t,i+1) for i,t in enumerate(angle_types)] )
	dihedral_type_numbers = dict( [(t,i+1) for i,t in enumerate(dihedral_types)] )
	
	f = open(run_name+'.data', 'w')
	f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0  impropers\n\n' % (len(atoms), len(bonds), len(angles), len(dihedrals)) )
	f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n0  improper types\n' % (len(atom_types), len(bond_types), len(angle_types), len(dihedral_types)) )
	f.write('''
 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses

'''+('\n'.join(["%d\t%f" % (atom_type_numbers[t], t.mass) for t in atom_types]))+'\n')

	if pair_coeffs_included: f.write('\nPair Coeffs\n\n'+('\n'.join(["%d\t%f\t%f" % (atom_type_numbers[t], t.vdw_e, t.vdw_r) for t in atom_types])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (bond_type_numbers[t], t.e, t.r) for t in bond_types]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (angle_type_numbers[t], t.e, t.angle) for t in angle_types]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((dihedral_type_numbers[t],)+tuple(t.e)+((0.0,) if len(t.e)==3 else ()) ) for t in dihedral_types]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in [a.index, 1, atom_type_numbers[a.type], a.charge, a.x, a.y, a.z]] ) for a in atoms]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, bond_type_numbers[b.type], b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, angle_type_numbers[a.type]]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, dihedral_type_numbers[d.type]]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()




def write_data_file2(atoms, bonds, angles, dihedrals, run_name):
	box_size = (100,100,100)
	
	os.chdir('lammps')
	f = open(run_name+'.data', 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals))+'''  dihedrals
0  impropers

'''+str(len(atoms))+'''  atom types
'''+str(len(bonds))+'''  bond types
'''+str(len(angles))+'''  angle types
'''+str(len(dihedrals))+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%f" % (atom.index, atom.type.mass) for atom in atoms]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%f\t%f" % (atom.index, atom.type.vdw_e, atom.type.vdw_r) for atom in atoms])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, bond.e, bond.d) for i,bond in enumerate(bonds)]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, angle.e, angle.theta) for i,angle in enumerate(angles)]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((i+1,)+tuple(dihedral.e)) for i,dihedral in enumerate(dihedrals) ]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in (a.index, 1, a.index, a.charge, a.x, a.y, a.z)] ) for i,a in enumerate(atoms) ]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1, b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()
	os.chdir('..')



def run_anneal(run_name, run_on='none'):
	steps = 1000
	T = 94
	P = 1.5
	os.chdir('lammps')
	f = open(run_name+".in", 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges
'''+("newton off\npackage gpu force/neigh 0 1 1" if run_on == 'gpu' else "")+'''

#pair_style	lj/cut/coul/long'''+("/gpu" if run_on == 'gpu' else "")+''' 10.0 8.0
pair_style lj/cut/coul/cut 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
#kspace_style pppm 1.0e-6
special_bonds lj/coul 0.0 0.0 0.5

read_data	'''+run_name+'''.data

thermo_style custom temp press vol pe etotal tpcpu
thermo		300#'''+str(max(100,steps/100))+'''

#dump	1 all xyz '''+str(max(10,steps/1000))+''' '''+run_name+'''.xyz
dump	1 all xyz 100 '''+run_name+'''.xyz

minimize 0.0 1.0e-8 1000 100000

velocity all create 5000 1 rot yes dist gaussian
timestep 1.0
fix anneal all nvt temp 5000 1 10
print "anneal"
run 30000
unfix anneal

minimize 0.0 1.0e-8 1000 100000

#velocity all create '''+str(T)+''' 1 rot yes dist gaussian
#timestep 2.0
#fix pressurize all npt temp '''+str(T)+''' '''+str(T)+''' 10 iso '''+str(P)+" "+str(P)+''' 100
#print "Pressurize"
#run 3000
#unfix pressurize

#velocity all create '''+str(T)+''' 1 rot yes dist gaussian
#timestep 2.0
#fix 1 all nvt temp '''+str(T)+" "+str(T)+''' 100

#print "Running"
#run	'''+str(steps))
	f.close()

	if run_on == 'cluster':
		f = open(run_name+".nbs", 'w')
		f.write('''#!/bin/bash
##NBS-nproc: 4
##NBS-speed: 3000
##NBS-queue: "batch"

$NBS_PATH/mpiexec /fs/home/jms875/bin/lammps -in %s.in > %s.log
''' % (run_name, run_name) )
		f.close()
		os.system("jsub "+run_name+".nbs")
	elif run_on == 'gpu':
		f = open(run_name+".nbs", 'w')
		f.write('''#!/bin/bash
##NBS-nproc: 2
##NBS-queue: "gpudev"

export PATH="$PATH:/opt/nvidia/cuda/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/nvidia/cuda/lib64"
$NBS_PATH/mpiexec /fs/home/jms875/bin/lammps_gpu -in %s.in > %s.log
''' % (run_name, run_name) )
		f.close()
		os.system("jsub "+run_name+".nbs")
	else:
		if os.system("lammps < "+run_name+".in") == 0:
			pass#os.system("python view.py "+run_name)
			
	os.chdir('..')
	
	
	
def run_minimize(run_name, restrained=None, restraint_value=None):
	os.chdir('lammps')
	f = open(run_name+".in", 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges

#pair_style	lj/cut/coul/long 10.0 8.0
pair_style lj/cut/coul/cut 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
#kspace_style pppm 1.0e-6
special_bonds lj/coul 0.0 0.0 0.5

read_data	'''+run_name+'''.data
'''+('fix holdit all restrain '+( {2:'bond',3:'angle',4:'dihedral'}[len(restrained.atoms)] )+' '+' '.join([str(a.index) for a in restrained.atoms])+' 2000.0 2000.0 '+str(restraint_value) if restrained else '')+'''
thermo_style custom pe
dump	1 all xyz 1000 '''+run_name+'''.xyz
fix_modify holdit energy yes
minimize 0.0 1.0e-8 1000 100000\n\n''')
	f.close()
	if os.system("lammps < "+run_name+".in") != 0:
		raise Exception('Minimize failed')
	os.chdir('..')
	

def anneal(atoms, bonds, angles, dihedrals, params, name=''):
	run_name = utils.unique_filename('lammps/', 'anneal_'+name, '.data')
	write_data_file(atoms, bonds, angles, dihedrals, params, run_name)
	run_anneal(run_name)
	jsub.wait(run_name)
	tail = subprocess.Popen('tail lammps/'+run_name+'.xyz -n '+str(len(atoms)), shell=True, stdout=subprocess.PIPE).communicate()[0]
	if 'No such file or directory' in tail:
		raise Exception('Anneal '+run_name+' failed')
	for i,line in enumerate(tail.splitlines()):
		atoms[i].x, atoms[i].y, atoms[i].z = [float(s) for s in line.split()[1:]]
	return atoms

def minimize(atoms, bonds, angles, dihedrals, params, name='', restrained=None, restraint_value=None):
	run_name = utils.unique_filename('lammps/', 'minimize_'+name, '.data')
	write_data_file(atoms, bonds, angles, dihedrals, params, run_name)
	run_minimize(run_name, restrained=restrained, restraint_value=restraint_value)
	tail = subprocess.Popen('tail lammps/'+run_name+'.xyz -n '+str(len(atoms)), shell=True, stdout=subprocess.PIPE).communicate()[0]
	if not tail:
		raise Exception('Minimize failed')
	return [[float(xyz) for xyz in line.split()[1:]] for line in tail.splitlines()]
	
def minimize2(atoms, bonds, angles, dihedrals, name='', restrained=None, restraint_value=None):
	run_name = utils.unique_filename('lammps/', 'minimize2_'+name, '.data')
	write_data_file2(atoms, bonds, angles, dihedrals, run_name)
	run_minimize(run_name, restrained=restrained, restraint_value=restraint_value)
	tail = subprocess.Popen('tail lammps/'+run_name+'.xyz -n '+str(len(atoms)), shell=True, stdout=subprocess.PIPE).communicate()[0]
	if not tail:
		raise Exception('Minimize failed')
	return [[float(xyz) for xyz in line.split()[1:]] for line in tail.splitlines()]

def get_dihedral_values(atoms, bonds, angles, dihedrals, params):
	run_name = 'dihedral_values'
	write_data_file2(atoms, bonds, angles, dihedrals, run_name)
	os.chdir('lammps')
	f = open(run_name+'.in', 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges
pair_style lj/cut/coul/cut 8.0 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
read_data '''+run_name+'''.data
compute 1 all dihedral/local phi
compute 2 all property/local dtype
dump d1 all local 1000 '''+run_name+'''.dat c_1 c_2
run 0
''')
	f.close()
	if subprocess.call('lammps < '+run_name+'.in', shell=True) != 0:
		raise Exception('Getting dihedral values failed')
	
	pairs = [ (int(line.split()[1]), float(line.split()[0])) for line in open(run_name+'.dat').readlines()[-len(dihedrals):] ]
	os.chdir('..')
	#print pairs
	return [pair[1] for pair in sorted(pairs)]

def set_internal_coordinates(atoms, bonds, angles, dihedrals, strong_dihedral=None):
	box_size = (100,100,100)
	os.chdir('lammps')
	f = open('internal_coords.data', 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals)*0)+'''	dihedrals
0  impropers

'''+str(len(atoms))+'''  atom types
'''+str(len(bonds))+'''  bond types
'''+str(len(angles))+'''  angle types
'''+str(len(dihedrals)*0)+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%f" % (atom.index, atom.type.mass) for atom in atoms]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%f\t%f" % (atom.index, atom.type.vdw_e, atom.type.vdw_r) for atom in atoms])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, max(10, bond.e), bond.d) for i,bond in enumerate(bonds)]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, max(10, angle.e), angle.theta) for i,angle in enumerate(angles)]))
	#if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d 10 1 %d 0.5" % (i+1, dihedral.theta+180) for i,dihedral in enumerate(dihedrals) ])) #energy 1 angle(degrees) weighting


	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in (a.index, 1, a.index, a.charge, a.x, a.y, a.z)] ) for i,a in enumerate(atoms) ]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1, b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	#if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()

	#'''+ '\n'.join(['fix dih'+str(i)+' all restrain dihedral '+' '.join([str(a.index) for a in d.atoms])+(' 0.0 0.0 ' if d!=strong_dihedral else ' 0.0 0.0 ')+str(d.theta+180)+'\nfix_modify dih'+str(i)+' energy yes' for i,d in enumerate(dihedrals)]) +'''

	f = open('internal_coords.in', 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges
pair_style lj/cut/coul/cut 10.0 8.0 #0.01 0.01
#pair_style lj/charmm/coul/charmm 8.0 10.0
bond_style harmonic
angle_style harmonic
#dihedral_style charmm
special_bonds lj/coul 0.0 0.0 0.5
read_data internal_coords.data

thermo_style custom epair emol etotal
dump 1 all xyz 1 internal_coords.xyz
#minimize 0.0 1.0e-8 1000 100000
#velocity all create 600.0 1337 mom yes rot yes dist gaussian
#fix NVE all nve/limit 1.0
#fix TFIX all langevin 600.0 0.0 100 1337
#run 10000
#fix TFIX all langevin 0.0 0.0 100 24601
#run 10000
minimize 0.0 1.0e-8 1000 100000
''')
	f.close()
	if subprocess.call('lammps < internal_coords.in', shell=True) != 0:
		raise Exception('Setting internal coordinates failed')
	tail = subprocess.Popen('tail internal_coords.xyz -n '+str(len(atoms)), shell=True, stdout=subprocess.PIPE).communicate()[0]
	os.chdir('..')
	return [[float(xyz) for xyz in line.split()[1:]] for line in tail.splitlines()]

def dynamics(atoms, bonds, angles, dihedrals, params, name=''):
	run_name = utils.unique_filename('lammps/', 'dynamics_'+name, '.data')
	#write_data_file(atoms, bonds, angles, dihedrals, params, run_name)
	#run_minimize(run_name, restrained=restrained, restraint_value=restraint_value)

def opls_energy(coords, atoms, bonds, angles, dihedrals, params):
	box_size = (100,100,100)
	os.chdir('lammps')
	f = open('energy.data', 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals))+'''  dihedrals
0  impropers

'''+str(len(atoms))+'''  atom types
'''+str(len(bonds))+'''  bond types
'''+str(len(angles))+'''  angle types
'''+str(len(dihedrals))+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%f" % (atom.index, atom.type.mass) for atom in atoms]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%f\t%f" % (atom.index, atom.type.vdw_e, atom.type.vdw_r) for atom in atoms])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, params.bond_e[i], bond.d) for i,bond in enumerate(bonds)]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, params.angle_e[i], angle.theta) for i,angle in enumerate(angles)]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((i+1,)+tuple(params.dihedral_e[i])) for i,dihedral in enumerate(dihedrals) ]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in (a.index, 1, a.index, a.charge)+tuple(coords[i])] ) for i,a in enumerate(atoms) ]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1, b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()

	f = open('energy.in', 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges

#pair_style	lj/cut/coul/long 10.0 8.0
pair_style lj/cut/coul/cut 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
#kspace_style pppm 1.0e-6
special_bonds lj/coul 0.0 0.0 0.5

read_data	energy.data

thermo_style custom etotal

run 0
''')
	f.close()
	result = subprocess.Popen('lammps < energy.in', shell=True, stdout=subprocess.PIPE).communicate()[0]
	energy = float( re.search('TotEng\s+(\S+)', result).group(1) )
	os.chdir('..')
	return energy



def opls_hbond_energy(coords, atoms, bonds, angles, dihedrals):
	run_name = 'energy'
	directory = 'lammps'
	box_size = (100,)*3
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = coords[i]
	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	
	hbonded = True
	polar_H_label = [a.type.index+1 for i,a in enumerate(atoms) if a.element=='H' and a.charge>0.3][0]
	electron_donor_type = [a.type for i,a in enumerate(atoms) if a.element=='N' and a.charge<-0.5][0]
	is_charged = any([a.charge!=0.0 for a in atoms])

	if not os.path.isdir(directory):
		os.mkdir(directory)
	os.chdir(directory)
	
	if hbonded:
		write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, pair_coeffs_included=False)
	else:
		write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name)
	
	f = open(run_name+'.in', 'w')
	f.write('''units	real\natom_style	full #bonds, angles, dihedrals, impropers, charges\n''')

	if hbonded:
		f.write('pair_style hybrid/overlay hbond/dreiding/morse 2 9.0 11.0 90 lj/cut/coul/cut 20.0\n')
	elif is_charged:
		f.write('pair_style lj/cut/coul/cut 20.0\n')
	else:
		f.write('pair_style lj/cut 10.0\n')
	if bonds: f.write('bond_style harmonic\n')
	if angles: f.write('angle_style harmonic\n')
	if dihedrals: f.write('dihedral_style opls\n')

	f.write('''special_bonds lj/coul 0.0 0.0 0.0\nread_data	'''+run_name+'''.data\n''')

	if hbonded:
		for ii,a in enumerate(atom_types):
			for jj,b in enumerate(atom_types[ii:]):
				i = ii+1; j = i+jj
				vdw_e = (a.vdw_e * b.vdw_e)**0.5
				vdw_r = (a.vdw_r * b.vdw_r)**0.5
				f.write('pair_coeff %d %d lj/cut/coul/cut %f %f\n' % (i, j, vdw_e, vdw_r))
				if a==electron_donor_type and b==electron_donor_type:
					#http://pubs.acs.org/doi/suppl/10.1021/ja8100227 supplementary materials
					D0 = 0.4300
					r0 = 3.4000
					alpha = 5/r0
					for donor_flag in ['i','j']: #donor = electronegative atom bonded to the H
						f.write('pair_coeff %d %d hbond/dreiding/morse %d %s %f %f %f 2\n' % (i, j, polar_H_label, donor_flag, D0, alpha, r0))

	f.write('''
compute   hb all pair hbond/dreiding/morse
variable    C_hbond equal c_hb[1] #number hbonds
variable    E_hbond equal c_hb[2] #hbond energy
thermo_style 	custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_E_hbond v_C_hbond
thermo_modify	line multi format float %14.6f
run 0
''')
	f.close()
	#os.system('lammps < '+run_name+'.in')
	result = subprocess.Popen('lammps < '+run_name+'.in', shell=True, stdout=subprocess.PIPE).communicate()[0]
	energy = float( re.search('TotEng\s+=\s+(\S+)', result).group(1) )
	os.chdir('..')
	return energy


def evaluate(mol_atoms, mol_bonds, mol_angles, mol_dihedrals, N=1):
	run_name = utils.unique_filename('lammps/', 'eval', '.data')
	
	import copy
	atoms = []
	bonds = []
	angles = []
	dihedrals = []
	for mol in range(N):
		for a in mol_atoms:
			aa = copy.copy(a)
			aa.y += mol*10.
			aa.index += mol*len(mol_atoms)
			atoms.append(aa)
		for b in mol_bonds:
			bb = copy.copy(b)
			bb.atoms = [atoms[ x.index-1+mol*len(mol_atoms) ] for x in b.atoms]
			bonds.append(bb)
		for b in mol_angles:
			bb = copy.copy(b)
			bb.atoms = [atoms[ x.index-1+mol*len(mol_atoms) ] for x in b.atoms]
			angles.append(bb)
		for b in mol_dihedrals:
			bb = copy.copy(b)
			bb.atoms = [atoms[ x.index-1+mol*len(mol_atoms) ] for x in b.atoms]
			dihedrals.append(bb)
	
	box_size = (max(atoms, key=lambda a: a.x).x*2+100, max(atoms, key=lambda a: a.y).y*2+100, max(atoms, key=lambda a: a.z).z*2+100)
	os.chdir('lammps')
	f = open(run_name+'.data', 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals))+'''  dihedrals
0  impropers

'''+str(len(atoms))+'''  atom types
'''+str(len(bonds))+'''  bond types
'''+str(len(angles))+'''  angle types
'''+str(len(dihedrals))+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%f" % (atom.index, atom.type.mass) for atom in atoms]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%f\t%f" % (atom.index, atom.type.vdw_e, atom.type.vdw_r) for atom in atoms])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, bond.e, bond.d) for i,bond in enumerate(bonds)]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (i+1, angle.e, angle.theta) for i,angle in enumerate(angles)]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((i+1,)+tuple(dihedral.e)) for i,dihedral in enumerate(dihedrals) ]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in (a.index, 1, a.index, a.charge, a.x, a.y, a.z)] ) for i,a in enumerate(atoms) ]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1, b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()

	f = open(run_name+'.in', 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges

#pair_style	lj/cut/coul/long 10.0 8.0
pair_style lj/cut/coul/cut 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
#kspace_style pppm 1.0e-6
special_bonds lj/coul 0.0 0.0 0.5

read_data	'''+run_name+'''.data

thermo_style custom temp press vol pe etotal tpcpu
thermo		300
dump	1 all xyz 100 '''+run_name+'''.xyz

minimize 0.0 1.0e-8 1000 100000

velocity all create 5000 1 rot yes dist gaussian
timestep 1.0
fix anneal all nvt temp 5000 1 10
print "anneal"
run 30000
unfix anneal

minimize 0.0 1.0e-8 1000 100000
''')
	f.close()
	os.system('lammps < '+run_name+'.in')
	os.chdir('..')

