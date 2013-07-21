import re, random, numpy, math
import utils

def parse_opls_parameters(parameter_file):
	elements = {}; atom_types = []; bond_types = []; angle_types = []; dihedral_types = []
	for line in open(parameter_file):
		columns = line.split()
		if not columns: continue
		if columns[0]=='atom':
			m = re.match('atom +(\d+) +(\d+) +(\S+) +"([^"]+)" +(\d+) +(\S+) +(\d+)', line)
			atom_type = utils.Struct(index=int(m.group(1)), index2=int(m.group(2)), element_name=m.group(3), notes=m.group(4), element=int(m.group(5)), mass=float(m.group(6)), bond_count=int(m.group(7)) )
			if atom_type.element not in elements: elements[atom_type.element] = []
			elements[atom_type.element].append(atom_type)
			if '(UA)' in atom_type.notes:
				atom_type.element = 0 #reject united-atom parameters
			atom_types.append(atom_type)
		elif columns[0]=='vdw':
			atom_types[int(columns[1])-1].vdw_r = float(columns[2])
			atom_types[int(columns[1])-1].vdw_e = float(columns[3])
		elif columns[0]=='charge':
			atom_types[int(columns[1])-1].charge = float(columns[2])
		elif columns[0]=='bond':
			bond_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:3]]),r=float(columns[3]),e=float(columns[3])) )
		elif columns[0]=='angle':
			angle_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:4]]),e=float(columns[4]),angle=float(columns[5])) )
		elif columns[0]=='torsion':
			dihedral_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:5]]),e=columns[5:]) )
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
		for b in bonds:
			if (b.atoms[0].type.index2, b.atoms[1].type.index2) not in bond_types_by_index2:
				missing_count += 100
		for a in angles:
			if tuple([atom.type.index2 for atom in a.atoms]) not in angle_types_by_index2:
				missing_count += 2
		for d in dihedrals:
			if tuple([a.type.index2 for a in d.atoms]) not in dihedral_types_by_index2:
				missing_count += 1
				if close_ok and (0,d.atoms[1].type.index2,d.atoms[2].type.index2,0) in dihedral_types_by_index2:
					missing_count -= 1
		return missing_count
	
	def anneal_params(T_max, steps, close_ok):
		best_error = missing_params(close_ok)
		for T in numpy.arange(T_max,0.,-T_max/steps):
			a = random.choice(atoms)
			old_type = a.type
			a.type = random.choice(atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))])
			error = missing_params(close_ok)
			#print (best_error-error)
			if error < best_error or random.random() < math.exp( (best_error-error)/T ):
				best_error = error
				#print error
			else:
				a.type = old_type
	
	anneal_params(1.,10000,False)
	anneal_params(0.1,10000,True)
	
	for a in atoms:
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
	
	#for a in angles:
	angle_types_used = [ angle_types_by_index2[(a.atoms[0].type.index2,a.atoms[1].type.index2,a.atoms[2].type.index2)] for a in angles]
	
	dihedral_types_used = []
	for d in dihedrals:
		index2s = tuple([a.type.index2 for a in d.atoms])
		if index2s in dihedral_types_by_index2:
			dihedrals.append(dihedral_types_by_index2[index2s])
		else:
			dihedrals.append(None)
	
	return bond_types_used, angle_types_used, dihedral_types_used


"""
def write_data_file(atoms, bonds, angles, dihedrals, starting_params, run_name):
	elements, atom_types, bond_types, angle_types, dihedral_types = starting_params
	
	f = open(run_name+".data", 'w')
	f.write('''LAMMPS Description

'''+str(len(atoms))+'''  atoms
'''+str(len(bonds))+'''  bonds
'''+str(len(angles))+'''  angles
'''+str(len(dihedrals))+'''  dihedrals
0  impropers

'''+str(len(atom_types_used.keys()))+'''  atom types
'''+str(len())+'''  bond types
'''+str(len())+'''  angle types
'''+str(len())+'''  dihedral types
0  improper types

 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses			

'''+('\n'.join(["%d\t%s" % (index, opls_atoms[number-1][0]) for number,index in sorted(atom_indices.items(), key=lambda x: x[1])]))+'''

Pair Coeffs

'''+('\n'.join(["%d\t%s\t%s" % (index, opls_atoms[number-1][3], opls_atoms[number-1][2]) for number,index in sorted(atom_indices.items(), key=lambda x: x[1])])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%s\t%s" % data for ID,data in sorted(bond_indices.items(), key=lambda x: x[1][0])]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%s\t%s" % data for ID,data in sorted(angle_indices.items(), key=lambda x: x[1][0])]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%s\t%s\t%s\t%s" % data for ID,data in sorted(dihedral_indices.items(), key=lambda x: x[1][0])]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in [i+1]+list(a[:-1])] ) for i,a in enumerate(atoms)]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1]+list(b)]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1]+list(a)]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1]+list(d)]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	if impropers: f.write('\n\nImpropers\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1]+list(p)]) for i,p in enumerate(impropers)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()


def write_input_file():
	f = open(run_name+".in", 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges
'''+("newton off\npackage gpu force/neigh 0 1 1" if run_on == 'gpu' else "")+'''

pair_style	lj/cut/coul/long'''+("/gpu" if run_on == 'gpu' else "")+''' 9 8
bond_style harmonic
angle_style harmonic
dihedral_style opls
kspace_style pppm 1.0e-4
special_bonds lj/coul 0.0 0.0 0.5

read_data	'''+run_name+'''.data

#thermo_style custom evdwl ecoul ebond eangle edihed etotal 
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
			os.system("python view.py "+run_name)

def anneal(atoms, bonds, angles, dihedrals, starting_params):
	write_data_file(atoms, bonds, angles, dihedrals, starting_params)
	write_input_file()
	run()
"""

