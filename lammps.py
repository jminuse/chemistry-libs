import re, random, numpy, math, os, subprocess
import utils

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
			atom_types[int(columns[1])-1].vdw_r = float(columns[2])
			atom_types[int(columns[1])-1].vdw_e = float(columns[3])
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
	
	angle_types_used = []
	for a in angles:
		index2s = tuple([atom.type.index2 for atom in a.atoms])
		if index2s in angle_types_by_index2:
			angle_types_used.append(angle_types_by_index2[index2s])
		else:
			angle_types_used.append( utils.Struct(angle=180., e=0.) )
	
	dihedral_types_used = []
	for d in dihedrals:
		index2s = tuple([a.type.index2 for a in d.atoms])
		if index2s in dihedral_types_by_index2:
			dihedral_types_used.append(dihedral_types_by_index2[index2s])
		else:
			dihedral_types_used.append( utils.Struct(e=(0.,0.,0.,0.)) )
			
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
	
	f = open(run_name+".data", 'w')
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




def run(run_name, run_type='anneal', run_on='none'):
	steps = 1000
	T = 94
	P = 1.5
	f = open(run_name+".in", 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges
'''+("newton off\npackage gpu force/neigh 0 1 1" if run_on == 'gpu' else "")+'''

pair_style	lj/cut/coul/long'''+("/gpu" if run_on == 'gpu' else "")+''' 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
kspace_style pppm 1.0e-6
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
			os.system("python view.py "+run_name)

def anneal(atoms, bonds, angles, dihedrals, starting_params):
	run_number = 0
	while True:
		run_name = 'anneal_'+str(run_number)
		if not os.path.exists(run_name+'.data'): break
		run_number += 1
	write_data_file(atoms, bonds, angles, dihedrals, starting_params, run_name)
	run(run_name)
	while(True):
		jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
		if run_name in jlist:
			os.sleep(60)
		else:
			break
	for i,line in enumerate(subprocess.Popen('tail '+run_name+'.xyz -n '+str(len(atoms)), shell=True, stdout=subprocess.PIPE).communicate()[0].splitlines()):
		atoms[i].x, atoms[i].y, atoms[i].z = [float(s) for s in line.split()[1:]]
	return atoms
	
