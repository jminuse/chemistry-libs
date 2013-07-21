import re, random, numpy, math
import utils

def parse_parameter_file(parameter_file):
	atom_types = []; bond_types = []; angle_types = []; dihedral_types = []; hbond_types = []; vdw_pair_types = []
	for line in open(parameter_file):
		columns = line.split()
		if not columns: continue
		if columns[0]=='atom':
			m = re.match('atom\s+(\S+)\s+(\S+)\s+"([^"]+)"\s+(\S+)\s+(\S+)\s+(\S+)', line)
			atom_types.append( utils.Struct(index=int(m.group(1)), element_name=m.group(2), notes=m.group(3), element=int(m.group(4)), mass=float(m.group(5)), bond_count=int(m.group(6)) ) )
		elif columns[0]=='vdw':
			atom_types[int(columns[1])-1].vdw_r = float(columns[2])
			atom_types[int(columns[1])-1].vdw_e = float(columns[3])
		elif columns[0]=='vdwpr':
			vdw_pair_types.append( utils.Struct(indices=tuple([int(s) for s in columns[1:3]]),vdw_r=float(columns[3]),vdw_e=float(columns[4])) )
		elif columns[0]=='hbond':
			hbond_types.append( utils.Struct(indices=tuple([int(s) for s in columns[1:3]]),r=float(columns[3]),e=float(columns[4])) )
		elif columns[0]=='bond':
			bond_types.append( utils.Struct(indices=tuple([int(s) for s in columns[1:3]]),e=float(columns[3]),r=float(columns[4])) )
		elif columns[0]=='angle':
			angle_types.append( utils.Struct(indices=tuple([int(s) for s in columns[1:4]]),e=float(columns[4]),angle=float(columns[5])) )
		elif columns[0]=='torsion':
			dihedral_types.append( utils.Struct(indices=tuple([int(s) for s in columns[1:5]]),e=columns[5:]) )
	
	return atom_types, bond_types, angle_types, dihedral_types, hbond_types, vdw_pair_types

def guess_parameters(atoms, bonds, angles, dihedrals, parameter_file):
	atom_types, bond_types, angle_types, dihedral_types, hbond_types, vdw_pair_types = parse_parameter_file(parameter_file)
	
	atom_types_by_element_and_bond_count = {}
	for a in atom_types:
		if (a.element, a.bond_count) not in atom_types_by_element_and_bond_count:
			atom_types_by_element_and_bond_count[(a.element, a.bond_count)] = []
		atom_types_by_element_and_bond_count[(a.element, a.bond_count)].append(a)

	atomic_number = {'H':1, 'C':6, 'N':7}
	for a in atoms:
		a.type = atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))][0]
	
	dihedral_types_by_index = dict([ (t.indices,t) for t in dihedral_types] + [ (t.indices[::-1],t) for t in dihedral_types])
	angle_types_by_index = dict([ (t.indices,t) for t in angle_types] + [ (t.indices[::-1],t) for t in angle_types])
	bond_types_by_index = dict([ (t.indices,t) for t in bond_types] + [ (t.indices[::-1],t) for t in bond_types])

	def missing_params(close_ok):
		missing_count = 0
		for b in bonds:
			if (b.atoms[0].type.index, b.atoms[1].type.index) not in bond_types_by_index:
				missing_count += 100
		for d in dihedrals:
			if tuple([a.type.index for a in d.atoms]) not in dihedral_types_by_index:
				missing_count += 1
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
	
	#anneal_params(1.,1000,False)
	
	for a in atoms:
		print a.type.notes
	print missing_params(True)
	
	missing_count = 0
	for b in bonds:
		if (b.atoms[0].type.index, b.atoms[1].type.index) not in bond_types_by_index:
			missing_count += 1
	print missing_count
	
	def get_bond_type(options):
		for o1 in options[0]:
			for o2 in options[1]:
				if (o1.index, o2.index) in bond_types_by_index:
					return o1, o2
	
	#for b in bonds:
	#	options = [atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))] for a in b.atoms]
	#	if not get_bond_type(options):
	#		print [a.index for a in b.atoms]
	
	for ii in range(len(bonds)*10):
		b = random.choice(bonds)
		if (b.atoms[0].type.index, b.atoms[1].type.index) not in bond_types_by_index:
			options = [atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))] for a in b.atoms]
			if get_bond_type(options):
				b.atoms[0].type, b.atoms[1].type = get_bond_type(options)
	
	for a in atoms:
		print a.type.notes
	print missing_params(True)
	
	#atom_types_used = dict([(a.type,True) for a in atoms]).keys()
	#bond_types_used = [ bond_types_by_index[(b.atoms[0].type.index,b.atoms[1].type.index)] for b in bonds]
	#angle_types_used = [ angle_types_by_index[(a.atoms[0].type.index,a.atoms[1].type.index,a.atoms[2].type.index)] for a in angles]
	#dihedral_types
	
	#return atom_types_used, bond_types_used, angle_types_used, dihedral_types_used
	
