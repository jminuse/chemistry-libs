import re, random
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
	
	dihedral_types_by_index2 = dict([ (t.index2s,t) for t in dihedral_types])
	#index2s = dict([ for a in atoms])
	#print dihedral_types_by_index2.keys()
	
	def missing_dihedrals():
		missing_count = 0
		for d in dihedrals:
			#print tuple([a.type.index2 for a in d.atoms])
			if tuple([a.type.index2 for a in d.atoms]) not in dihedral_types_by_index2:
				missing_count += 1
				#if (0,d.atoms[1].type.index2,d.atoms[2].type.index2,0) in dihedral_types_by_index2:
				#	missing_count -= 1
		return missing_count
	
	taboo_list = [0 for a in atoms]
	best_error = missing_dihedrals()
	for i in range(100):
		best_error_this_time = len(dihedrals)*10
		best_i = 0
		best_type = None
		for i,a in enumerate(atoms):
			old_type = a.type
			a.type = random.choice(atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))])
			error = missing_dihedrals()
			if error < best_error_this_time:
				best_error_this_time = error
				best_i = i
				best_type = a.type
			if taboo_list[i]>0:
				taboo_list[i] -= 1
			a.type = old_type
		
		if taboo_list[best_i]==0 or best_error_this_time < best_error:
			taboo_list[best_i] = 5
			if best_error_this_time < best_error: best_error = best_error_this_time
			atoms[best_i].type = best_type
		else:
			best_error_this_time = len(dihedrals)*10
			best_i = 0
			best_type = None
			for i,a in enumerate(atoms):
				if taboo_list[i]>0: continue
				old_type = a.type
				a.type = random.choice(atom_types_by_element_and_bond_count[(atomic_number[a.element],len(a.bonded))])
				error = missing_dihedrals()
				if error < best_error_this_time:
					best_error_this_time = error
					best_i = i
					best_type = a.type
				a.type = old_type
			taboo_list[best_i] = 5
			atoms[best_i].type = best_type
		print best_error_this_time, best_error
	
	'''
	atoms_by_index2 = {}
	for a in atom_types:
		if a.index2 not in atoms_by_index2: atoms_by_index2[a.index2] = []
		atoms_by_index2[a.index2].append(a)
	for index2, atom_list in atoms_by_index2: #check that all atoms of the same index2 have the same element, bond_count
		element, bond_count = atom_list[0].element, atom_list[0].bond_count
		for a in atom_list[1:]:
			if a.element != element and a.bond_count != bond_count:
				raise Exception('Atoms with the same index2 do not have the same element, bond_count', atom_list[0], a)
	dihedral_types_by_atom_type = {}
	for d in dihedral_types:
		#assumes all atoms of the same index2 have the same element, bond_count
		types = tuple([(atoms_by_index2[index2][0].element, atoms_by_index2[index2][0].bond_count) for index2 in d.atoms])
		if types not in dihedral_types_by_atom_type: dihedral_types_by_atom_type[types] = []
		if reversed(types) not in dihedral_types_by_atom_type: dihedral_types_by_atom_type[reversed(types)] = []
		dihedral_types_by_atom_type[types].append(d)
		dihedral_types_by_atom_type[reversed(types)].append(d)
	
	atom_types_by_element_and_bond_count = dict([((a.element, a.bond_count), atom_type) for a in atom_types])
	atom_type_options = [ atom_types_by_element_and_bond_count[(atom.element, atom.bond_count)] for atom in atoms ]
	dihedral_type_options = []
	for d in dihedrals:
		#assumes all atoms of the same index2 have the same element, bond_count
		types = tuple([atoms_by_index2[index2].element, atoms_by_index2[index2].bond_count for index2 in d.atoms])
		options = dihedral_types_by_atom_type[types]
		#rule out conflicting atom types
		
		
		dihedral_type_options.append( options )
	'''

