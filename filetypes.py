import math, random, sys, copy, re
import utils, lammps

def get_bonds(atoms):
	bonds = []
	for a in atoms: a.bonded = []
	for i,a in enumerate(atoms):
		for b in atoms[i+1:]:
			d = utils.dist_squared(a,b)**0.5
			if (a.element not in [1,'H'] and b.element not in [1,'H'] and d<2.) or (d < 1.2 and (a.element in [1,'H'])!=(b.element in [1,'H']) ):
				bonds.append( utils.Struct(atoms=(a,b), d=d, e=None) ) #offset from current, distance
				a.bonded.append(b)
				b.bonded.append(a)
	return bonds

def get_angles_and_dihedrals(atoms, bonds):
	angles = [];
	for center in atoms:
		if len(center.bonded)<2: continue
		for i,a in enumerate(center.bonded):
			for b in center.bonded[i+1:]:
				A = math.sqrt((center.z-b.z)**2+(center.x-b.x)**2+(center.y-b.y)**2)
				N = math.sqrt((a.z-b.z)**2+(a.x-b.x)**2+(a.y-b.y)**2)
				B = math.sqrt((center.z-a.z)**2+(center.x-a.x)**2+(center.y-a.y)**2)
				try:
					theta = 180/math.pi*math.acos((A**2+B**2-N**2)/(2*A*B))
				except: theta = 0.0
				angles.append( utils.Struct( atoms=(a,center,b), theta=theta, e=None, type=None) )
	dihedral_set = {}
	for angle in angles:
		for a in angle.atoms[0].bonded:
			if a is angle.atoms[1]: continue
			dihedral = (a,) + angle.atoms
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True
		
		for b in angle.atoms[2].bonded:
			if b is angle.atoms[1]: continue
			dihedral = angle.atoms + (b,)
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True
	dihedrals = [utils.Struct( atoms=d, theta=None, e=None, type=None) for d in dihedral_set.keys()]
	
	return angles, dihedrals

def parse_tinker_arc(molecule_file, parameter_file=None):
	if parameter_file:
		elements, atom_types, bond_types, angle_types, dihedral_types = lammps.parse_opls_parameters(parameter_file)

	atoms = []
	for line in open(molecule_file):
		columns = line.split()
		if len(columns)>3:
			atoms.append( utils.Struct(index=int(columns[0]), element=columns[1], x=float(columns[2]), y=float(columns[3]), z=float(columns[4]), bonded=[int(s) for s in columns[6:]], type=([t for t in atom_types if t.index==int(columns[5][-3:])][0] if parameter_file else None), charge=None) )
			if len(columns[5])>3:
				atom_types.append( copy.deepcopy(atoms[-1].type) )
				atoms[-1].type.index = int(columns[5])
	bond_set = {}
	for a in atoms:
		a.bonded = [atoms[i-1] for i in a.bonded]
		for b in a.bonded:
			if (b,a) not in bond_set:
				bond_set[(a,b)] = True
	bonds = [utils.Struct(atoms=b, d=utils.dist_squared(b[0],b[1])**0.5, e=None, type=None) for b in bond_set.keys()]
	angles, dihedrals = get_angles_and_dihedrals(atoms, bonds)
	
	return atoms, bonds, angles, dihedrals

def parse_sdf(molecule_file):
	atoms = []
	for line in open(molecule_file):
		columns = line.split()
		if len(columns)==16: #atom line
			x,y,z = [float(s) for s in columns[:3]]
			element = columns[3]
			atoms.append( utils.Struct(index=len(atoms)+1, element=element, x=x, y=y, z=z, bonded=[], type=utils.Struct(index=0)) )
		elif len(columns)==7: #bond line
			i1, i2 = [int(s) for s in columns[:2]]
			atoms[i1-1].bonded.append( atoms[i2-1] )
			atoms[i2-1].bonded.append( atoms[i1-1] )
	return atoms

def compare_structures(molecule_files):
	atoms, bonds, angles, dihedrals = [],[],[],[]
	for molecule_file in molecule_files:
		a,b,c,d = parse_tinker(molecule_file)
		atoms.append(a); bonds.append(b); angles.append(c); dihedrals.append(d)

	print 'Bonds'
	for i in range(len(bonds[0])):
		error = (bonds[0][i][2]-bonds[1][i][2])/bonds[1][i][2]
		if abs(error)>0.05:
			print ("%6d"*2 + "%10.3f"*3) % (bonds[0][i][:2]+(bonds[0][i][2], bonds[1][i][2], error))

	print 'Angles'
	for i in range(len(angles[0])):
		error = (angles[0][i][3]-angles[1][i][3])
		if angles[0][i][:3]!=angles[1][i][:3]:
			print angles[1][i][:3]
		if abs(error)>5:
			print ("%6d"*3 + "%10.3f"*3) % (angles[0][i][:3]+(angles[0][i][3], angles[1][i][3], error))

	print 'Dihedrals'
	for i in range(len(dihedrals[0])):
		error = (dihedrals[0][i][4]-dihedrals[1][i][4])
		if abs(error)>5:
			print ("%6d"*4 + "%10.3f"*3) % (dihedrals[0][i][:4]+(dihedrals[0][i][4], dihedrals[1][i][4], error))

def write_xyz(name, atoms, f=None):
	if not f:
		f = open(name+'.xyz', 'w')
		f.write(str(len(atoms))+'\nAtoms\n')
		for a in atoms:
			f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		f.close()
	else:
		f.write(str(len(atoms))+'\nAtoms\n')
		for a in atoms:
			f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
			
def parse_xyz(name):
	atoms = []
	for line in open(name):
		columns = line.split()
		if len(columns)>=4:
			x,y,z = [float(s) for s in columns[1:4]]
			atoms.append( utils.Struct(element=columns[0], x=x,y=y,z=z) )
	return atoms
	
def write_arc(name, atoms, index_offset=0):
	f = open(name+'.arc', 'w')
	f.write(str(len(atoms))+' Atoms\n')
	for a in atoms:
		f.write('%d %s %f %f %f %d ' % (a.index, a.element, a.x, a.y, a.z, a.type.index+index_offset) + (' '.join([str(a.index) for a in a.bonded]))+'\n' )
	f.close()

def gaussian_to_xyz(input_file, output_file):
	f = open(output_file, 'w')
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		print 'Job did not finish'
	else:
		m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
		time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))
		print m.group(0)
		print "%.2e s" % time

	#a = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
	#print a
	print '\n'.join(re.findall('SCF Done: +\S+ += +(\S+)', contents))
	
	if 'Counterpoise: corrected energy =' in contents:
		print 'Counterpoise E = ',
		print '\n'.join(re.findall('Counterpoise: corrected energy = +(\S+)', contents))

	start = 0
	while True:
		try:
			next_coordinates = contents.index('Coordinates (Angstroms)', start)
		except: break
		start = contents.index('---\n', next_coordinates)+4
		end = contents.index('\n ---', start)
		lines = contents[start:end].splitlines()
		start = end

		f.write("%d\nAtoms\n" % len(lines))

		for line in lines:
			columns = line.split()
			element = columns[1]
			x,y,z = columns[3:6]
			f.write( '\t'.join([element, x, y, z]) + '\n' )

	f.close()
	import os
	os.system('/fs/europa/g_pc/vmd-1.9 '+output_file+' > /dev/null')

