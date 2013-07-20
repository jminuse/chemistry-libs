import math, random, sys
import utils

def dist_squared(a,b):
	return (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2

def dihedral_angle(a,b,c,d):
	dd = dist_squared
	dd12, dd23, dd34, dd24, dd13, dd14 = dd(a,b), dd(b,c), dd(c,d), dd(b,d), dd(a,c), dd(a,d)

	dd12, dd23, dd34, dd24, dd13, dd14 = [i+random.random()*1e-10 for i in (dd12, dd23, dd34, dd24, dd13, dd14)]

	P = dd12*(dd23+dd34-dd24) + dd23*(-dd23+dd34+dd24) + dd13*(dd23-dd34+dd24) - 2*dd23*dd14
	d12,d23,d34,d24,d13 = [math.sqrt(i) for i in dd12,dd23,dd34,dd24,dd13]
	Q = (d12 + d23 + d13) * ( d12 + d23 - d13) * (d12 - d23 + d13) * (-d12 + d23 + d13 ) * (d23 + d34 + d24) * ( d23 + d34 - d24 ) * (d23 - d34 + d24) * (-d23 + d34 + d24 )
	
	return 180/math.pi*math.acos(P/math.sqrt(Q))

def get_bonds(atoms):
	bonds = []
	for i,a in enumerate(atoms):
		for b in atoms[i+1:]:
			d = dist_squared(a,b)**0.5
			if (a.element!='H' and b.element!='H' and d<2.) or d < 1.2:
				bonds.append( utils.Struct(atoms=(a,b), d=d) ) #offset from current, distance
				a.bonds.append(b)
				b.bonds.append(a)
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
				theta = 180/math.pi*math.acos((A**2+B**2-N**2)/(2*A*B))
				angles.append( utils.Struct( atoms=(a,center,b), theta=theta) )
	dihedral_set = {}
	for angle in angles:
		for a in angle.atoms[0].bonded:
			if a is angle.atoms[1]: continue
			dihedral = (a,) + angle.atoms
			if reversed(dihedral) not in dihedral_set:
				dihedral_set[dihedral] = True
		
		for b in angle.atoms[2].bonded:
			if b is angle.atoms[1]: continue
			dihedral = angle.atoms + (b,)
			if reversed(dihedral) not in dihedral_set:
				dihedral_set[dihedral] = True
	dihedrals = [utils.Struct( atoms=d, theta=dihedral_angle(d[0],d[1],d[2],d[3]) ) for d in dihedral_set.keys()]
	
	return angles, dihedrals

def parse_tinker_arc(molecule_file):
	atoms = []
	for line in open(molecule_file):
		columns = line.split()
		if len(columns)>3:
			atoms.append( utils.Struct(index=int(columns[0]), element=columns[1], x=float(columns[2]), y=float(columns[3]), z=float(columns[4]), bonded=[int(s) for s in columns[6:]], type=None) )
	bond_set = {}
	for a in atoms:
		a.bonded = [atoms[i-1] for i in a.bonded]
		for b in a.bonded:
			if (b,a) not in bond_set:
				bond_set[(a,b)] = True
	bonds = [utils.Struct(atoms=b, d=dist_squared(b[0],b[1])**0.5) for b in bond_set.keys()]
	angles, dihedrals = get_angles_and_dihedrals(atoms, bonds)
	return atoms, bonds, angles, dihedrals


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

