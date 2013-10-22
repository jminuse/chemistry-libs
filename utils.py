import os, math, numpy, copy
import lammps, filetypes

class Memoize:
	def __init__ (self, f):
		self.f = f
		self.mem = {}
	def __call__ (self, *args, **kwargs):
		if (args, str(kwargs)) in self.mem:
			return self.mem[args, str(kwargs)]
		else:
			tmp = self.f(*args, **kwargs)
			self.mem[args, str(kwargs)] = tmp
			return tmp

class Struct:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
	def __repr__(self):
		return str( dict([ (a,None) if type(self.__dict__[a]) in (list,dict) else (a,self.__dict__[a]) for a in self.__dict__]) )


def frange(low, high, step):
	while low < high:
		yield low
		low += step
'''
def quat_to_mat(quat):	#quat = [w i j k]
	w2 = quat[0]*quat[0]
	i2 = quat[1]*quat[1]
	j2 = quat[2]*quat[2]
	k2 = quat[3]*quat[3]
	twoij = 2.0*quat[1]*quat[2]
	twoik = 2.0*quat[1]*quat[3]
	twojk = 2.0*quat[2]*quat[3]
	twoiw = 2.0*quat[1]*quat[0]
	twojw = 2.0*quat[2]*quat[0]
	twokw = 2.0*quat[3]*quat[0]

	return [(w2+i2-j2-k2, twoij-twokw, twojw+twoik),
	(twoij+twokw, w2-i2+j2-k2, twojk-twoiw),
	(twoik-twojw, twojk+twoiw, w2-i2-j2+k2)]
'''
'''
def quat_to_mat(q): #quat = [w i j k]
	q[0],q[3] = q[3],q[0]
	return [ [1.0 - 2*q[2]*q[2] - 2*q[3]*q[3], 2*q[1]*q[2] - 2*q[3]*q[0], 2*q[1]*q[3] + 2*q[2]*q[0]],
	[2*q[1]*q[2] + 2*q[3]*q[0], 1.0 - 2*q[1]*q[1] - 2*q[3]*q[3], 2*q[2]*q[3] - 2*q[1]*q[0]],
	[2*q[1]*q[3] - 2*q[2]*q[0], 2*q[2]*q[3] + 2*q[1]*q[0], 1.0 - 2*q[1]*q[1] - 2*q[2]*q[2]] ]
'''
def quat_to_mat(q): #quat = [w i j k]
	d,b,c,a = q
	return [ [a**2+b**2-c**2-d**2, 2*b*c-2*a*d, 2*b*d+2*a*c],
	[2*b*c+2*a*d, a**2-b**2+c**2-d**2, 2*c*d-2*a*b],
	[2*b*d-2*a*c, 2*c*d-2*a*b, a**2-b**2-c**2+d**2]
	]


def transform_difference(P, Q): #two sets of points on a rigid body
	import numpy
	center1, center2 = [[(max(pp, key=lambda p:p[i])[i]+min(pp, key=lambda p:p[i])[i])/2 for i in range(3)] for pp in (P, Q)]
	#P = [[p[i]-center1[i] for i in range(3)] for p in P]
	#Q = [[p[i]-center2[i] for i in range(3)] for p in Q]
	
	for s in (P,Q):
		s.append( numpy.add(s[0], numpy.cross(numpy.subtract(s[1],s[0]), numpy.subtract(s[2],s[0])) ) )
	
	PP, QQ = [numpy.transpose([ numpy.subtract(s[1],s[0]), numpy.subtract(s[2],s[0]), numpy.subtract(s[3],s[0]) ]) for s in (P,Q)]
	
	return lambda v: numpy.dot(numpy.dot(QQ,numpy.linalg.inv(PP)), v) + Q[0] - numpy.dot(numpy.dot(QQ,numpy.linalg.inv(PP)), P[0])
	

def matvec(m,v):
	return (m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2], m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2], m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2])
	
def rand_rotation(): #http://tog.acm.org/resources/GraphicsGems/, Ed III
	import random
	x = (random.random(), random.random(), random.random())
	theta = x[0] * 2*math.pi
	phi   = x[1] * 2*math.pi
	z = x[2] * 2
	#Compute a vector V used for distributing points over the sphere via the reflection I - V Transpose(V).  This formulation of V will guarantee that if x[1] and x[2] are uniformly distributed, the reflected points will be uniform on the sphere.  Note that V has length sqrt(2) to eliminate the 2 in the Householder matrix.
	r = math.sqrt(z)
	Vx = math.sin( phi ) * r
	Vy = math.cos( phi ) * r
	Vz = math.sqrt( 2.0 - z )
	#Compute the row vector S = Transpose(V) * R, where R is a simple rotation by theta about the z-axis.  No need to compute Sz since it's just Vz.
	st = math.sin( theta )
	ct = math.cos( theta )
	Sx = Vx * ct - Vy * st
	Sy = Vx * st + Vy * ct

	#Construct the rotation matrix  ( V Transpose(V) - I ) R, which is equivalent to V S - R.

	M = [ [0.,0.,0.], [0.,0.,0.], [0.,0.,0.] ]
	
	M[0][0] = Vx * Sx - ct
	M[0][1] = Vx * Sy - st
	M[0][2] = Vx * Vz

	M[1][0] = Vy * Sx + st
	M[1][1] = Vy * Sy - ct
	M[1][2] = Vy * Vz

	M[2][0] = Vz * Sx
	M[2][1] = Vz * Sy
	M[2][2] = 1.0 - z	# This equals Vz * Vz - 1.0
	
	return M

elements_by_atomic_number = ['','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo']

class Molecule():
	elements, atom_types, bond_types, angle_types, dihedral_types = None,None,None,None,None
	@staticmethod
	def set_params(param_file):
		Molecule.elements, Molecule.atom_types, Molecule.bond_types, Molecule.angle_types, Molecule.dihedral_types = lammps.parse_opls_parameters(param_file)
	
	def __init__(self, atoms_or_filename_or_all, bonds=None, angles=None, dihedrals=None): #set atoms, bonds, etc, or assume 'atoms' contains all those things if only one parameter is passed in
		if type(atoms_or_filename_or_all)==type('string'):
			filename = atoms_or_filename_or_all
			atoms = []
			for line in open(filename):
				columns = line.split()
				if len(columns)>3:
					atoms.append( Struct(index=int(columns[0]), element=columns[1], x=float(columns[2]), y=float(columns[3]), z=float(columns[4]), bonded=[int(s) for s in columns[6:]], type=[t for t in Molecule.atom_types if t.index==int(columns[5])][0], charge=None) )
			bond_set = {}
			for a in atoms:
				a.bonded = [atoms[i-1] for i in a.bonded]
				for b in a.bonded:
					if (b,a) not in bond_set:
						bond_set[(a,b)] = True
			bonds = [Struct(atoms=b, d=dist_squared(b[0],b[1])**0.5, e=None, type=None) for b in bond_set.keys()]
			angles, dihedrals = filetypes.get_angles_and_dihedrals(atoms, bonds)
		elif not bonds:
			atoms, bonds, angles, dihedrals = atoms_or_filename_or_all
		else:
			atoms = atoms_or_filename_or_all
		self.atoms = atoms
		self.bonds = bonds
		self.angles = angles
		self.dihedrals = dihedrals
		average_position = [(max(atoms, key=lambda a:a.x).x+min(atoms, key=lambda a:a.x).x)/2, (max(atoms, key=lambda a:a.y).y+min(atoms, key=lambda a:a.y).y)/2, (max(atoms, key=lambda a:a.z).z+min(atoms, key=lambda a:a.z).z)/2]
		for a in atoms: #center atoms
			a.x -= average_position[0]
			a.y -= average_position[1]
			a.z -= average_position[2]
		
		if type(atoms_or_filename_or_all)==type('string'):
			self.set_types(Molecule.bond_types, Molecule.angle_types, Molecule.dihedral_types)

	def set_types(self, bond_types, angle_types, dihedral_types): #given the atom types, find all the other types using the provided lists
		net_charge = sum([x.type.charge for x in self.atoms])
		count_positive = len([x for x in self.atoms if x.type.charge>0.0])
		charge_adjustment = net_charge/count_positive if net_charge>0.0 else net_charge/(len(self.atoms)-count_positive)
		for x in self.atoms:
			x.charge = x.type.charge-charge_adjustment if (x.type.charge>0.0)==(net_charge>0.0) else x.type.charge
		for x in self.bonds+self.angles+self.dihedrals:
			index2s = tuple([a.type.index2 for a in x.atoms])
			try:
				x.type = [t for t in bond_types+angle_types+dihedral_types if t.index2s==index2s or t.index2s==tuple(reversed(index2s))][0]
			except:
				print 'No params for', index2s
				if x in self.bonds: x.type = bond_types[0]
				if x in self.angles: x.type = angle_types[0]; x.type.e = 0.0
				if x in self.dihedrals: x.type = dihedral_types[0]; x.type.e = [0.0]*3
	def add_to(self, x, y, z, atoms, bonds, angles, dihedrals): #add an instance of this molecule to the provided lists at the provided position
		atom_offset = len(atoms)
		for a in self.atoms:
			atoms.append( Struct(index=a.index+atom_offset, element=a.element, x=a.x+x, y=a.y+y, z=a.z+z, bonded=None, type=a.type, charge=a.charge) )
		for t in self.bonds:
			bonds.append( Struct(atoms=tuple([atoms[a.index+atom_offset-1] for a in t.atoms]), type=t.type) )
		for t in self.angles:
			angles.append( Struct(atoms=tuple([atoms[a.index+atom_offset-1] for a in t.atoms]), type=t.type) )
		for t in self.dihedrals:
			dihedrals.append( Struct(atoms=tuple([atoms[a.index+atom_offset-1] for a in t.atoms]), type=t.type) )
	def rotate(self, m):
		for a in self.atoms:
			a.x, a.y, a.z = matvec(m, (a.x, a.y, a.z))
	def rand_rotate(self):
		rand_m = rand_rotation()
		self.rotate(rand_m)

def unique_filename(directory, name, filetype):
	number = 0
	while True:
		unique_name = name+'_'+str(number)
		number += 1
		if not os.path.exists(directory+unique_name+filetype): return unique_name

#@Memoize	
def dist_squared(a,b):
	return (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2

#@Memoize
def dist(a,b):
	return dist_squared(a,b)**0.5

#@Memoize
def angle_size(a,center,b):
	A = math.sqrt((center.z-b.z)**2+(center.x-b.x)**2+(center.y-b.y)**2)
	N = math.sqrt((a.z-b.z)**2+(a.x-b.x)**2+(a.y-b.y)**2)
	B = math.sqrt((center.z-a.z)**2+(center.x-a.x)**2+(center.y-a.y)**2)
	return 180/math.pi*math.acos((A**2+B**2-N**2)/(2*A*B))

#@Memoize
dihedral_angle_cache = {}
def dihedral_angle(a,b,c,d):
	cache_key = a.x+b.x+c.x+d.x
	#if cache_key in dihedral_angle_cache:
	#	return dihedral_angle_cache[cache_key]
	
	sqrt = math.sqrt
	cos = math.cos
	
	#a,b,c,d = d,c,b,a
	
	vb1x, vb1y, vb1z = a.x-b.x, a.y-b.y, a.z-b.z
	vb2x, vb2y, vb2z = c.x-b.x, c.y-b.y, c.z-b.z
	vb3x, vb3y, vb3z = d.x-c.x, d.y-c.y, d.z-c.z
	
	# c0 calculation

	sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
	sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
	sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

	rb1 = sqrt(sb1)
	rb3 = sqrt(sb3)

	c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3

	# 1st and 2nd angle

	b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
	b1mag = sqrt(b1mag2)
	b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
	b2mag = sqrt(b2mag2)
	b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
	b3mag = sqrt(b3mag2)

	ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
	r12c1 = 1.0 / (b1mag*b2mag)
	c1mag = ctmp * r12c1

	ctmp = vb2x*vb3x + vb2y*vb3y + vb2z*vb3z
	r12c2 = 1.0 / (b2mag*b3mag)
	c2mag = -ctmp * r12c2

	# cos and sin of 2 angles and final c

	sin2 = max(1.0 - c1mag*c1mag,0.0)
	sc1 = sqrt(sin2)
	#if sc1 < SMALL: sc1 = SMALL
	sc1 = 1.0/sc1

	sin2 = max(1.0 - c2mag*c2mag,0.0)
	sc2 = sqrt(sin2)
	#if sc2 < SMALL: sc2 = SMALL
	sc2 = 1.0/sc2

	s1 = sc1 * sc1
	s2 = sc2 * sc2
	s12 = sc1 * sc2
	c = (c0 + c1mag*c2mag) * s12
	
	cx = vb1y*vb2z - vb1z*vb2y
	cy = vb1z*vb2x - vb1x*vb2z
	cz = vb1x*vb2y - vb1y*vb2x
	cmag = sqrt(cx*cx + cy*cy + cz*cz)
	dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag
	
	
	if c>1.0: c = 1.0
	if c<-1.0: c = -1.0
	
	phi = math.acos(c)
	if dx < 0.0:
		phi *= -1.0
	phi *= -1.0
	
	dihedral_angle_cache[cache_key] = phi, math.cos(phi), math.cos(2*phi), math.cos(3*phi)
	
	return phi, math.cos(phi), math.cos(2*phi), math.cos(3*phi), math.cos(4*phi)
	
	'''
	vb1x = a.x - b.x
	vb1y = a.y - b.y
	vb1z = a.z - b.z
	#domain->minimum_image(vb1x,vb1y,vb1z)

	vb2x = c.x - b.x
	vb2y = c.y - b.y
	vb2z = c.z - b.z
	#domain->minimum_image(vb2x,vb2y,vb2z)

	vb2xm = -vb2x
	vb2ym = -vb2y
	vb2zm = -vb2z
	#domain->minimum_image(vb2xm,vb2ym,vb2zm)

	vb3x = d.x - c.x
	vb3y = d.y - c.y
	vb3z = d.z - c.z
	#domain->minimum_image(vb3x,vb3y,vb3z)

	ax = vb1y*vb2zm - vb1z*vb2ym
	ay = vb1z*vb2xm - vb1x*vb2zm
	az = vb1x*vb2ym - vb1y*vb2xm
	bx = vb3y*vb2zm - vb3z*vb2ym
	by = vb3z*vb2xm - vb3x*vb2zm
	bz = vb3x*vb2ym - vb3y*vb2xm

	rasq = ax*ax + ay*ay + az*az
	rbsq = bx*bx + by*by + bz*bz
	rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm
	rg = sqrt(rgsq)

	rginv = ra2inv = rb2inv = 0.0
	if (rg > 0): rginv = 1.0/rg
	if (rasq > 0): ra2inv = 1.0/rasq
	if (rbsq > 0): rb2inv = 1.0/rbsq
	rabinv = sqrt(ra2inv*rb2inv)

	c = (ax*bx + ay*by + az*bz)*rabinv
	s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z)

	if (c > 1.0): c = 1.0
	if (c < -1.0): c = -1.0
	return 180.0*math.atan2(s,c)/math.pi
	'''

def rotate_about_dihedral(atoms, dihedral, angle):
	angle = -angle
	atoms_to_rotate = {}
	starting_atom = dihedral.atoms[2]
	def recurse(atom):
		if atom==starting_atom or atom==dihedral.atoms[1]: return
		if atom not in atoms_to_rotate:
			#print atom.index
			atoms_to_rotate[atom] = True
			[recurse(a) for a in atom.bonded]
	[recurse(a) for a in starting_atom.bonded]
	#print [a.index for a in dihedral.atoms]
	#print [a.index for a in atoms_to_rotate], len(atoms_to_rotate)
	#raise Exception()
	origin = dihedral.atoms[2].x, dihedral.atoms[2].y, dihedral.atoms[2].z
	axis = dihedral.atoms[1].x-dihedral.atoms[2].x, dihedral.atoms[1].y-dihedral.atoms[2].y, dihedral.atoms[1].z-dihedral.atoms[2].z
	length = sum([x**2 for x in axis])**0.5
	axis = [x/length for x in axis]
	angle *= math.pi/180
	for atom in atoms_to_rotate:
		v = atom.x, atom.y, atom.z
		v = numpy.subtract(v, origin)
		v2 = numpy.dot(v, math.cos(angle)) + numpy.dot( numpy.cross(axis, v), math.sin(angle) ) + numpy.dot(axis, numpy.dot(axis,v)*(1-math.cos(angle)))
		atom.x, atom.y, atom.z = numpy.add(v2, origin)


def lj_energy(atoms, bonds, angles, dihedrals, bonded=None, angled=None, dihedraled=None, lj_eps=None):
	if not bonded or not angled or not dihedraled:
		bonded = [[] for a in atoms]
		for b in bonds:
			bonded[b.atoms[0].index-1].append( b.atoms[1] )
			bonded[b.atoms[1].index-1].append( b.atoms[0] )
		angled = [[] for a in atoms]
		for a in angles:
			angled[a.atoms[0].index-1].append( a.atoms[2] )
			angled[a.atoms[2].index-1].append( a.atoms[0] )
		dihedraled = [[] for a in atoms]
		for d in dihedrals:
			dihedraled[d.atoms[0].index-1].append( d.atoms[3] )
			dihedraled[d.atoms[3].index-1].append( d.atoms[0] )
	
	if not lj_eps:
		lj_eps = {}
		for i,a in enumerate(atoms):
			for b in atoms[i+1:]:
				eps = math.sqrt(a.type.vdw_e*b.type.vdw_e)
				lj_eps[(a.type,b.type)] = eps
	
	K = 332.063708371
	
	vdw_e = 0.
	#coul_e = 0.
	for i,a in enumerate(atoms):
		for b in atoms[i+1:]:
			d_sq = dist_squared(a,b)
			if d_sq>10.0**2: continue
			bd = a in bonded[ b.index-1 ]
			ad = a in angled[ b.index-1 ]
			dd = a in dihedraled[ b.index-1 ]
			eps = lj_eps[(a.type,b.type)]
			sigma_sq = (a.type.vdw_r*b.type.vdw_r)
			if not bd and not ad:
				vdw_e += (0.5 if dd else 1.0) * 4.*eps*( (sigma_sq/d_sq)**6. - (sigma_sq/d_sq)**3.)
				#coul_e += (0.5 if dd else 1.0) * K * a.charge*b.charge / math.sqrt(d_sq)

	return vdw_e

def opls_energy(coords, atoms, bonds, angles, dihedrals, bonded=None, angled=None, dihedraled=None, no_dihedrals=False, dihedrals_only=False, list_components=False):
	starting_coords = [(a.x, a.y, a.z) for a in atoms]
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = coords[i]
	
	if not bonded or not angled or not dihedraled:
		bonded = [[] for a in atoms]
		for b in bonds:
			bonded[b.atoms[0].index-1].append( b.atoms[1] )
			bonded[b.atoms[1].index-1].append( b.atoms[0] )
		angled = [[] for a in atoms]
		for a in angles:
			angled[a.atoms[0].index-1].append( a.atoms[2] )
			angled[a.atoms[2].index-1].append( a.atoms[0] )
		dihedraled = [[] for a in atoms]
		for d in dihedrals:
			dihedraled[d.atoms[0].index-1].append( d.atoms[3] )
			dihedraled[d.atoms[3].index-1].append( d.atoms[0] )
	
	K = 332.063708371
	
	vdw_e = 0.
	coul_e = 0.
	bond_e = 0.
	angle_e = 0.
	if not dihedrals_only:
		for i,a in enumerate(atoms):
			for b in atoms[i+1:]:
				d_sq = dist_squared(a,b)
				if d_sq>10.0**2: continue
				bd = a in bonded[ b.index-1 ]
				ad = a in angled[ b.index-1 ]
				dd = a in dihedraled[ b.index-1 ]
				eps = math.sqrt(a.type.vdw_e*b.type.vdw_e)
				sigma_sq = (a.type.vdw_r*b.type.vdw_r)
				if not bd and not ad:
					vdw_e += (0.5 if dd else 1.0) * 4.*eps*( (sigma_sq/d_sq)**6. - (sigma_sq/d_sq)**3.)
					coul_e += (0.5 if dd else 1.0) * K * a.charge*b.charge / math.sqrt(d_sq)

		for i,b in enumerate(bonds):
			bond_e += b.type.e*(b.type.r-dist(b.atoms[0], b.atoms[1]))**2
	
		for i,a in enumerate(angles):
			theta = angle_size(a.atoms[0], a.atoms[1], a.atoms[2])
			dif = min(( (a.type.angle-theta)**2, (a.type.angle-theta+360)**2, (a.type.angle-theta-360)**2))
			angle_e += a.type.e*(math.pi/180)**2 * dif
	
	dihedral_e = 0.
	if not no_dihedrals:
		for i,d in enumerate(dihedrals):
			try:
				k1,k2,k3 = d.type.e
			except:
				k1,k2,k3,k4 = d.type.e
			psi, cos1, cos2, cos3, cos4 = dihedral_angle(*d.atoms)
			dihedral_e += 0.5 * ( k1*(1+cos1) + k2*(1-cos2) + k3*(1+cos3) )
	
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = starting_coords[i]
	
	#print vdw_e, coul_e, bond_e, angle_e, dihedral_e
	if list_components:
		return vdw_e, coul_e, bond_e, angle_e, dihedral_e
	else:
		return vdw_e + coul_e + bond_e + angle_e + dihedral_e

import subprocess,re
def energy_compare(coords, atoms, bonds, angles, dihedrals, params):
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = coords[i]
	
	bonded = [[] for a in atoms]
	for b in bonds:
		bonded[b.atoms[0].index-1].append( b.atoms[1] )
		bonded[b.atoms[1].index-1].append( b.atoms[0] )
	angled = [[] for a in atoms]
	for a in angles:
		angled[a.atoms[0].index-1].append( a.atoms[2] )
		angled[dihedral_anglea.atoms[2].index-1].append( a.atoms[0] )
	dihedraled = [[] for a in atoms]
	for d in dihedrals:
		dihedraled[d.atoms[0].index-1].append( d.atoms[3] )
		dihedraled[d.atoms[3].index-1].append( d.atoms[0] )
	
	K = 332.063708371
	
	vdw_e = 0.
	coul_e = 0.
	for i,a in enumerate(atoms):
		for b in atoms[i+1:]:
			d_sq = dist_squared(a,b)
			if d_sq>10.0**2: continue
			bd = a in bonded[ b.index-1 ]
			ad = a in angled[ b.index-1 ]
			dd = a in dihedraled[ b.index-1 ]
			eps = (a.type.vdw_e*b.type.vdw_e)**0.5
			sigma_sq = (a.type.vdw_r*b.type.vdw_r)
			if not bd and not ad:
				vdw_e += (0.5 if dd else 1.0) * 4.*eps*( (sigma_sq/d_sq)**6. - (sigma_sq/d_sq)**3.)
				if d_sq < 8.0**2:
					coul_e += (0.5 if dd else 1.0) * K * a.charge*b.charge / d_sq**0.5
	
	bond_e = 0.
	for i,b in enumerate(bonds):
		bond_e += params.bond_e[i]*(b.d-dist(b.atoms[0], b.atoms[1]))**2
	
	angle_e = 0.
	for i,a in enumerate(angles):
		theta = angle_size(a.atoms[0], a.atoms[1], a.atoms[2])
		dif = min(( (a.theta-theta)**2, (a.theta-theta+360)**2, (a.theta-theta-360)**2))
		angle_e += params.angle_e[i]*(math.pi/180)**2 * dif
	
	cos = math.cos
	dihedral_e = 0.
	#dihedral_e2 = 0.
	
	#for i,theta in enumerate(lammps.get_dihedral_values(atoms, bonds, angles, dihedrals, starting_params)):
	#	dihedrals[i].theta = theta
	
	for i,d in enumerate(dihedrals):
		k1,k2,k3,k4 = params.dihedral_e[i]
		psi = dihedral_angle(*d.atoms, k1=k1,k2=k2,k3=k3,k4=k4)
		#print psi, d.theta
		psi /= 180/math.pi
		dihedral_e += 0.5 * (  k1*(1+cos(psi)) + k2*(1-cos(2*psi)) + k3*(1+cos(3*psi)) + k4*(1-cos(4*psi))  )
		#dihedral_e2 += e
	
	box_size = (100,100,100)
	os.chdir('lammps')
	f = open('test.data', 'w')
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
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type a. b. c.
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, i+1]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()

	f = open('test.in', 'w')
	f.write('''units	real
atom_style	full #bonds, angles, dihedrals, impropers, charges

#pair_style	lj/cut/coul/long 10.0 8.0
pair_style lj/cut/coul/cut 10.0 8.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
#kspace_style pppm 1.0e-6
special_bonds lj/coul 0.0 0.0 0.5

read_data	test.data

thermo_style custom etotal

run 0
''')
	f.close()
	result = subprocess.Popen('~/Documents/grok/lammps/trunk2/lmp_mycomp < test.in', shell=True, stdout=subprocess.PIPE).communicate()[0]
	print result
	e = float( re.search('TotEng\s+(\S+)', result).group(1) )
	os.chdir('..')
	raise Exception(vdw_e + coul_e + bond_e + angle_e + dihedral_e, e)



def residual_E_coeffs(coords, atoms, bonds, angles, dihedrals, dihedral_types_list=None, coeffs_per_dihedral=3):
	starting_coords = [(a.x, a.y, a.z) for a in atoms]
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = coords[i]
	
	bonded = [[] for a in atoms]
	for b in bonds:
		bonded[b.atoms[0].index-1].append( b.atoms[1] )
		bonded[b.atoms[1].index-1].append( b.atoms[0] )
	angled = [[] for a in atoms]
	for a in angles:
		angled[a.atoms[0].index-1].append( a.atoms[2] )
		angled[a.atoms[2].index-1].append( a.atoms[0] )
	dihedraled = [[] for a in atoms]
	for d in dihedrals:
		dihedraled[d.atoms[0].index-1].append( d.atoms[3] )
		dihedraled[d.atoms[3].index-1].append( d.atoms[0] )
	
	K = 332.063708371
	
	residual_E_coeffs = []
	
	'''for i,b in enumerate(bonds):
		residual_E_coeffs.append( (b.d-dist(b.atoms[0], b.atoms[1]))**2 )

	for i,a in enumerate(angles):
		theta = angle_size(a.atoms[0], a.atoms[1], a.atoms[2])
		dif = min(( (a.theta-theta)**2, (a.theta-theta+360)**2, (a.theta-theta-360)**2))
		residual_E_coeffs.append( (math.pi/180)**2 * dif )
	'''
	if not dihedral_types_list: #each dihedral has its own type
		for i,d in enumerate(dihedrals):
			psi, cos1, cos2, cos3, cos4 = dihedral_angle(*d.atoms)
			residual_E_coeffs += [0.5*(1+cos1),  0.5*(1-cos2), 0.5*(1+cos3)] + ([] if coeffs_per_dihedral==3 else [0.5*(1-cos4)])
	else:
		for d_type in dihedral_types_list:
			terms = [0.0]*coeffs_per_dihedral
			for d in dihedrals:
				this_type = tuple([a.type.index2 for a in d.atoms])
				if this_type==d_type or this_type==tuple(reversed(d_type)):
					psi, cos1, cos2, cos3, cos4 = dihedral_angle(*d.atoms)
					terms[0] += 0.5*(1+cos1)
					terms[1] += 0.5*(1-cos2)
					terms[2] += 0.5*(1+cos3)
					if coeffs_per_dihedral>3:
						terms[3] += 0.5*(1-cos4)
			residual_E_coeffs += terms
			
	for i,a in enumerate(atoms):
		a.x, a.y, a.z = starting_coords[i]
	return residual_E_coeffs

