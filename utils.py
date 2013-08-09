import os, math, numpy

class Struct:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
	def __repr__(self):
		return str( dict([ (a,None) if type(self.__dict__[a]) in (list,dict) else (a,self.__dict__[a]) for a in self.__dict__]) )

def unique_filename(directory, name, filetype):
	number = 0
	while True:
		unique_name = name+'_'+str(number)
		number += 1
		if not os.path.exists(directory+unique_name+filetype): return unique_name
		
def dist_squared(a,b):
	return (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2

def rotate_about_dihedral(atoms, dihedral, angle):
	angle = -angle
	atoms_to_rotate = {}
	starting_atom = dihedral.atoms[2]
	def recurse(atom):
		if atom==starting_atom or atom==dihedral.atoms[1]: return
		if atom not in atoms_to_rotate:
			print atom.index
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

def dihedral_angle(a,b,c,d):
	'''dd = dist_squared
	dd12, dd23, dd34, dd24, dd13, dd14 = dd(a,b), dd(b,c), dd(c,d), dd(b,d), dd(a,c), dd(a,d)

	dd12, dd23, dd34, dd24, dd13, dd14 = [i for i in (dd12, dd23, dd34, dd24, dd13, dd14)]

	P = dd12*(dd23+dd34-dd24) + dd23*(-dd23+dd34+dd24) + dd13*(dd23-dd34+dd24) - 2*dd23*dd14
	d12,d23,d34,d24,d13 = [math.sqrt(i) for i in dd12,dd23,dd34,dd24,dd13]
	Q = (d12 + d23 + d13) * ( d12 + d23 - d13) * (d12 - d23 + d13) * (-d12 + d23 + d13 ) * (d23 + d34 + d24) * ( d23 + d34 - d24 ) * (d23 - d34 + d24) * (-d23 + d34 + d24 )
	
	print 180/math.pi * math.acos(P/Q**0.5)
	'''
	u = (0, b.x-a.x, b.y-a.y, b.z-a.z)
	v = (0, c.x-a.x, c.y-a.y, c.z-a.z)
	Ua = (u[2]*v[3] - u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1])
	
	u = (0, c.x-b.x, c.y-b.y, c.z-b.z)
	v = (0, d.x-b.x, d.y-b.y, d.z-b.z)
	Ub = (u[2]*v[3] - u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1])
	
	return 180/math.pi * math.acos(min(( Ua[0]*Ub[0] + Ua[1]*Ub[1] + Ua[2]*Ub[2] )/(sum([x**2 for x in Ua])**0.5 * sum([x**2 for x in Ub])**0.5), 1.))

