import os

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
		
