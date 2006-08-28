#
# Convert Dyall basis set files
#
# Andre S.P. Gomes <a.gomes@few.vu.nl>, Ch. R. Jacob <jacob@few.vu.nl>, Dirac meeting 2006
#

import re

class basis_set :

	def __init__ (self) :
		self.exponents   = {}
		self.coefficients = {}

	def add_exponent (self, sym, exp) :

		if sym not in self.exponents :
			self.exponents[sym] = []
			self.coefficients[sym] = []
		self.exponents[sym].append(exp)
		self.coefficients[sym].append([]) 

	def add_coefficients (self, sym, exp_index, coefs) :
		"""
		coefs: list of coefficents 
		"""

		self.coefficients[sym][exp_index] += coefs

	def get_number_of_exponents (self, sym) :
		return len(self.exponents[sym])

	def get_number_of_orbitals (self, sym) :
		return len(self.coefficients[sym][0])

	def get_symmetries (self) :
		syms = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k']

		for i, s in enumerate(syms) :
			if s in self.exponents :
				max_i = i

		return syms[:max_i+1]


	def dirac_format (self) :
		string = ""

		for sym in self.get_symmetries() :
			if sym in self.exponents :
				string += "$ "+sym+" functions \n"
			
				string += "%5i%5i%5i\n" % (self.get_number_of_exponents(sym),
							   self.get_number_of_orbitals(sym),
							   0)                                # always zero, nobody knows why

				for exp_index, exp, coefs in zip(range(self.get_number_of_exponents(sym)),
								 self.exponents[sym],
								 self.coefficients[sym]) :

					string += "%15.8e" % exp

					i = 0
					for i, c in enumerate(coefs) :
						if ((i+2) % 7 == 1) :
							string+='   '
						string += " %11.4e" % c
						if ((i+2) % 7 == 0) :
							string += '\n'
					if not ((i+2) % 7 == 0) :
						string += "\n"
			else:
				string += "$ "+sym+" functions \n"
				string += "%5i%5i%5i\n" % (0, 0, 0)
		return string


class dyall_set :

	def __init__ (self, element, filename) :
		self.element  = element
		self.filename = filename
		
		self.pol_sets = {}

	def scf_reader (self, type, functions) :
		f = file (self.filename, 'r')
		lines = f.readlines ()

		start1 = re.compile (type+" basis sets?")
		start2 = re.compile (functions+" functions")
		atom   = re.compile (r"\*\* "+self.element+" \*\*")
		functions = re.compile (r".* functions\s*$")

		symline = re.compile (r"^\s+\d([spdfghik])[+-]?")
		index   = re.compile (r"^\s*([-+]?\d+)\s") 
		number  = re.compile (r'(?P<num>[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?)')

		for i, l in enumerate(lines) :
			if start1.match (l) :
				break
		lines = lines[i+1:]

		for i, l in enumerate(lines) :
			if atom.match (l) :
				break
		lines = lines[i+1:]

		for i, l in enumerate(lines) :
			if start2.search (l) :
				break
		lines = lines[i+1:]

		bs = basis_set ()

		sym = ''
		for l in lines :

			if functions.match(l) :
				break

			symline_match = symline.match(l)
			index_match   = index.match(l)
			
			if symline_match :
				new_sym = symline_match.group(1)
				if new_sym == sym:
					coeffs_only = True
				else:
					coeffs_only = False
				sym = new_sym

			elif index_match :

				i = int(index_match.group(1))-1

				l = l[index_match.end(1)+1:]

				line_nums = [ float(m.group("num")) for m in number.finditer(l) ]
				if not coeffs_only :
					bs.add_exponent(sym, line_nums[0])
				bs.add_coefficients(sym, i, line_nums[1:])

		f.close
		return bs

	def pol_reader (self, type) :

		start      = re.compile (r'^Exponents of '+type.replace('(','\(').replace(')','\)')+
					  ' functions?\w*$', re.IGNORECASE)
		symmetries = re.compile (r'\s+([spdfghik])')
		atom       = re.compile (self.element)
		number     = re.compile (r'(?P<num>[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?)')

		f = file(self.filename, 'r')
		lines = f.readlines()
		
		for i, l in enumerate(lines) :
			if start.match (l) :
				break

		syms  = []
		exps = []
		for l in lines[i+1:] :
			if symmetries.match(l) :
				for m in symmetries.finditer(l) :
					syms.append (m.group(1))
			if atom.match(l) :
				for m in number.finditer(l) :
					exps.append (float(m.group('num')))
				break
		
		f.close

		bs = basis_set ()
		for s, e in zip(syms, exps) :
			bs.add_exponent (s, e)

		return bs

	def get_pol_set_types (self) :
		pol_set_types = []
		
		pol_header = re.compile(r'^Exponents of (.+) functions?\w*$', re.IGNORECASE | re.MULTILINE)

		f = file(self.filename, 'r')
		for m in pol_header.finditer(f.read()) :
			pol_set_types.append(m.group(1)) 
		f.close()
		
		return pol_set_types

	def read (self) :
		self.scf_set_DC_large       = self.scf_reader ('Dirac-Hartree-Fock', 'Large component')
		self.scf_set_DC_small       = self.scf_reader ('Dirac-Hartree-Fock', 'Small component')
		self.scf_set_SF_large       = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Large component')
#		self.scf_set_SF_pseudolarge = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Pseudo-large component')
		self.scf_set_SF_FW          = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Foldy-Wouthuysen')

		for t in self.get_pol_set_types() :
			self.pol_sets[t] = self.pol_reader(t)


