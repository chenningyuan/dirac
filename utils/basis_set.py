#
# Convert Dyall basis set files
#
# Ch. R. Jacob <jacob@few.vu.nl>, Andre S.P. Gomes <a.gomes@few.vu.nl>, Dirac meeting 2006
#

import re
import copy
import periodic_table

class basis_set :

	def __init__ (self) :
		self.exponents   = {}
		self.coefficients = {}

	def add_exponent (self, sym, exp) :

		if sym not in self.exponents :
			self.exponents[sym] = []
			self.coefficients[sym] = []
		exp_already_there = False
		for e in self.exponents[sym] :
			if (abs((e-exp)/e) < 1e-6) :
				exp_already_there = True
		if exp_already_there == False :
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


	def dirac_format (self,uncontracted=False) :
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

	def __add__ (self, other) :
		# WARNING: This only works for uncontracted sets

		sum = copy.copy (self)
		for s in other.get_symmetries() :
			if s in other.exponents :
				for e in other.exponents[s] :
					sum.add_exponent (s, e)
		return sum		


class dyall_set :

	def __init__ (self, element, filename) :
		self.element  = element
		self.filename = filename
		
		self.pol_sets = {}

	def comment_reader (self, read=False):

		comment_string = ""

		if read :
			f = file (self.filename, 'r')
			lines = f.readlines()
                
			comment_end = re.compile (r"^Dirac-Hartree-Fock basis sets")

			for l in lines:
				if comment_end.match (l) :
					break
				else :
					comment_string = comment_string+"$ "+l

			f.close()

		return comment_string

	def scf_reader (self, type, functions, uncontracted=False) :
		f = file (self.filename, 'r')
		lines = f.readlines ()

		start1 = re.compile (type+" basis sets?")
		start2 = re.compile (functions+" functions")
		atom   = re.compile (r"\*\* "+self.element+"\s{1,2}\*\*")
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
				if uncontracted == False :
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

		syms = []
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

	def read (self, uncontracted=False) :
		self.scf_set_DC_large       = self.scf_reader ('Dirac-Hartree-Fock', 'Large component', uncontracted)
		self.scf_set_DC_small       = self.scf_reader ('Dirac-Hartree-Fock', 'Small component', uncontracted)
		self.scf_set_SF_large       = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Large component', uncontracted)
#		self.scf_set_SF_pseudolarge = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Pseudo-large component', uncontracted)
		self.scf_set_SF_FW          = self.scf_reader ('Spin-free Dirac-Hartree-Fock', 'Foldy-Wouthuysen', uncontracted)

		for t in self.get_pol_set_types() :
			self.pol_sets[t] = self.pol_reader(t)

	def read_header (self, show=False):
		comment = self.comment_reader(show)
		return comment




#
def parse_cmdline():
	from optparse import OptionParser
	
        usage   = "basis_set.py [options]"
        version = "0.01"
	parser  = OptionParser(usage,version=version)

        parser.add_option("--file",dest="raw_basis_file", default=None,
                     help="processes the raw basis set file (default=None)")
        parser.add_option("--output",dest="library_basis_file", default="my-dyall-lib", 
                     help="outputs the dirac basis set library to the given filename (default=my-dyall-lib)")
        parser.add_option("--list_elements",action="store_true", default=False, dest="listelem",
                     help="lists the elements avaiable in the raw file (default=False)")
        parser.add_option("--list_polarization",action="store_true", default=False,dest="listpol", 
                     help="lists the types of extra polarization functions in the raw file (default=False)")
        parser.add_option("--uncontracted",action="store_true", default=False,dest="uncontracted", 
                     help="prints only the exponents, which is equivalent to an uncontracted set (default=False)")

        parser.add_option("--print_scfbas",action="store_true", default=False,dest="print_scf_set", 
                     help="specifies the scf basis set should be printed (default=False)")
        parser.add_option("--print_polarization",action="store_true", default=False,dest="print_pol_set", 
                     help="specifies polarization sets should be printed (default=False)")
        parser.add_option("--polarization_type",default="valence correlating",dest="selected_pol_set", 
                     help="Indicates which polarization set for printing. a few important things here: "+
                          "first, the arguments here are those obtained from the output of "+
                          "--print_polarization, or \"all\" to include all types of functions "+
                          "available; and third, the diffuse sets should also be specified via "+
                          "this flag. (default='valence correlating')")

        parser.add_option("--DC-L",action="store_true", default=True,dest="dc_large", 
                     help="prints large component Dirac-Coulomb basis set (default=True)")
        parser.add_option("--DC-S",action="store_true", default=False,dest="dc_small", 
                     help="prints small component Dirac-Coulomb basis set (default=False)")
        parser.add_option("--SFDC",action="store_true", default=False,dest="sf_large", 
                     help="prints Spin-Free Dirac-Coulomb basis set (default=False)")

        (options, args) = parser.parse_args()
        return options


def get_elements(filename):

	element = re.compile(r'^\*\* ([a-zA-Z]+)\s{1,2}\*\*', re.MULTILINE)
	elementlist = []
        f = file(filename,'r')
        for m in element.finditer(f.read()):
        	elementlist.append(m.group(1)) 
        f.close()

	return set(elementlist)
     

def main():

	options = parse_cmdline()

       	file = options.raw_basis_file
       	elements = get_elements(options.raw_basis_file)

	element_list = ""
	for e in elements:
		element_list = element_list+" "+e

	polarization_list = options.selected_pol_set.split(',') 

        if options.listelem :
		print "The following elements are present in "+options.raw_basis_file
	if options.listpol :
		print "The following polarization sets are present in "+options.raw_basis_file+":"

	show_header = True
	for e in elements:
                if options.listelem:
			print e
		else:
			d = dyall_set(e,file)
                	if options.listpol :
				print e+" = "+str(d.get_pol_set_types())
			else :
#
# first we print the header for this file, containing original comments, elements supported etc 
#
				header = d.read_header(show_header)
				if show_header :
					print "$ "+options.library_basis_file
					print "$"
					print "$  REFERENCE"
					print "$ see dirac website"
					print "$"
					print "$ Elements Supported\n$"+element_list
					print "$"
					print "$ Description"
					print header

				show_header = False
#
# now processing the basis set themselves
#
				d.read(options.uncontracted)
                		print "$ "+e
                		print "a "+str(periodic_table.get_element_number(e))
				set = None
				if options.print_scf_set :
					if options.dc_large:
						set = d.scf_set_DC_large
					if options.dc_small:
						set = d.scf_set_DC_small
					if options.sf_large:
						set = d.scf_set_SF_large
				if options.print_pol_set :
					for s in d.get_pol_set_types() :
						if options.selected_pol_set == "all" or s in polarization_list :
#						print "adding "+s
							if not set :
								set = d.pol_sets[s]
							else :
								set = set + d.pol_sets[s]

				print set.dirac_format()

main()
