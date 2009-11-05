#!/usr/bin/env python

from __future__ import generators
import sys
import string

#-------------------------------------------------------------------------------

variable_l = ['n', 'z', 'g']
max_order = 4

is_scalar_product = {}

is_scalar_product['n'] = 0
is_scalar_product['z'] = 1
is_scalar_product['g'] = 0

#-------------------------------------------------------------------------------

class term_c:
      def __init__(self, l, s):
          self.index_partition   = l[:]
          self.available_indices = s

#-------------------------------------------------------------------------------

def xcombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcombinations(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc

def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc

def xselections(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss

def xpermutations(items):
    return xcombinations(items, len(items))

def remove_from_s(s_in, s_remove):
    s_out = ''
    for s in s_in:
        if not s in s_remove:
           s_out += s
    return s_out

def get_perturbation_l(l, order):
    x = []
    for i in range(1, order + 1):
        for c in xuniqueCombinations(l[0:order], i):
            x.append(''.join(c))
    return x

def get_partitions(index_s, partition_l):

    term_l = []
    term = term_c([], index_s)
    term_l.append(term)

    for partition in partition_l:
        if partition == 1:
           for term in term_l:
               term.index_partition.append(term.available_indices[0])
               term.available_indices = term.available_indices[1:]
        else:
           new_term_l = []
           for term in term_l:
               for c in xuniqueCombinations(term.available_indices, partition):
                   new_index_partition = term.index_partition[:]
                   new_index_partition.append(''.join(c))
                   new_available_indices = remove_from_s(term.available_indices, c)
                   new_term = term_c(new_index_partition, new_available_indices)
                   new_term_l.append(new_term)
           term_l = new_term_l[:]

    l = []
    for term in term_l:
        l.append(term.index_partition)
    return l

def get_partitions_scalar_product(index_s):
    l = []
    for i in range(len(index_s) + 1):
        for c in xuniqueCombinations(index_s, i):
            left  = ''.join(c)
            if left == '':
               left = '0'
            right = remove_from_s(index_s, ''.join(c))
            if right == '':
               right = '0'
            l.append([left, right])
    half_len = len(l)/2
    return l[0:half_len]

def count_variables(l_ref, s):
    l = []
    for i in range(len(l_ref)):
        l.append(0)
        for x in s:
            if x == l_ref[i]:
               l[i] += 1
    return tuple(l)

def get_variable_combinations(order):
    l1 = []
    for x in xselections(variable_l, order):
        l1.append(x)
#   remove selections which contain more than one gradient
    l2 = []
    for x in l1:
        if count_variables(['g'], x)[0] < 2:
           l2.append(x)
    return l2

def _get_partition_l(n, l):
    if not n:
       yield []
    for i in range(min(n, l), 0, -1):
        for p in _get_partition_l(n - i, i):
            yield [i] + p

def get_partition_l(n):
    l = []
    for partition in _get_partition_l(n, n):
        l.append(partition)
    return l

#-------------------------------------------------------------------------------

def print_f90_code(term_d, order):

    s = ''
    for term in term_d:

       #s = ''
       #print term
        if 'g' in term_d[term]:
          s += '\n    tv = ' + term_d[term]
        else:
          s += '\n    tu = ' + term_d[term]

        nr_n = term[0]
        nr_z = term[1]

        if 'g' in term_d[term]:
           s += '\n    v  = v + derv(taylor_index(5, (/%i, 0, %i, 0, 0/)))*tv(:)'      % (nr_n,     nr_z + 1)
        else:
           s += '\n    v  = v + derv(taylor_index(5, (/%i, 0, %i, 0, 0/)))*tu*g(:, 0)' % (nr_n,     nr_z + 1)
           s += '\n    u  = u + derv(taylor_index(5, (/%i, 0, %i, 0, 0/)))*tu'         % (nr_n + 1, nr_z)
       #print s

    s = string.replace(s, 'g', 'gn')

    print '''  subroutine xc_response_%i(u, v, derv, n, gn, z) 
    
!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
!   ----------------------------------------------------------------------------
    integer, external      :: taylor_index
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------
%s
    
  end subroutine
''' % (order, s)

def collect_terms(term_d, term_variable_l, term_index_l, index_to_matrix):

    s = ''

    nr_n = 0
    nr_z = 0
    nr_g = 0

    for i in range(len(term_variable_l)):
        if term_variable_l[i] == 'z':
           nr_z += 1
           s += '*'
           if len(term_index_l[i]) > 1:
              s += '('
           j = 0
           for x in term_index_l[i]:
               j += 1
               if j > 1:
                  s += ' + '
               s += '%s(%s, %s)' % (term_variable_l[i], index_to_matrix[x[0]], index_to_matrix[x[1]])
           if len(term_index_l[i]) > 1:
              s += ')'
        elif term_variable_l[i] == 'g':
           nr_g += 1
           s += '*g(:, %s)' % index_to_matrix[term_index_l[i]]
        else:
           nr_n += 1
           s += '*%s(%s)' % (term_variable_l[i], index_to_matrix[term_index_l[i]])
    t = tuple([nr_n, nr_z, nr_g])
    if t in term_d.keys():
       term_d[t] += ' &\n       + ' + s[1:]
    else:
       term_d[t]  = s[1:]
    return term_d

#-------------------------------------------------------------------------------

def generate_code(order):

    index_s = '12345678'
    assert (order <= len(index_s))

#   index_to_matrix dictionary is used to find where a specific matrix is
#   there are 2^order matrices
#   example: order           3
#            index           0 1 2 3 12 13 23 123
#            index_to_matrix 0 1 2 3 4  5  6  7
    i = 0
    index_to_matrix = {}
    matrix_to_index = {}
    for perturbation in get_perturbation_l(index_s, order):
        i += 1
        index_to_matrix[perturbation] = str(i)
        matrix_to_index[i] = perturbation
    index_to_matrix['0'] = '0'
    matrix_to_index[0]   = '0'

    term_d = {}

    for partition_l in get_partition_l(order):
        variable_combination_ll = get_variable_combinations(len(partition_l))
        for variable_combination_l in variable_combination_ll:
            for partition in get_partitions(index_s[:order], partition_l):
                term_variable_l = []
                term_index_l    = []
                for i in range(len(partition)):
                    selection = partition[i]
                    if is_scalar_product[variable_combination_l[i]]:
                       term_variable_l.append(variable_combination_l[i])
                       temp_l = []
                       for x in get_partitions_scalar_product(selection):
#                         except linear response highest order terms do not contribute
                          if (order < 2) or (x != ['0', str(matrix_to_index[2**order - 1])]):
                             temp_l.append(x)
                       term_index_l.append(temp_l)
                    else:
                       if (order < 2) or (len(variable_combination_l) > 1):
#                         except linear response highest order terms do not contribute
                          term_variable_l.append(variable_combination_l[i])
                          term_index_l.append(selection)
#               print term_variable_l, term_index_l
                if len(term_variable_l) > 0:
                   term_d = collect_terms(term_d, term_variable_l, term_index_l, index_to_matrix)

    print_f90_code(term_d, order)

#-------------------------------------------------------------------------------

# generate code up to max_order
for i in range(1, max_order + 1):
    generate_code(i)
