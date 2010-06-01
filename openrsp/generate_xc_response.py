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

def flatten(l_in):
#   example:
#   l_in = [[[1], 2], 3, [[4]]]
#   returns [1, 2, 3, 4]
    l = []
    for x in l_in:
        if hasattr(x, "__iter__") and not isinstance(x, basestring):
           l.extend(flatten(x))
        else:
           l.append(x)
    return l

def multiply_list(l):
#   example:
#   l = [[1, 2], [3, 4, 5]]
#   returns [[1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5]]
    l_old = [[]]
    for i in range(len(l)):
        l_new = []
        for j in range(len(l_old)):
            for k in range(len(l[i])):
                l_new.append([l_old[j], l[i][k]])
        l_old = l_new[:]
    l_out = []
    for l in l_old:
        l_out.append(flatten(l))
    return l_out

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
#   example: 
#   index_s = '234' 
#   returns ['0-234', '2-34', '3-24', '4-23']
#   note that the "other" combinations are missing on purpose:
#                 '234-0', '34-2', ...
#   this is compensated by the factor 2 in front of all z variables
    l = []
    for i in range(len(index_s) + 1):
        for c in xuniqueCombinations(index_s, i):
            left  = ''.join(c)
            if left == '':
               left = '0'
            right = remove_from_s(index_s, ''.join(c))
            if right == '':
               right = '0'
            l.append(left + '-' + right)
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

def collect_terms(term_d, term_d_swap, term_variable_l, term_index_l):

    t = count_variables(variable_l, term_variable_l)

#   avoid double counting:
    accept_term = 1
    for i in range(len(term_variable_l)):
        for j in range(i + 1, len(term_variable_l)):
            term_variable_l_swap    = term_variable_l[:]
            term_variable_l_swap[i] = term_variable_l[j]
            term_variable_l_swap[j] = term_variable_l[i]
            term_index_l_swap       = term_index_l[:]
            term_index_l_swap[i]    = term_index_l[j]
            term_index_l_swap[j]    = term_index_l[i]
            if t + tuple(term_variable_l_swap + term_index_l_swap) in term_d_swap.keys():
               accept_term = 0

    if accept_term:
       term_d_swap[t + tuple(term_variable_l + term_index_l)] = 1
       l = []
       for i in range(len(term_variable_l)):
           l.append([term_variable_l[i], term_index_l[i]])
       if t in term_d.keys():
          term_d[t].append(l)
       else:
          term_d[t] = [l]

    return term_d, term_d_swap

def print_f90_code(term_d, order, index_to_matrix):

    s = ''
    for term in term_d:
        if term_d[term] != []:

           nr_n = term[0]
           nr_z = term[1]
           nr_g = term[2]
          
           if nr_g == 1:
             s += '\n    tv = '
           else:
             s += '\n    tu = '
          
           for i in range(len(term_d[term])):
               for j in range(len(term_d[term][i])):
                   if i > 0 and j == 0:
                      s += ' &\n       + '
                   if j > 0:
                      s += '*'
                   if term_d[term][i][j][0] == 'z':
                      left  = string.split(term_d[term][i][j][1], '-')[0]
                      right = string.split(term_d[term][i][j][1], '-')[1]
                      s += 'z(%s, %s)' % (index_to_matrix[left], index_to_matrix[right])
                   elif term_d[term][i][j][0] == 'g':
                      s += 'g(:, %s)' % index_to_matrix[term_d[term][i][j][1]]
                   else:
                      s += 'n(%s)' % index_to_matrix[term_d[term][i][j][1]]
          
           if nr_g == 1:
              s += '\n    v  = v + derv(xc_index(XC_RSZXY, (/%i, 0, %i, 0, 0/)))*tv(:)'      % (nr_n,     nr_z + 1)
           else:
              s += '\n    v  = v + derv(xc_index(XC_RSZXY, (/%i, 0, %i, 0, 0/)))*tu*g(:, 0)' % (nr_n,     nr_z + 1)
              s += '\n    u  = u + derv(xc_index(XC_RSZXY, (/%i, 0, %i, 0, 0/)))*tu'         % (nr_n + 1, nr_z)

    s = string.replace(s, 'g', 'gn')


    print '''  subroutine xc_response_%i(u, v, derv, n, gn, z, use_alda) 
    
!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_alda
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------
%s
    
  end subroutine
''' % (order + 1, s)

#-------------------------------------------------------------------------------

def generate_contributions(order):

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
    term_d_swap = {}

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
                          if (order < 2) or (x != '0' + '-' + str(matrix_to_index[2**order - 1])):
                             temp_l.append(x)
                       term_index_l.append(temp_l)
                    else:
                       if (order < 2) or (len(variable_combination_l) > 1):
#                         except linear response highest order terms do not contribute
                          term_variable_l.append(variable_combination_l[i])
                          term_index_l.append([selection])

                if len(term_variable_l) > 0:
                   for l in multiply_list(term_index_l):
                       term_d, term_d_swap = collect_terms(term_d, term_d_swap, term_variable_l, l)

    print_f90_code(term_d, order, index_to_matrix)

#-------------------------------------------------------------------------------

# generate contributions up to max_order
for i in range(1, max_order + 1):
    generate_contributions(i)
