#!/usr/bin/env python

from __future__ import generators
import sys

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

def get_variable_combinations(variable_l, order):
    l = []
    for s in xselections(variable_l, order):
        l.append(s)
    for s in xselections(variable_l, order - 1):
        s.append('g')
        l.append(s)
    return l

#-------------------------------------------------------------------------------

index_s = '12345678'
variable_l = ['n', 'Z']

#-------------------------------------------------------------------------------

# a partition of an integer n is an unordered collection
# of positive integers (not necessarily distinct) which add up to n

partition_ll = []

partition_ll.append([[0]])
partition_ll.append([[1]])
partition_ll.append([[2],
                     [1, 1]])
partition_ll.append([[3],
                     [2, 1],
                     [1, 1, 1]])
partition_ll.append([[4],
                     [3, 1],
                     [2, 2],
                     [2, 1, 1],
                     [1, 1, 1, 1]])
partition_ll.append([[5],
                     [4, 1],
                     [3, 2],
                     [3, 1, 1],
                     [2, 2, 1],
                     [2, 1, 1, 1],
                     [1, 1, 1, 1, 1]])
partition_ll.append([[6],
                     [5, 1],
                     [4, 2],
                     [4, 1, 1],
                     [3, 3],
                     [3, 2, 1],
                     [3, 1, 1, 1],
                     [2, 2, 2],
                     [2, 2, 1, 1],
                     [2, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1]])
partition_ll.append([[7],
                     [6, 1],
                     [5, 2],
                     [5, 1, 1],
                     [4, 3],
                     [4, 2, 1],
                     [4, 1, 1, 1],
                     [3, 3, 1],
                     [3, 2, 2],
                     [3, 2, 1, 1],
                     [3, 1, 1, 1, 1],
                     [2, 2, 2, 1],
                     [2, 2, 1, 1, 1],
                     [2, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1]])
partition_ll.append([[8],
                     [7, 1],
                     [6, 2],
                     [6, 1, 1],
                     [5, 3],
                     [5, 2, 1],
                     [5, 1, 1, 1],
                     [4, 4],
                     [4, 3, 1],
                     [4, 2, 2],
                     [4, 2, 1, 1],
                     [4, 1, 1, 1, 1],
                     [3, 3, 2],
                     [3, 3, 1, 1],
                     [3, 2, 2, 1],
                     [3, 2, 1, 1, 1],
                     [3, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2],
                     [2, 2, 2, 1, 1],
                     [2, 2, 1, 1, 1, 1],
                     [2, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1, 1]])

# check that partitions sum up to n
# it does not check whether some partitions are missing
for i in range(len(partition_ll)):
    for j in range(len(partition_ll[i])):
        if (i != sum(partition_ll[i][j])):
           sys.exit('partition %i is wrong' % i)

#-------------------------------------------------------------------------------

def get_nr_terms(order):
    nr_terms = 0
    for partition_l in partition_ll[order]:
        variable_combination_ll = get_variable_combinations(variable_l, len(partition_l))
        for variable_combination_l in variable_combination_ll:
            for partition in get_partitions(index_s[:order], partition_l):
                term_variable_l = []
                term_index_l    = []
                for i in range(len(partition)):
                    selection = partition[i]
                    term_variable_l.append(variable_combination_l[i])
                    if variable_combination_l[i] == 'Z':
                       term_index_l.append(get_partitions_scalar_product(selection))
                    else:
                       term_index_l.append(selection)
    
           #    print term_variable_l, term_index_l
                t = 1
                for i in range(len(term_variable_l)):
                    if term_variable_l[i] == 'Z':
                       t *= len(term_index_l[i])
           #    print t
                nr_terms += t
    print order, nr_terms

get_nr_terms(1)
get_nr_terms(2)
get_nr_terms(3)
get_nr_terms(4)
get_nr_terms(5)
get_nr_terms(6)
get_nr_terms(7)
get_nr_terms(8)
