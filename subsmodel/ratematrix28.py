# this function takes rates c (gene conversion), g (gain) and d (deletion)
# and returns a matrix showing transitions from each of the 28 states to each
# of the others.  States are in the canonical order:
# null, A, B, AA, AB, BB, AAA, AAB, ABB, BBB, AAAA, AAAB, AABB, ABBB, BBBB
# AAAAA, AAAAB, AAABB, AABBB, ABBBB, BBBBB, AAAAAA, AAAAAB, AAAABB, AAABBB,
# AABBBB, ABBBBB, BBBBBB

numstates = 28
c = 0.9
d = 0.2
g = 0.3

def rates(c, d, g):
  r = [[0.0 for x in xrange(numstates)] for x in xrange(numstates)]
  # convenience names for rows and columns

#Original parameterization
#  null = 0
#  a = 1
#  b = 2
#  aa = 3 
#  ab = 4
#  bb = 5
#  aaa = 6
#  aab = 7
#  abb = 8
#  bbb = 9
#  aaaa = 10
#  aaab = 11
#  aabb = 12
#  abbb = 13
#  bbbb = 14
#  aaaaa = 15
#  aaaab = 16
#  aaabb = 17
#  aabbb = 18
#  abbbb = 19
#  bbbbb = 20
#  aaaaaa = 21
#  aaaaab = 22
#  aaaabb = 23
#  aaabbb = 24
#  aabbbb = 25
#  abbbbb = 26
#  bbbbbb = 27
  
  null=0
  b=1
  bb=2
  bbb=3
  bbbb=4
  bbbbb=5
  bbbbbb=6
  a=7
  ab=8
  abb=9
  abbb=10
  abbbb=11
  abbbbb=12
  aa=13
  aab=14
  aabb=15
  aabbb=16
  aabbbb=17
  aaa=18
  aaab=19
  aaabb=20
  aaabbb=21
  aaaa=22
  aaaab=23
  aaaabb=24
  aaaaa=25
  aaaaab=26
  aaaaaa=27

  # fill in non-zero entries only, omitting diagonals
  # syntax is r[to][from], watch out for this!

  # null goes to nothing

  # naked A goes to AA, null
  r[null][a] = r[null][b] = d
  r[aa][a] = r[bb][b] = g

  # AA goes to A, AAA
  r[a][aa] = r[b][bb] = 2.0*d
  r[aaa][aa] = r[bbb][bb] = 2.0*g

  # AAA goes to AA, AAAA
  r[aa][aaa] = r[bb][bbb] = 3.0*d
  r[aaaa][aaa] = r[bbbb][bbb] = 3.0*g

  # AAAA goes to AAA, AAAAA
  r[aaa][aaaa] = r[bbb][bbbb] = 4.0*d
  r[aaaaa][aaaa] = r[bbbbb][bbbb] = 4.0*g   # new

  # AAAAA goes to AAAA, AAAAAA
  r[aaaa][aaaaa] = r[bbbb][bbbbb] = 5.0*d   # new
  r[aaaaaa][aaaaa] = r[bbbbbb][bbbbb] = 5.0*g  # new

  # AAAAAA goes to AAAAA
  r[aaaaa][aaaaaa] = r[bbbbb][bbbbbb] = 6.0*d  # new

  # AB goes to AA, BB, A, B, AAB, ABB
  r[aa][ab] = r[bb][ab] = c
  r[a][ab] = r[b][ab] = d
  r[aab][ab] = r[abb][ab] = g

  # AAB goes to AA, AB, AAA, ABB, AAAB, AABB
  r[aa][aab] = r[bb][abb] = d
  r[ab][aab] = r[ab][abb] = 2.0*d
  r[aaa][aab] = r[bbb][abb] = c
  r[abb][aab] = r[aab][abb] = 2.0 * 0.5 * c
  r[aaab][aab] = r[abbb][abb] = 2.0 * g
  r[aabb][aab] = r[aabb][abb] = g
 
  # AAAB goes to AAA, AAB, AAAA, AABB, AAAAB, AAABB
  r[aaa][aaab] = r[bbb][abbb] = d
  r[aab][aaab] = r[abb][abbb] = 3.0*d
  r[aaaa][aaab] = r[bbbb][abbb] = c
  r[aabb][aaab] = r[aabb][abbb] = 3.0 * 1.0/3.0 * c
  r[aaaab][aaab] = r[abbbb][abbb] = 3.0 * g    # new
  r[aaabb][aaab] = r[aabbb][abbb] = g          # new

  # AABB goes to AAB, ABB, AAAB, ABBB, AAABB, AABBB
  r[aab][aabb] = r[abb][aabb] = 2.0 * d
  r[aaab][aabb] = r[abbb][aabb] = 2.0 * 2.0/3.0 * c
  r[aaabb][aabb] = r[aabbb][aabb] = 2.0 * g    # new

  # AAAAB goes to AAAA, AAAB, AAAAA, AAABB, AAAAAB, AAAABB  # all new
  r[aaaa][aaaab] = r[bbbb][abbbb] = d         
  r[aaab][aaaab] = r[abbb][abbbb] = 4.0 * d
  r[aaaaa][aaaab] = r[bbbbb][abbbb] = c
  r[aaabb][aaaab] = r[aabbb][abbbb] = 4.0 * 1.0/4.0 * c
  r[aaaaab][aaaab] = r[abbbbb][abbbb] = 4.0 * g
  r[aaaabb][aaaab] = r[aabbbb][abbbb] = g

  # AAABB goes to AAAB, AABB, AAAAB, AABBB, AAAABB, AAABBB  # all new
  r[aaab][aaabb] = r[abbb][aabbb] = 2.0 * d
  r[aabb][aaabb] = r[aabb][aabbb] = 3.0 * d
  r[aaaab][aaabb] = r[abbbb][aabbb] = 2.0 * 3.0/4.0 * c
  r[aabbb][aaabb] = r[aaabb][aabbb] = 3.0 * 2.0/4.0 * c
  r[aaaabb][aaabb] = r[aabbbb][aabbb] = 3.0 * g
  r[aaabbb][aaabb] = r[aaabbb][aabbb] = 2.0 * g

  # AAAAAB goes to AAAAA, AAAAB, AAAAAA, AAAABB  # all new
  r[aaaaa][aaaaab] = r[bbbbb][abbbbb] = d
  r[aaaab][aaaaab] = r[abbbb][abbbbb] = 5.0 * d
  r[aaaaaa][aaaaab] = r[bbbbbb][abbbbb] = c
  r[aaaabb][aaaaab] = r[aabbbb][abbbbb] = 5.0 * 1.0/5.0 * c

  # AAAABB goes to AAAAB, AAABB, AAAAAB, AAABBB  # all new
  r[aaaab][aaaabb] = r[abbbb][aabbbb] = 2.0 * d
  r[aaabb][aaaabb] = r[aabbb][aabbbb] = 4.0 * d
  r[aaaaab][aaaabb] = r[abbbbb][aabbbb] = 2.0 * 4.0/5.0 * c
  r[aaabbb][aaaabb] = r[aaabbb][aabbbb] = 4.0 * 2.0/5.0 * c

  # AAABBB goes to AAABB, AABBB, AAAABB, AABBBB  # all new
  r[aaabb][aaabbb] = r[aabbb][aaabbb] = 3.0 * d
  r[aaaabb][aaabbb] = r[aabbbb][aaabbb] = 3.0 * 3.0/5.0 * c

  # fill in diagonals as -(sum of column)
  for column in xrange(numstates):
    columntot = 0.0
    for row in xrange(numstates):
      if row != column:
        columntot += r[row][column]
    r[column][column] = -columntot
  
  return r


# main program for testing
import sys
import numpy.linalg
import math

#c = float(sys.argv[1])
#d = float(sys.argv[2])
#g = float(sys.argv[3])

m = rates(c,d,g)

#for row in m:
#  for entry in row:
#    print "{0:.4f}".format(entry),
#  print

A = numpy.matrix(m)
(w,v) = numpy.linalg.eig(A)
B = numpy.diag(w)
T = numpy.matrix(v)
Tinv = numpy.linalg.inv(T)

#for i in xrange(len(m)):
#  for j in xrange(len(m)):
#    print "{0:0.4f}".format(m[i][j]),
#  print

import pickle

picklefile = open("rates.pkl","w")
pickle.dump(B, picklefile)
pickle.dump(T, picklefile)
pickle.dump(Tinv, picklefile)
picklefile.close()

numpy.savetxt('B.txt',B,delimiter=',')
numpy.savetxt('T.txt',T,delimiter=',')
numpy.savetxt('Tinv.txt',Tinv,delimiter=',')
numpy.savetxt('A.txt',A,delimiter=',')
