# This program uses the 28-state model to evaluate the likelihood
# of a tiny tree at various branch lengths, demonstrating how
# the evaluations work.  It relies on a rate matrix made by
# program ratematrix28.py, and uses the eigenvalue/eigenvector 
# approach to compute the likelihoods.

epsilon = 0.00000000001 # 1e-11

null = 0
a = 7
b = 1
aa = 13
ab = 8
bb = 2
aaa = 18
aab = 14
abb = 9
bbb = 3
aaaa = 22
aaab = 19
aabb = 15
abbb = 10
bbbb = 4
aaaaa = 25
aaaab = 23
aaabb = 20
aabbb = 16
abbbb = 11
bbbbb = 5
aaaaaa = 27
aaaaab = 26
aaaabb = 24
aaabbb = 21
aabbbb = 17
abbbbb = 12
bbbbbb = 6


def notzero(x):
  return x < epsilon

# unpickle the components of the rate matrices
import pickle
picklefile = open("rates.pkl","r")
# b matrix
bmatrix = pickle.load(picklefile)
# t matrix
tmatrix = pickle.load(picklefile)
# t inverse
tinverse = pickle.load(picklefile)
picklefile.close()

# read in the data
infile = open("data","r")
seq1 = infile.readline()
seq1 = seq1.rstrip()
seq1 = seq1.split()
seq2 = infile.readline()
seq2 = seq2.rstrip()
seq2 = seq2.split()
numsites = len(seq1)
numstates = 28
assert numsites == len(seq2)
infile.close()

# iterate over t1 and t2
testvals = [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3,0.325,0.35,0.375,0.4]
scores = []
for t1 in testvals:
  myscores = []
  for t2 in testvals:
    # set up three dlcells
    dl1 = [[0.0 for x in xrange(numstates)] for x in xrange(numsites)]
    dl2 = [[0.0 for x in xrange(numstates)] for x in xrange(numsites)]
    dl3 = [[0.0 for x in xrange(numstates)] for x in xrange(numsites)]
    
    for n in xrange(numsites):
      site = int(seq1[n])
      dl1[n][site] = 1.0
      site = int(seq2[n])
      dl2[n][site] = 1.0
    
    # WATCH OUT probably need logs here!
    # compute probabilities down branches
    # exponentiate for t1
    
    import numpy
    import copy
    import math
    
    b1 = copy.deepcopy(bmatrix)
    for i in xrange(len(b1)):
      b1[i][i] = math.exp(b1[i][i] * t1) 
    
    p1 = tmatrix * b1 * tinverse
    p1 = p1.tolist()
    
    b2 = copy.deepcopy(bmatrix)
    for i in xrange(len(b2)):
      b2[i][i] = math.exp(b2[i][i] * t2) 
    
    p2 = tmatrix * b2 * tinverse
    p2 = p2.tolist()
    
    # trying without logs for now, we'll see....
    #for i in xrange(len(newprobs)):
    #  for j in xrange(len(newprobs)):
    #    if notzero(newprobs[i][j]):  # this line zeroes any element < epsilon
    #      newprobs[i][j] = math.log(newprobs[i][j])
    #    else:
    #      newprobs[i][j] = 0.0
      
    # compute dl3 values based on dl1, dl2, and newprobs 
    for site in xrange(numsites):
      for top1 in xrange(numstates):
        for top2 in xrange(numstates):
          for bottom in xrange(numstates):
            dl3[site][bottom] += dl1[site][top1] * p1[top1][bottom] * dl2[site][top2] * p2[top2][bottom]
            
    # compute data likelihood assuming ancestor was state 4 with 100%
    # probability; here we go to logs

    lnlike = 0.0
    for site in xrange(numsites):
      lnlike += math.log(dl3[site][ab])
    myscores.append(lnlike)
  scores.append(myscores[:])
  
#import matplotlib.pyplot as plt

bestval = scores[0][0]
best1 = testvals[0]
best2 = testvals[0]
for t1 in xrange(len(testvals)):
  for t2 in xrange(len(testvals)):
    if scores[t1][t2] > bestval:
      bestval = scores[t1][t2]
      best1 = testvals[t1]
      best2 = testvals[t2]

print bestval, best1, best2
#plt.imshow(scores)
#plt.show()
