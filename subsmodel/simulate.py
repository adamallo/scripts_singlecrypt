# simulate evolution down a tree of two tips according to
# our SCA model.  Uses horserace logic to provide a quasi
# independent check of the matrix math in the evaluator.
# Considers states up to 6 copies

import random
import math
import dendropy
import argparse
import os
import re


###Constants/Variables
######################

failnum = -1.0
bigtime = 10000000.0
maxit=1000

null='@'
b='A'
bb='B'
bbb='C'
bbbb='D'
bbbbb='E'
bbbbbb='F'
a='G'
ab='H'
abb='I'
abbb='J'
abbbb='K'
abbbbb='L'
aa='M'
aab='N'
aabb='O'
aabbb='P'
aabbbb='Q'
aaa='R'
aaab='S'
aaabb='T'
aaabbb='U'
aaaa='V'
aaaab='W'
aaaabb='X'
aaaaa='Y'
aaaaab='Z'
aaaaaa='[' 

charDict={'@':(0,0),'A':(0,1),'B':(0,2),'C':(0,3),'D':(0,4),'E':(0,5),'F':(0,6),'G':(1,0),'H':(1,1),'I':(1,2),'J':(1,3),'K':(1,4),'L':(1,5),'M':(2,0),'N':(2,1),'O':(2,2),'P':(2,3),'Q':(2,4),'R':(3,0),'S':(3,1),'T':(3,2),'U':(3,3),'V':(4,0),'W':(4,1),'X':(4,2),'Y':(5,0),'Z':(5,1),'[':(6,0)}

arrayChar=(('@','A','B','C','D','E','F'),('G','H','I','J','K','L'),('M','N','O','P','Q'),('R','S','T','U'),('V','W','X'),('Y','Z'),('[')) #arrayChar[na][nb]

##FUNCTIONS
##########

##Obtains the scaling factor (beta) to scale the branches given in expected substitutions per site, instead of scaling the whole subsmatrix
def getScalingFactor():
  
  numstates=28

  #Hardcoded indexes!
  iin = 0
  ia = 7
  ib = 1
  iaa = 13
  iab = 8
  ibb = 2
  iaaa = 18
  iaab = 14
  iabb = 9
  ibbb = 3
  iaaaa = 22
  iaaab = 19
  iaabb = 15
  iabbb = 10
  ibbbb = 4
  iaaaaa = 25
  iaaaab = 23
  iaaabb = 20
  iaabbb = 16
  iabbbb = 11
  ibbbbb = 5
  iaaaaaa = 27
  iaaaaab = 26
  iaaaabb = 24
  iaaabbb = 21
  iaabbbb = 17
  iabbbbb = 12
  ibbbbbb = 6

  r=[[0.0 for x in xrange(numstates) ] for y in xrange(numstates)]

  #naked A goes to AA, null
  r[ia][iin] = r[ib][iin] = d
  r[ia][iaa] = r[ib][ibb] = g
  
  # AA goes to A, AAA
  r[iaa][ia] = r[ibb][ib] = 2.0*d
  r[iaa][iaaa] = r[ibb][ibbb] = 2.0*g
  
  # AAA goes to AA, AAAA
  r[iaaa][iaa] = r[ibbb][ibb] = 3.0*d
  r[iaaa][iaaaa] = r[ibbb][ibbbb] = 3.0*g
  
  # AAAA goes to AAA, AAAAA
  r[iaaaa][iaaa] = r[ibbbb][ibbb] = 4.0*d
  r[iaaaa][iaaaaa] = r[ibbbb][ibbbbb] = 4.0*g
  
  # AAAAA goes to AAAA, AAAAAA
  r[iaaaaa][iaaaa] = r[ibbbbb][ibbbb] = 5.0*d
  r[iaaaaa][iaaaaaa] = r[ibbbbb][ibbbbbb] = 5.0*g
  
  # AAAAAA goes to AAAAA
  r[iaaaaaa][iaaaaa] = r[ibbbbbb][ibbbbb] = 6.0*d
  
  # AB goes to AA, BB, A, B, AAB, ABB
  r[iab][iaa] = r[iab][ibb] = c
  r[iab][ia] = r[iab][ib] = d
  r[iab][iaab] = r[iab][iabb] = g
  
  # AAB goes to AA, AB, AAA, ABB, AAAB, AABB
  r[iaab][iaa] = r[iabb][ibb] = d
  r[iaab][iab] = r[iabb][iab] = 2.0*d
  r[iaab][iaaa] = r[iabb][ibbb] = c
  r[iaab][iabb] = r[iabb][iaab] = 2.0 * 0.5 * c
  r[iaab][iaaab] = r[iabb][iabbb] = 2.0 * g
  r[iaab][iaabb] = r[iabb][iaabb] = g
  
  # AAAB goes to AAA, AAB, AAAA, AABB, AAAAB, AAABB
  r[iaaab][iaaa] = r[iabbb][ibbb] = d
  r[iaaab][iaab] = r[iabbb][iabb] = 3.0*d
  r[iaaab][iaaaa] = r[iabbb][ibbbb] = c
  r[iaaab][iaabb] = r[iabbb][iaabb] = 3.0 * 1.0/3.0 * c
  r[iaaab][iaaaab] = r[iabbb][iabbbb] = 3.0 * g
  r[iaaab][iaaabb] = r[iabbb][iaabbb] = g
  
  # AABB goes to AAB, ABB, AAAB, ABBB, AAABB, AABBB
  r[iaabb][iaab] = r[iaabb][iabb] = 2.0 * d
  r[iaabb][iaaab] = r[iaabb][iabbb] = 2.0 * 2.0/3.0 * c
  r[iaabb][iaaabb] = r[iaabb][iaabbb] = 2.0 * g
  
  # AAAAB goes to AAAA, AAAB, AAAAA, AAABB, AAAAAB, AAAABB
  r[iaaaab][iaaaa] = r[iabbbb][ibbbb] = d
  r[iaaaab][iaaab] = r[iabbbb][iabbb] = 4.0 * d
  r[iaaaab][iaaaaa] = r[iabbbb][ibbbbb] = c
  r[iaaaab][iaaabb] = r[iabbbb][iaabbb] = 4.0 * 1.0/4.0 * c
  r[iaaaab][iaaaaab] = r[iabbbb][iabbbbb] = 4.0 * g
  r[iaaaab][iaaaabb] = r[iabbbb][iaabbbb] = g
  
  # AAABB goes to AAAB, AABB, AAAAB, AABBB, AAAABB, AAABBB
  r[iaaabb][iaaab] = r[iaabbb][iabbb] = 2.0 * d
  r[iaaabb][iaabb] = r[iaabbb][iaabb] = 3.0 * d
  r[iaaabb][iaaaab] = r[iaabbb][iabbbb] = 2.0 * 3.0/4.0 * c
  r[iaaabb][iaabbb] = r[iaabbb][iaaabb] = 3.0 * 2.0/4.0 * c
  r[iaaabb][iaaaabb] = r[iaabbb][iaabbbb] = 3.0 * g
  r[iaaabb][iaaabbb] = r[iaabbb][iaaabbb] = 2.0 * g
  
  # AAAAAB goes to AAAAA, AAAAB, AAAAAA, AAAABB
  r[iaaaaab][iaaaaa] = r[iabbbbb][ibbbbb] = d
  r[iaaaaab][iaaaab] = r[iabbbbb][iabbbb] = 5.0 * d
  r[iaaaaab][iaaaaaa] = r[iabbbbb][ibbbbbb] = c
  r[iaaaaab][iaaaabb] = r[iabbbbb][iaabbbb] = 5.0 * 1.0/5.0 * c
  
  # AAAABB goes to AAAAB, AAABB, AAAAAB, AAABBB
  r[iaaaabb][iaaaab] = r[iaabbbb][iabbbb] = 2.0 * d
  r[iaaaabb][iaaabb] = r[iaabbbb][iaabbb] = 4.0 * d
  r[iaaaabb][iaaaaab] = r[iaabbbb][iabbbbb] = 2.0 * 4.0/5.0 * c
  r[iaaaabb][iaaabbb] = r[iaabbbb][iaaabbb] = 4.0 * 2.0/5.0 * c
  
  # AAABBB goes to AAABB, AABBB, AAAABB, AABBBB
  r[iaaabbb][iaaabb] = r[iaaabbb][iaabbb] = 3.0 * d
  r[iaaabbb][iaaaabb] = r[iaaabbb][iaabbbb] = 3.0 * 3.0/5.0 * c

  # fill in diagonals as -(sum of column)
  for column in xrange(numstates):
    columntot = 0.0
    for row in xrange(numstates):
      if row != column:
        columntot += r[row][column]
    r[column][column] = -columntot
  
  sumdiag=0.0
  for i in xrange(numstates):
    sumdiag+=-r[i][i]

  return numstates/sumdiag ##Unweighted mean
 
#Converts str sequences to lists of tuples
def strToTouples(seq):
  seqTouples=list()
  for char in seq:
    seqTouples.append(charDict[char])
  return seqTouples

#Duplicate list of tuples
def applyGdTouple(seq):
  newseq=list()
  for tuple in seq:
    a=tuple[0]*2
    b=tuple[1]*2
    newseq.append((a,b))
  return newseq

#Converts tuple sequences to string sequences
def tuplesToStr(seq):
  newseq=""
  for tuple in seq:
    oa=tuple[0]
    ob=tuple[1]
    a=oa
    b=ob
    while a + b > 6: ##Keeps decreasing the state of the highest allele until the initial condition is met. If the a==b, it uses the original states to give priority, if oa==ob, reduces the two of them
      if a>b:
        a-=1
      elif b>a:
        b-=1
      else:# b==a
        if ob>oa:
          a-=1
        elif oa>ob:
          b-=1
        else:
            a-=1
            b-=1
#      if a > oa:
#        a-=1
#      if b > ob:
#        b-=1
    if oa + ob > 6:
      print("WARNING: Original state %d,%d, modified to fit in our model %d,%d" % (oa,ob,a,b))
    #if a + b > 6 :
      #print("DEBUG:This sequence cannot undergo GD, this state would be %d,%d" % (a,b)) 
      #raise Exception()
    newseq+=arrayChar[a][b]
  return newseq

#Duplicates Sequence
def applyGdSeq(seq):
  newseq=""
  for char in seq:
    oa=charDict[char][0]
    ob=charDict[char][1]
    a=oa*2
    b=ob*2
    while a + b > 6: ##Keeps decreasing the state of the highest allele until the initial condition is met. If the a==b, it uses the original states to give priority, if oa==ob, reduces the two of them
      if a>b:
        a-=1
      elif b>a:
        b-=1
      else:
        if ob>oa:
          a-=1
        elif oa>ob:
          b-=1
        else:
            a-=1
            b-=1
#      if a > oa:
#        a-=1
#      if b > ob:
#        b-=1
    if oa*2 + ob*2 > 6:
      print("WARNING: Original state %d,%d, duplicated state %d,%d, modified to fit in our model %d,%d" % (oa,ob,oa*2,ob*2,a,b))
    #if a + b > 6 :
      #print("DEBUG:This sequence cannot undergo GD, this state would be %d,%d" % (a,b)) 
      #raise Exception()
    newseq+=arrayChar[a][b]
  return newseq

#Calculate baseline mean
def baselineStr(seq):
  total=0
  for char in seq:
    total+=charDict[char][0] + charDict[char][1]
  return total/float(len(seq))

def baselineTouple(seq):
  total=0
  for tuple in seq:
    total+=tuple[0] + tuple[1]
  return total/float(len(seq))


def devFromBaselineModification(seq):
#  baseL=round(baselineStr(seq)/2.0,0)
  baseL=round(baselineTouple(seq)/2.0,0)
#  newseq=""
  newseq=list()
#  for char in seq:
  for tuple in seq:
    na=1
    nb=1
#    if charDict[char][0] > baseL :
    if tuple[0] > baseL :
      na=2
#   elif charDict[char][0] < baseL :
    elif tuple[0] < baseL :
      na=0
#   if charDict[char][1] > baseL :
    if tuple[1] > baseL :
      nb=2
#   elif charDict[char][1] < baseL :
    elif tuple[1] < baseL :
      nb=0
#    newseq+=()arrayChar[na][nb]
    newseq.append((na,nb))
#    if arrayChar[na][nb] != char:
#      print "DEBUG: %d %s %s" % (baseL,char,arrayChar[na][nb])
  return newseq

def relativeFromBaselineModification(seq,rtuples):
#   baseL=int(round(baselineStr(seq)/2.0,0))
  baseL=int(round(baselineTouple(seq)/2.0,0))
#  newseq=""
  newseq=list()
#  for char in seq:
  for position in xrange(len(seq)):
    tuple=seq[position]
    rtup=rtuples[position]
#    na=charDict[char][0]
#    nb=charDict[char][1]
    na=tuple[0]
    nb=tuple[1]
    ra=rtup[0]
    rb=rtup[1]
    #A copy
    if na % baseL == 0:
      na=na/baseL
    elif ra < 0.5:
      na=int(math.floor(na/float(baseL)))
      if na==0:
        na=1 #There was signal, and therefore we assume that at least there have to be one copy
    else:
      na=int(math.ceil(na/float(baseL)))
    #B copy
    if nb % baseL == 0:
      nb=nb/baseL
    elif rb < 0.5:
      nb=int(math.floor(nb/float(baseL)))
      if nb==0:
        nb=1 #There was signbl, and therefore we assume that at least there have to be one copy
    else:
      nb=int(math.ceil(nb/float(baseL)))

#    newseq+=arrayChar[na][nb]
    newseq.append((na,nb))
#    if arrayChar[na][nb] != char:
#      print "DEBUG: %d %s %s" % (baseL,char,arrayChar[na][nb])
  return newseq 

# draw an exponential waiting time
def drawtime(prob):
  if prob == 0.0:
    return failnum
  u = random.random()
  tyme = -math.log(u)/prob
  return tyme

#draw an exponential waiting time with an exponentially growing rate
def drawtimeExp(prob,growthRate,initTime):
  if prob == 0.0 or growthRate == 0.0:
    return failnum
  u = random.random()
  absoluteTime = math.log(-math.log(u)/prob * growthRate + math.exp(growthRate * initTime))/growthRate
#  print "DEBUG: normal WGD time would be %f, but with the exponential growth rate it is %f\n" % (-math.log(u)/prob,absoluteTime-initTime)
  return absoluteTime - initTime

# one contender in the "horserace"; knows its outcome if chosen

class Horse:
  def __init__(self, prob, result):
    self.prob = prob
    self.result = result
  def gettime(self, endtime):
    trytime = drawtime(self.prob)
    if trytime > endtime: 
      return failnum
    else:
      return trytime
  def getresult(self):
    return self.result

def besthorse(horses,endtime,curstate):
  times = []
  besttime = bigtime
  newstate = failnum
  random.shuffle(horses)   # random tiebreaks!
  for h in horses:
    newtime = h.gettime(endtime)
    if newtime != failnum:
      if newtime < besttime:
        besttime = newtime
        newstate = h.getresult()
  if besttime == bigtime:
    return (endtime, curstate)
  else:
    return (besttime, newstate)
  
def nextstate(curstate, endtime):
  horses = []

  if curstate == null:   # null; no horserace as outcome is fixed
    return (endtime,curstate)

  elif curstate == a: # A
    horses.append(Horse(d,null))   # null
    horses.append(Horse(g,aa))   # AA

  elif curstate == b: # B
    horses.append(Horse(d,null))   # null  
    horses.append(Horse(g,bb))   # BB

  elif curstate == aa: # AA
    horses.append(Horse(2.0*d,a))  # A
    horses.append(Horse(2.0*g,aaa))  # AAA

  elif curstate == ab: # AB
    horses.append(Horse(d,a))   # A
    horses.append(Horse(d,b))   # B
    horses.append(Horse(c,aa))   # AA
    horses.append(Horse(c,bb))   # BB
    horses.append(Horse(g,aab))   # AAB
    horses.append(Horse(g,abb))   # ABB
      
  elif curstate == bb: # BB
    horses.append(Horse(2.0*d,b))  # B
    horses.append(Horse(2.0*g,bbb))  # BBB

  elif curstate == aaa: # AAA
    horses.append(Horse(3.0*d,aa))   # AA
    horses.append(Horse(3.0*g,aaaa))  # AAAA

  elif curstate == aab: # AAB
    horses.append(Horse(d,aa))       # AA
    horses.append(Horse(2.0*d,ab))   # AB
    horses.append(Horse(c,aaa))       # AAA
    horses.append(Horse(2.0*0.5*c,abb)) # ABB 
    horses.append(Horse(2.0*g,aaab))  # AAAB
    horses.append(Horse(g,aabb))      # AABB

  elif curstate == abb: # ABB
    horses.append(Horse(2.0*d,ab))   # AB
    horses.append(Horse(d,bb))       # BB
    horses.append(Horse(2.0*0.5*c,aab))  # AAB
    horses.append(Horse(c,bbb))       # BBB
    horses.append(Horse(g,aabb))      # AABB
    horses.append(Horse(2.0*g,abbb))  # ABBB

  elif curstate == bbb: # BBB
    horses.append(Horse(3.0*d,bb))   # BB
    horses.append(Horse(3.0*g,bbbb))  # BBBB

  elif curstate == aaaa: # AAAA
    horses.append(Horse(4.0*d,aaa))   # AAA
    horses.append(Horse(4.0*g,aaaaa))  # AAAAA
    
  elif curstate == aaab: # AAAB
    horses.append(Horse(d,aaa))       # AAA
    horses.append(Horse(3.0*d,aab))   # AAB
    horses.append(Horse(c,aaaa))      # AAAA
    horses.append(Horse(3.0*1.0/3.0*c,aabb))  # AABB
    horses.append(Horse(3.0*g,aaaab))  # AAAAB
    horses.append(Horse(g,aaabb))      # AAABB
    
  elif curstate == aabb: # AABB
    horses.append(Horse(2.0*d,aab))   # AAB
    horses.append(Horse(2.0*d,abb))   # ABB
    horses.append(Horse(2.0*2.0/3.0*c,aaab))  # AAAB
    horses.append(Horse(2.0*2.0/3.0*c,abbb))  # ABBB
    horses.append(Horse(2.0*g,aaabb))  # AAABB
    horses.append(Horse(2.0*g,aabbb))  # AABBB
    
  elif curstate == abbb: # ABBB
    horses.append(Horse(3.0*d,abb))   # ABB
    horses.append(Horse(d,bbb))       # BBB
    horses.append(Horse(3.0*1.0/3.0*c,aabb))  # AABB
    horses.append(Horse(c,bbbb))      # BBBB
    horses.append(Horse(g,aabbb))      # AABBB
    horses.append(Horse(3.0*g,abbbb))  # ABBBB

  elif curstate == bbbb: # BBBB
    horses.append(Horse(4.0*d,bbb))   # BBB
    horses.append(Horse(4.0*g,bbbbb))  # BBBBB

  elif curstate == aaaaa: # AAAAA
    horses.append(Horse(5.0*d,aaaa))  # AAAA
    horses.append(Horse(5.0*g,aaaaaa))  # AAAAAA

  elif curstate == aaaab: # AAAAB
    horses.append(Horse(d,aaaa))      # AAAA
    horses.append(Horse(4.0*d,aaab))  # AAAB
    horses.append(Horse(c,aaaaa))      # AAAAA
    horses.append(Horse(4.0*1.0/4.0*c,aaabb))  # AAABB
    horses.append(Horse(4.0*g,aaaaab))  # AAAAAB
    horses.append(Horse(g,aaaabb))      # AAAABB

  elif curstate == aaabb: # AAABB
    horses.append(Horse(2.0*d,aaab))  # AAAB
    horses.append(Horse(3.0*d,aabb))  # AABB
    horses.append(Horse(2.0*3.0/4.0*c,aaaab))  # AAAAB
    horses.append(Horse(3.0*2.0/4.0*c,aabbb))  # AABBB
    horses.append(Horse(3.0*g,aaaabb))  # AAAABB
    horses.append(Horse(2.0*g,aaabbb))  # AAABBB

  elif curstate == aabbb: # AABBB
    horses.append(Horse(3.0*d,aabb))  # AABB
    horses.append(Horse(2.0*d,abbb))  # ABBB
    horses.append(Horse(3.0*2.0/4.0*c,aaabb))  # AAABB
    horses.append(Horse(2.0*3.0/4.0*c,abbbb))  # ABBBB
    horses.append(Horse(2.0*g,aaabbb))  # AAABBB
    horses.append(Horse(3.0*g,aabbbb))  # AABBBB

  elif curstate == abbbb: # ABBBB
    horses.append(Horse(4.0*d,abbb))  # ABBB
    horses.append(Horse(d,bbbb))      # BBBB
    horses.append(Horse(4.0*1.0/4.0*c,aabbb))  # AABBB
    horses.append(Horse(c,bbbbb))      # BBBBB
    horses.append(Horse(g,aabbbb))      # AABBBB
    horses.append(Horse(4.0*g,abbbbb))  # ABBBBB

  elif curstate == bbbbb: # BBBBB
    horses.append(Horse(5.0*d,bbbb))  # BBBB
    horses.append(Horse(5.0*g,bbbbbb))  # BBBBBB

  elif curstate == aaaaaa: # AAAAAA
    horses.append(Horse(6.0*d,aaaaa))  # AAAAA

  elif curstate == aaaaab: # AAAAAB
    horses.append(Horse(d,aaaaa))      # AAAAA
    horses.append(Horse(5*d,aaaab))    # AAAAB
    horses.append(Horse(c,aaaaaa))      # AAAAAA
    horses.append(Horse(5.0*1.0/5.0*c,aaaabb))  # AAAABB

  elif curstate == aaaabb: # AAAABB
    horses.append(Horse(2.0*d,aaaab))  # AAAAB
    horses.append(Horse(4.0*d,aaabb))  # AAABB
    horses.append(Horse(2.0*4.0/5.0*c,aaaaab))  # AAAAAB
    horses.append(Horse(4.0*2.0/5.0*c,aaabbb))  # AAABBB

  elif curstate == aaabbb: # AAABBB
    horses.append(Horse(3.0*d,aaabb))  # AAABB
    horses.append(Horse(3.0*d,aabbb))  # AABBB
    horses.append(Horse(3.0*3.0/5.0*c,aaaabb))  # AAAABB
    horses.append(Horse(3.0*3.0/5.0*c,aabbbb))  # AABBBB

  elif curstate == aabbbb: # AABBBB
    horses.append(Horse(4.0*d,aabbb))  # AABBB
    horses.append(Horse(2.0*d,abbbb))  # ABBBB
    horses.append(Horse(4.0*2.0/5.0*c,aaabbb))  # AAABBB
    horses.append(Horse(2.0*4.0/5.0*c,abbbbb))  # ABBBBB

  elif curstate == abbbbb: # ABBBBB
    horses.append(Horse(5.0*d,abbbb))  # ABBBB
    horses.append(Horse(d,bbbbb))      # BBBBB
    horses.append(Horse(5.0*1.0/5.0*c,aabbbb))  # AABBBB
    horses.append(Horse(c,bbbbbb))      # BBBBBB

  elif curstate == bbbbbb: # BBBBBB
    horses.append(Horse(6.0*d,bbbbb))  # BBBBB

  else:
    print "Illegal state",curstate,"found"
    exit()

  return besthorse(horses,endtime,curstate)

def simbranch(endtime,numloci):
  answers = []
  for i in xrange(numloci):
    curtime = 0.0
    curstate = 4
    while curtime < endtime:
      (curtime,curstate) = nextstate(curstate,endtime)
    answers.append(curstate)
  return answers

def gendouble(curstate):
  if curstate == null:  return null
  if curstate == a:  return aa
  if curstate == b:  return bb
  if curstate == aa:  return aaaa
  if curstate == ab:  return aabb
  if curstate == bb:  return bbbb
  if curstate == aaa:  return aaaaaa
  if curstate == aab:  return aaaabb
  if curstate == abb:  return aabbbb
  if curstate == aaaa:  return aaaaaa
  if curstate == aaab: 
    if random.random() < 0.5:  return aaaabb
    else:                      return aaaaab
  if curstate == aabb:  return aaabbb
  if curstate == abbb:
    if random.random() < 0.5:  return aabbbb
    else:                      return abbbbb
  if curstate == bbbb:  return bbbbbb
  if curstate == aaaaa:  return aaaaaa
  if curstate == aaaab:  
    if random.random() < 0.8:  return aaaaab
    else:                      return aaaabb
  if curstate == aaabb:
    if random.random() < 0.6:  return aaaabb
    else:                      return aaabbb
  if curstate == aabbb:
    if random.random() < 0.4:  return aaabbb
    else:                      return aabbbb
  if curstate == abbbb:
    if random.random() < 0.2:  return aabbbb
    else:                      return abbbbb
  if curstate == bbbbb:  return bbbbbb
  return curstate   # no adjustment for 6-copy states

# simulate evolution along a branch (from root end
# to tip end).  
# brlen is the length of the branch
# seq is the parental sequence
# alreadydoubled indicates whether the parent had a genome
# duplication already (in which case no more can occur).
def simbranchinseq(brlen,seq,alreadydoubled,initTime,beta):
  wasdoubled = alreadydoubled
  if not alreadydoubled and pGD>0.0:
    gdtime = drawtimeExp(pGD,gGD,initTime) ###Working here, I need to keep accumulating branchlenghts to know the initTime
  else:
    gdtime = bigtime
  if gdtime < brlen:   # genome doubling happens
    #print "Genome Doubling!"
    wasdoubled = True
    answers = simbranchinseq_segment(gdtime*beta,seq) ##Beta scales the branch lengths for the mutational process since the matrix has not been scaled
    for i in xrange(len(answers)):
      answers[i] = gendouble(answers[i])
    answers = simbranchinseq_segment((brlen - gdtime)*beta,answers)
  else:
    answers = simbranchinseq_segment(brlen*beta,seq)
  return (answers, wasdoubled)
    
def simbranchinseq_segment(endtime,seq):
  answers = []
  for curstate in seq:
    curtime = 0.0
    while curtime < endtime:
      (curtime,curstate) = nextstate(curstate,endtime)
    answers.append(curstate)
  return answers

def simseqtree(numloci,tree,beta):
  root_states=[ab]*numloci
  tree.calc_node_root_distances(return_leaf_distances_only=False)
  initTime=0.0
  for edge in tree.preorder_edge_iter():
    node=edge.head_node
    if not hasattr(node, "sequence"):
      setattr(node, "sequence", [])
    seq_list = getattr(node, "sequence")
    if edge.tail_node:
      par = edge.tail_node
      if len(seq_list) != n_prev_seq:
      	raise ValueError("'%s' length varies among nodes" % par.sequence)
      par_seq = getattr(par, "sequence")[-1]#-1??
      length = getattr(edge, "length")
      #mutation_rate = getattr(edge, self.edge_rate_attr, None) or self.mutation_rate
      assert hasattr(par, "doubled")
      wasdoubled = getattr(par, "doubled")
      initTime=getattr(par,"root_distance")
#      print "DEBUG: branch start %f, branchlength %f\n" %(initTime,length)
      (newseq,wasdoubled) = simbranchinseq(length,par_seq,wasdoubled,initTime,beta)
      setattr(node, "doubled", wasdoubled)
      seq_list.append(newseq)
    else:
       # no tail node: root
      setattr(node, "doubled", False)
      n_prev_seq = len(seq_list)
      if root_states is not None:
        seq_list.append(root_states)
      else:
        assert n_prev_seq > 0
        n_prev_seq -= 1
    
def readnexus(input):
  names = []
  seqs = []
  beforematrix = True
  inmatrix = False
  aftermatrix = False
  header = []
  footer = []
  for line in input:
    if beforematrix:   # save the header
      header.append(line)
      if line.strip() == "MATRIX":
        beforematrix = False
        inmatrix = True
      continue
    elif inmatrix:     # parse the data
      sline = line.rstrip().split()
      if sline[0] == ";":  
        footer.append(line)
        inmatrix = False
        aftermatrix = True
        continue
      else:
        names.append(sline[0])
        seqs.append(sline[1])
    else:              # save the footer
      footer.append(line)
  return (names, seqs, header, footer)


#Config XML output
beast_outname="beast.out"

# MAIN PROGRAM
##############

parser = argparse.ArgumentParser(description="Simulates evolution down a tree in newick format according to our SCA model. Up to 6 copies")
parser.add_argument("-c",type=float,default=0,required=True,help="Relative conversion rate")
parser.add_argument("-d",type=float,default=0,required=True,help="Relative loss rate")
parser.add_argument("-g",type=float,default=0,required=True,help="Relative gain rate")
parser.add_argument("-r",type=float,default=0,required=True,help="Mutation rate per allele and fragment")
parser.add_argument("-pgd",type=float,default=0,required=True,help="Whole genome duplication rate")
parser.add_argument("-ggd",type=float,default=0,required=True,help="Exponential growth rate of the whole genome duplication rate")
parser.add_argument("-i",type=str,default="infile.tree",required=True,help="Input Newick tree file")
parser.add_argument("-o",type=str,default="outtree.nex",required=False,help="Output nexus file")
parser.add_argument("-n",type=int,default=1,required=True,help="Number of loci")
parser.add_argument("--xml",action='store_true',help="XML alignment element for BEAST")
parser.set_defaults(xml=False)
parser.add_argument("--randomGD",type=int,default=0,required=False,help="Number of random whole genome doubling eevents to add in the tips")
parser.add_argument("--seed",type=float,default=15,required=False,help="Seed for the pseudo-random number generator")
parser.add_argument("--ngen",type=int,default=100000000,required=False,help="Number of MCMC generations in the XML BEAST output file")
parser.add_argument("--period",type=int,default=10000,required=False,help="Sampling period for the MCMC chain")
parser.add_argument("-a",type=float,default=40,required=False,help="Assumed age of the patient when the cenancestor arose")
parser.add_argument("--mod",type=str,default="none",required=False,help="Modification of the sequences. So far we implemented, none: no modification, baseline: relative to the baseline, max2: 9 states, with a maximum of 2 copies per cromosome")
parser.add_argument("--exp",action="store_true",help="Coalescent model with exponential growth")
parser.set_defaults(exp=False)
args = parser.parse_args()

#Reclycing variables
trate=args.r
trelrate=args.c+args.d+args.g
 
#Normalized relative rates to use in the matrix
c=args.c/trelrate*trate
d=args.d/trelrate*trate
g=args.g/trelrate*trate
numloci = args.n

pGD = args.pgd # probability of genome doubling per unit time
gGD = args.ggd # exponential growth rate of the genome doubling rate

#Random init
random.seed(args.seed)

#Sequence simulation using our model along a rooted input newick tree
tree=dendropy.Tree.get(path=args.i,rooting="force-rooted",schema="newick")

##Calculate the scaling factor B
beta=getScalingFactor()
#print(beta)
#exit()

#simseqtree(numloci,tree,beta)
simseqtree(numloci,tree,1) ###We are not scaling the matrix since a) we are not using a common rate, b) if we did, this may generate a lot of estimation error since we would be assuming equal equilibrium frequencies, where we don't even have.

#Matrix which will obtain the data still pending from the tree
char_matrix= dendropy.StandardCharacterMatrix(taxon_namespace=tree.taxon_namespace,default_state_alphabet=None)
char_matrix.taxon_namespace = tree.taxon_namespace

#Object of class with methods to harvest the sequences from the tree
seq_evolver= dendropy.model.discrete.DiscreteCharacterEvolver(seq_attr="sequence")
seq_evolver.extend_char_matrix_with_characters_on_tree(char_matrix=char_matrix,tree=tree)

#Modification of sequences a posteriori
#######################################
names=char_matrix.taxon_namespace.labels()
seqs=char_matrix.sequences()
doublednames = []
for tip in tree.leaf_node_iter():
  if getattr(tip,"doubled"):
    doublednames.append(tip.taxon.label)

##Removing the outgroup
re_outgroup=re.compile('outgroup')
for i in xrange(len(names)):
  if re_outgroup.match(names[i]):
    del seqs[i]
    del names[i]

#Conversion of sequences to lists of tuples
for i in xrange(len(seqs)):
  strSeq=seqs[i].symbols_as_string(sep="")
  seqs[i]=strToTouples(strSeq)

#Adding random GD a posteriori
if args.randomGD > 0:
  it=0
  ti=0
  while it<args.randomGD:
    i=random.randint(0,len(names)-1)
#    print("DEBUG: name %d, len names %d" % (i,len(names)))
    if names[i] in doublednames:
      it-=1
      print "Random duplication in a duplicated leaf, trying again"
    else:
      try:
        seqs[i]=applyGdTouple(seqs[i])
#        seqs[i]=applyGdSeq(seqs[i])
        doublednames.append(names[i])
      except Exception as inst:
        it-=1
        print "Random duplication in a leaf that cannot undergo GD"
    if ti >= maxit:
      raise Exception("Infinite loop")
    it+=1
    ti+=1

#Recoding the states to tackle GD
#Sequences as list of tuples are the real states. Even the original-sequence option needs some recoding, since the random duplication of tips generates states not included in our model
newseqs=list()
if args.mod == "none":
  print "Original sequences"
  for i in xrange(len(seqs)):##Unnecessary, just for debuggin purposes
    #seqsstringSeq=seqs[i].symbols_as_string(sep="")
    newseqs.append(tuplesToStr(seqs[i]))
elif args.mod == "max2":
  print "Sequences modified to the 9-state model"
  for i in xrange(len(seqs)):
    #stringSeq=seqs[i].symbols_as_string(sep="")
    newseqs.append(tuplesToStr(devFromBaselineModification(seqs[i])))
elif args.mod == "baseline":
  print "States relative to the baseline"
  r=[(random.random(),random.random()) for i in xrange(len(seqs[0]))]
  for i in xrange(len(seqs)):
    #stringSeq=seqs[i].symbols_as_string(sep="")
    newseqs.append(tuplesToStr(relativeFromBaselineModification(seqs[i],r)))
else:
  raise Exception('Incorrect --mod option')
seqs=newseqs

##Removing invariable sites (H columns). Quite inefficient
toremove=list()
for i in xrange(len(seqs[0])):
  remove=1
  for seq in seqs:
    if seq[i]!='H':
      remove=0
      break
  toremove.append(remove)
newseqs=list()
for seq in seqs:
  newseq=str()
  for i in xrange(len(seq)):
    if toremove[i]==0:
      newseq+=seq[i]
  newseqs.append(newseq)

##Obtainint dates for the XML output
if args.xml==True:
  diff_date=0.0
  max_date=0.0
  min_date=bigtime
  heights=tree.calc_node_root_distances()
  for date in heights:
    currHeight=round(date + args.a,0)
    if currHeight>max_date:
      max_date=currHeight
    elif currHeight<min_date:
      min_date=currHeight
  diff_date=max_date-min_date
#print("DEBUG: max_date=" + str(max_date) + ", min_date=" + str(min_date) + "\n")

##Output
outfile = open(args.o,"w")
if args.xml==False:
  outfile.write("#NEXUS\n")
  outfile.write("\n")
  outfile.write("BEGIN TAXA;\n")
  outline = "\tDIMENSIONS NTAX="
  outline += str(len(names)) + ";\n"
  outfile.write(outline)
  outfile.write("\tTAXLABELS\n")
  for name in names:
    if name in doublednames:
      name += "*"
    outfile.write("\t" + name + "\n")
  outfile.write("  ;\n")
  outfile.write("END")
  outfile.write("\n")
  outfile.write("BEGIN CHARACTERS;\n")
  outline = "\tDIMENSIONS NCHAR="
  outline += str(len(newseqs[0])) + ";\n"
  outfile.write(outline)
  outfile.write('\tFORMAT DATATYPE=STANDARD SYMBOLS="";\n')
  outfile.write("\tMATRIX\n")
else:
  #
  namespace=tree.taxon_namespace
  labels=namespace.labels()
  labels.remove("outgroup")
  outfile.write("<?xml version=\"1.0\" standalone=\"yes\"?>\n<beast>\n<taxa id=\"taxa\">\n")
  #Get labels and discard the outgroup label
  for label in labels:
    if label in doublednames:
      label_to_print=label + "*"
    else:
      label_to_print=label
    date= args.a + tree.find_node_for_taxon(namespace.get_taxon(label)).distance_from_root()#Date equals the distance from the root plus a constant
    outfile.write("\t<taxon id=\"" + label_to_print + "\">\n\t\t<date value=\"" + str(round(date,0)) + "\" direction=\"forwards\" units=\"years\"/>\n\t</taxon>\n")
  outfile.write("</taxa>")
  outfile.write('\n\n<generalDataType id="cnv">\n\t<state code="@"/> <!-- Genotype: 0,0 ; Beast State: 0 -->\n\t<state code="A"/> <!-- Genotype: 0,1 ; Beast State: 1 -->\n\t<state code="B"/> <!-- Genotype: 0,2 ; Beast State: 2 -->\n\t<state code="C"/> <!-- Genotype: 0,3 ; Beast State: 3 -->\n\t<state code="D"/> <!-- Genotype: 0,4 ; Beast State: 4 -->\n\t<state code="E"/> <!-- Genotype: 0,5 ; Beast State: 5 -->\n\t<state code="F"/> <!-- Genotype: 0,6 ; Beast State: 6 -->\n\t<state code="G"/> <!-- Genotype: 1,0 ; Beast State: 7 -->\n\t<state code="H"/> <!-- Genotype: 1,1 ; Beast State: 8 -->\n\t<state code="I"/> <!-- Genotype: 1,2 ; Beast State: 9 -->\n\t<state code="J"/> <!-- Genotype: 1,3 ; Beast State: 10 -->\n\t<state code="K"/> <!-- Genotype: 1,4 ; Beast State: 11 -->\n\t<state code="L"/> <!-- Genotype: 1,5 ; Beast State: 12 -->\n\t<state code="M"/> <!-- Genotype: 2,0 ; Beast State: 13 -->\n\t<state code="N"/> <!-- Genotype: 2,1 ; Beast State: 14 -->\n\t<state code="O"/> <!-- Genotype: 2,2 ; Beast State: 15 -->\n\t<state code="P"/> <!-- Genotype: 2,3 ; Beast State: 16 -->\n\t<state code="Q"/> <!-- Genotype: 2,4 ; Beast State: 17 -->\n\t<state code="R"/> <!-- Genotype: 3,0 ; Beast State: 18 -->\n\t<state code="S"/> <!-- Genotype: 3,1 ; Beast State: 19 -->\n\t<state code="T"/> <!-- Genotype: 3,2 ; Beast State: 20 -->\n\t<state code="U"/> <!-- Genotype: 3,3 ; Beast State: 21 -->\n\t<state code="V"/> <!-- Genotype: 4,0 ; Beast State: 22 -->\n\t<state code="W"/> <!-- Genotype: 4,1 ; Beast State: 23 -->\n\t<state code="X"/> <!-- Genotype: 4,2 ; Beast State: 24 -->\n\t<state code="Y"/> <!-- Genotype: 5,0 ; Beast State: 25 -->\n\t<state code="Z"/> <!-- Genotype: 5,1 ; Beast State: 26 -->\n\t<state code="["/> <!-- Genotype: 6,0 ; Beast State: 27 -->\n\t<ambiguity code="-" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>\n\t<ambiguity code="?" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>\n</generalDataType>\n')
  beast_outname=os.path.basename(args.o)
  beast_outname=re.sub(".temp","",beast_outname)
  outfile.write("<alignment id=\"alignment\">\n\t<dataType idref=\"cnv\"/>\n")

for name,seq in zip(names,newseqs):
  if name in doublednames:
    name += "*"
  if args.xml==False:
    outline = "\t\t" + name + "\t" + seq + "\n"
  else:
    outline = "\t<sequence>\n\t\t<taxon idref=\"" + name + "\"/>\n\t\t" + seq + "\n\t</sequence>\n"
  outfile.write(outline)

if args.xml==False:
  outfile.write("\t;\n")
  outfile.write("END;\n")
else:
  outfile.write('''</alignment>
<ascertainedCharacterPatterns id="patterns">
	<alignment idref="alignment"/>
	<state code='H'/>
</ascertainedCharacterPatterns>

<!--
<patterns id="patterns" from="1" strip="false">
	<alignment idref="alignment"/>
</patterns>
-->''' + ('''
<!-- A prior assumption that the population size has remained constant       -->
<!-- throughout the time spanned by the genealogy.                           -->
<constantSize id="constant" units="years">
	<populationSize>
		<parameter id="constant.popSize" value="1" lower="0.0"/>
	</populationSize>
</constantSize>
''','''
<exponentialGrowth id="exponential" units="years">
        <populationSize>
                <parameter id="exponential.popSize" value="1" lower="0.0"/>
        </populationSize>
        <growthRate>
                <parameter id="exponential.growthRate" value="3.0E-4"/>
        </growthRate>
</exponentialGrowth>
''')[args.exp==True] + '''
<!-- Generate a random starting tree under the coalescent process            -->
<coalescentSimulator id="startingTree">
	<taxa idref="taxa"/>
	''' + ('''<constantSize idref="constant"/>''','''<exponentialGrowth idref="exponential"/>''')[args.exp==True] +'''
</coalescentSimulator>

<!-- Generate a tree model                                  -->
<treeModel id="treeModel">
	<coalescentTree idref="startingTree"/>
	<rootHeight>
		<parameter id="treeModel.rootHeight"/>
	</rootHeight>
	<nodeHeights internalNodes="true">
		<parameter id="treeModel.internalNodeHeights"/>
	</nodeHeights>
	<nodeHeights internalNodes="true" rootNode="true">
		<parameter id="treeModel.allInternalNodeHeights"/>
	</nodeHeights>
</treeModel>

<!-- Generate a coalescent likelihood                                        -->
<coalescentLikelihood id="coalescent">
	<model>
		''' +('''<constantSize idref="constant"/>''','''<exponentialGrowth idref="exponential"/>''')[args.exp==True] + '''
	</model>
	<populationTree>
		<treeModel idref="treeModel"/>
	</populationTree>
</coalescentLikelihood>

<!-- The strict clock (Uniform rates across branches)                        -->

<strictClockCenancestorBranchRates id="branchRates">
	<rate>
		<parameter id="clock.rate" value="1"/>
	</rate>
</strictClockCenancestorBranchRates>

<frequencyModel id="frequencies">
	<dataType idref="cnv"/>
	<frequencies>
		<parameter id="cnv.frequencies" value="0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"/>
	</frequencies>
</frequencyModel>

<CNVModel id="cnv_subsmodel">
	<frequencies>
		<frequencyModel idref="frequencies"/>
	</frequencies>
	<gain_rate>
		<parameter id="cnv.gain" value="1" lower="0"/>
	</gain_rate>
	<loss_rate>
		<parameter id="cnv.loss" value="1" lower="0"/>
	</loss_rate>
	<conversion_rate>
		<parameter id="cnv.conversion" value="1" lower="0"/>
	</conversion_rate>
</CNVModel>

<siteModel id="siteModel">
	<substitutionModel>
		<CNVModel idref="cnv_subsmodel"/>
	</substitutionModel>
</siteModel>

<cenancestorTreeLikelihood id="treeLikelihood" useAmbiguities="false">
	<patterns idref="patterns"/>
	<treeModel idref="treeModel"/>
	<siteModel idref="siteModel"/>
	<cenancestorHeight>
		<parameter id="luca_height" lower="''' + str(diff_date) + '''" upper="''' + str(max_date) + '''"/>
	</cenancestorHeight>
	<cenancestorBranch>
		<parameter id="luca_branch" value="1" upper="''' + str(min_date) + '''" lower="0.0"/>
		<!-- Value 1 as a safe starting value -->
	</cenancestorBranch>
	<strictClockCenancestorBranchRates idref="branchRates"/>
</cenancestorTreeLikelihood>

<operators id="operators" optimizationSchedule="default">
	<scaleOperator scaleFactor="0.25" weight="0.25">
                <parameter idref="cnv.loss"/>
	</scaleOperator>
	<scaleOperator scaleFactor="0.25" weight="0.25">
                <parameter idref="cnv.conversion"/>
	</scaleOperator>
	<scaleOperator scaleFactor="0.5" weight="10.0">
                <parameter idref="clock.rate"/>
	</scaleOperator>
	<subtreeSlide size="2.5" gaussian="true" weight="15.0"> <!-- 2.5 years. They will be automatically optimized by BEAST though -->
		<treeModel idref="treeModel"/>
	</subtreeSlide>
	<narrowExchange weight="15.0">
		<treeModel idref="treeModel"/>
	</narrowExchange>
	<wideExchange weight="3.0">
		<treeModel idref="treeModel"/>
	</wideExchange>
	<wilsonBalding weight="3.0">
		<treeModel idref="treeModel"/>
	</wilsonBalding>
	<scaleOperator scaleFactor="0.75" weight="5.0">
		<parameter idref="treeModel.rootHeight"/>
	</scaleOperator>
	<uniformOperator weight="30.0">
		<parameter idref="treeModel.internalNodeHeights"/>
	</uniformOperator>
	
	<scaleOperator scaleFactor="0.2" weight="1.0"> <!-- We operate the branch since it is relative to the root. Operating luca_height is error prone, since it depends on the root -->
                <parameter idref="luca_branch"/>
        </scaleOperator>
''' + ('''
	<scaleOperator scaleFactor="0.5" weight="3.0">
		<parameter idref="constant.popSize"/>
	</scaleOperator>
''','''
	<scaleOperator scaleFactor="0.75" weight="3">
                <parameter idref="exponential.popSize"/>
        </scaleOperator>

        <randomWalkOperator windowSize="1.0" weight="3">
                <parameter idref="exponential.growthRate"/>
        </randomWalkOperator>
''')[args.exp==True] + '''
        <upDownOperator scaleFactor="0.75" weight="5.0">
                <up>
                        <parameter idref="clock.rate"/>
                </up>
                <down>
                        <parameter idref="treeModel.allInternalNodeHeights"/>
                </down>
        </upDownOperator>

</operators>

<!-- Define MCMC                                                             -->
<mcmc id="mcmc" chainLength="''' + str(args.ngen) + '''" autoOptimize="true" operatorAnalysis="''' + beast_outname + '''.ops">
	<posterior id="posterior">
		<prior id="prior">
                        <coalescentLikelihood idref="coalescent"/>''' + ('''
			<oneOnXPrior>
				<parameter idref="constant.popSize"/>
			</oneOnXPrior>''','''
			<oneOnXPrior>
                                <parameter idref="exponential.popSize"/>
                        </oneOnXPrior>
                        <exponentialPrior mean="0.01" offset="0.0">
                                <parameter idref="exponential.growthRate"/>
                        </exponentialPrior>
''')[args.exp==True] + '''
		
			<!-- Clock (gain) Rate Prior. More than 50 SGAs/breakpoint/year seems an unreasonable enough value to use as upper bound-->
			<!-- 
			<uniformPrior lower="0.0" upper="50">
				<parameter idref="clock.rate"/>
			</uniformPrior>
			-->
			<!-- Clock (gain) Rate Prior. -->
			<logNormalPrior mean="-4.0" stdev="2.5" offset="0.0" meanInRealSpace="false">
                        	<parameter idref="clock.rate"/>
                        </logNormalPrior>

			<!-- Loss and conversion (relative to gain) rate priors-->
			<exponentialPrior mean="1.0" offset="0.0">
				<parameter idref="cnv.loss"/>
			</exponentialPrior>
			<exponentialPrior mean="1.0" offset="0.0">
				<parameter idref="cnv.conversion"/>
			</exponentialPrior>

                        <!-- Cenancestor Prior on the height, since it is easier to have a meaningful prior on it (time of the initial development of the BE fragment) -->
                        <uniformPrior lower="''' + str(diff_date) + '''" upper="''' + str(max_date) + '''">
                        	<parameter idref="luca_height"/>
                        </uniformPrior>
		</prior>
		<likelihood id="likelihood">
			<cenancestorTreeLikelihood idref="treeLikelihood"/>
		</likelihood>
	</posterior>
	<operators idref="operators"/>

	<!-- write log to screen                                                     -->
	<log id="screenLog" logEvery="''' + str(args.period) + '''">
		<column label="Posterior" dp="4" width="12">
			<posterior idref="posterior"/>
		</column>
		<column label="Prior" dp="4" width="12">
			<prior idref="prior"/>
		</column>
		<column label="Likelihood" dp="4" width="12">
			<likelihood idref="likelihood"/>
		</column>
		<column label="rel_loss_rate" sf="6" width="12">
			<parameter idref="cnv.loss"/>
		</column>
		<column label="rel_conv_rate" sf="6" width="12">
			<parameter idref="cnv.conversion"/>
		</column>
		<column label="gain_rate" sf="6" width="12">
			<parameter idref="clock.rate"/>
		</column>

		<column label="rootHeight" sf="6" width="12">
			<parameter idref="treeModel.rootHeight"/>
		</column>
		
		<column label="luca_height" sf="6" width="12">
			<parameter idref="luca_height"/>
		</column>
		
		<column label="luca_branch" sf="6" width="12">
			<parameter idref="luca_branch"/>
		</column>
	</log>

	<!-- write log to file                                                       -->
	<log id="fileLog" logEvery="''' + str(args.period) + '''" fileName="''' + beast_outname + '''.log" overwrite="false">
		<posterior idref="posterior"/>
		<prior idref="prior"/>
		<likelihood idref="likelihood"/>
		<parameter idref="cnv.loss"/>
		<parameter idref="cnv.conversion"/>
		<parameter idref="treeModel.rootHeight"/>
		<parameter idref="luca_height"/>
		<parameter idref="luca_branch"/>''' + ('''
		<parameter idref="constant.popSize"/>''','''
		<parameter idref="exponential.popSize"/>
                <parameter idref="exponential.growthRate"/>''')[args.exp==True] + '''
		<parameter idref="clock.rate"/>
		<coalescentLikelihood idref="coalescent"/>
	</log>

	<!-- write tree log to file                                                  -->
	<logTree id="treeFileLog" logEvery="''' + str(args.period) + '''" nexusFormat="true" fileName="''' + beast_outname + ".trees" + '''" sortTranslationTable="true">
		<treeModel idref="treeModel"/>
		<trait name="rate" tag="rate">
			<strictClockCenancestorBranchRates idref="branchRates"/>
		</trait>
		<posterior idref="posterior"/>
	</logTree>
</mcmc>
<report>
	<property name="timer">
		<mcmc idref="mcmc"/>
	</property>
</report>
</beast>
''')
outfile.close()
