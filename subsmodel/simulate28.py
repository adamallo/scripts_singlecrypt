# simulate evolution down a tree of two tips according to
# our SCA model.  Uses horserace logic to provide a quasi
# independent check of the matrix math in the evaluator.
# Considers states up to 6 copies

import random
import math
c = 0.1
d = 0.2
g = 0.4

failnum = -1.0
bigtime = 10000000.0

# draw an exponential waiting time
def drawtime(prob):
  if prob == 0.0:
    return failnum
  u = random.random()
  tyme = -math.log(u)/prob
  return tyme

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

  if curstate == 0:   # null; no horserace as outcome is fixed
    return (endtime,curstate)

  elif curstate == 1: # A
    horses.append(Horse(d,0))   # null
    horses.append(Horse(g,3))   # AA

  elif curstate == 2: # B
    horses.append(Horse(d,0))   # null  
    horses.append(Horse(g,5))   # BB

  elif curstate == 3: # AA
    horses.append(Horse(2.0*d,1))  # A
    horses.append(Horse(2.0*g,6))  # AAA

  elif curstate == 4: # AB
    horses.append(Horse(d,1))   # A
    horses.append(Horse(d,2))   # B
    horses.append(Horse(c,3))   # AA
    horses.append(Horse(c,5))   # BB
    horses.append(Horse(g,7))   # AAB
    horses.append(Horse(g,8))   # ABB
      
  elif curstate == 5: # BB
    horses.append(Horse(2.0*d,2))  # B
    horses.append(Horse(2.0*g,9))  # BBB

  elif curstate == 6: # AAA
    horses.append(Horse(3.0*d,3))   # AA
    horses.append(Horse(3.0*g,10))  # AAAA

  elif curstate == 7: # AAB
    horses.append(Horse(d,3))       # AA
    horses.append(Horse(2.0*d,4))   # AB
    horses.append(Horse(c,6))       # AAA
    horses.append(Horse(2.0*0.5*c,8)) # ABB 
    horses.append(Horse(2.0*g,11))  # AAAB
    horses.append(Horse(g,12))      # AABB

  elif curstate == 8: # ABB
    horses.append(Horse(2.0*d,4))   # AB
    horses.append(Horse(d,5))       # BB
    horses.append(Horse(2.0*0.5*c,7))  # AAB
    horses.append(Horse(c,9))       # BBB
    horses.append(Horse(g,12))      # AABB
    horses.append(Horse(2.0*g,13))  # ABBB

  elif curstate == 9: # BBB
    horses.append(Horse(3.0*d,5))   # BB
    horses.append(Horse(3.0*g,14))  # BBBB

  elif curstate == 10: # AAAA
    horses.append(Horse(4.0*d,6))   # AAA
    horses.append(Horse(4.0*g,15))  # AAAAA
    
  elif curstate == 11: # AAAB
    horses.append(Horse(d,6))       # AAA
    horses.append(Horse(3.0*d,7))   # AAB
    horses.append(Horse(c,10))      # AAAA
    horses.append(Horse(3.0*1.0/3.0*c,12))  # AABB
    horses.append(Horse(3.0*g,16))  # AAAAB
    horses.append(Horse(g,17))      # AAABB
    
  elif curstate == 12: # AABB
    horses.append(Horse(2.0*d,7))   # AAB
    horses.append(Horse(2.0*d,8))   # ABB
    horses.append(Horse(2.0*2.0/3.0*c,11))  # AAAB
    horses.append(Horse(2.0*2.0/3.0*c,13))  # ABBB
    horses.append(Horse(2.0*g,17))  # AAABB
    horses.append(Horse(2.0*g,18))  # AABBB
    
  elif curstate == 13: # ABBB
    horses.append(Horse(3.0*d,8))   # ABB
    horses.append(Horse(d,9))       # BBB
    horses.append(Horse(3.0*1.0/3.0*c,12))  # AABB
    horses.append(Horse(c,14))      # BBBB
    horses.append(Horse(g,18))      # AABBB
    horses.append(Horse(3.0*g,19))  # ABBBB

  elif curstate == 14: # BBBB
    horses.append(Horse(4.0*d,9))   # BBB
    horses.append(Horse(4.0*g,27))  # BBBBB

  elif curstate == 15: # AAAAA
    horses.append(Horse(5.0*d,10))  # AAAA
    horses.append(Horse(5.0*g,21))  # AAAAAA

  elif curstate == 16: # AAAAB
    horses.append(Horse(d,10))      # AAAA
    horses.append(Horse(4.0*d,11))  # AAAB
    horses.append(Horse(c,15))      # AAAAA
    horses.append(Horse(4.0*1.0/4.0*c,17))  # AAABB
    horses.append(Horse(4.0*g,22))  # AAAAAB
    horses.append(Horse(g,23))      # AAAABB

  elif curstate == 17: # AAABB
    horses.append(Horse(2.0*d,11))  # AAAB
    horses.append(Horse(3.0*d,12))  # AABB
    horses.append(Horse(2.0*3.0/4.0*c,16))  # AAAAB
    horses.append(Horse(3.0*2.0/4.0*c,18))  # AABBB
    horses.append(Horse(3.0*g,23))  # AAAABB
    horses.append(Horse(2.0*g,24))  # AAABBB

  elif curstate == 18: # AABBB
    horses.append(Horse(3.0*d,12))  # AABB
    horses.append(Horse(2.0*d,13))  # ABBB
    horses.append(Horse(3.0*2.0/4.0*c,17))  # AAABB
    horses.append(Horse(2.0*3.0/4.0*c,19))  # ABBBB
    horses.append(Horse(2.0*g,24))  # AAABBB
    horses.append(Horse(3.0*g,25))  # AABBBB

  elif curstate == 19: # ABBBB
    horses.append(Horse(4.0*d,13))  # ABBB
    horses.append(Horse(d,14))      # BBBB
    horses.append(Horse(4.0*1.0/4.0*c,18))  # AABBB
    horses.append(Horse(c,20))      # BBBBB
    horses.append(Horse(g,25))      # AABBBB
    horses.append(Horse(4.0*g,26))  # ABBBBB

  elif curstate == 20: # BBBBB
    horses.append(Horse(5.0*d,14))  # BBBB
    horses.append(Horse(5.0*g,27))  # BBBBBB

  elif curstate == 21: # AAAAAA
    horses.append(Horse(6.0*d,15))  # AAAAA

  elif curstate == 22: # AAAAAB
    horses.append(Horse(d,15))      # AAAAA
    horses.append(Horse(5*d,16))    # AAAAB
    horses.append(Horse(c,21))      # AAAAAA
    horses.append(Horse(5.0*1.0/5.0*c,23))  # AAAABB

  elif curstate == 23: # AAAABB
    horses.append(Horse(2.0*d,16))  # AAAAB
    horses.append(Horse(4.0*d,17))  # AAABB
    horses.append(Horse(2.0*4.0/5.0*c,22))  # AAAAAB
    horses.append(Horse(4.0*2.0/5.0*c,24))  # AAABBB

  elif curstate == 24: # AAABBB
    horses.append(Horse(3.0*d,17))  # AAABB
    horses.append(Horse(3.0*d,18))  # AABBB
    horses.append(Horse(3.0*3.0/5.0*c,23))  # AAAABB
    horses.append(Horse(3.0*3.0/5.0*c,25))  # AABBBB

  elif curstate == 25: # AABBBB
    horses.append(Horse(4.0*d,18))  # AABBB
    horses.append(Horse(2.0*d,19))  # ABBBB
    horses.append(Horse(4.0*2.0/5.0*c,24))  # AAABBB
    horses.append(Horse(2.0*4.0/5.0*c,26))  # ABBBBB

  elif curstate == 26: # ABBBBB
    horses.append(Horse(5.0*d,19))  # ABBBB
    horses.append(Horse(d,20))      # BBBBB
    horses.append(Horse(5.0*1.0/5.0*c,25))  # AABBBB
    horses.append(Horse(c,27))      # BBBBBB

  elif curstate == 27: # BBBBBB
    horses.append(Horse(6.0*d,20))  # BBBBB

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

# main program

numloci = 500
t1 = 0.2
t2 = 0.3
branch1 = simbranch(t1,numloci)
branch2 = simbranch(t2,numloci)

outfile = open("data","w")

data1 = ""
for site in branch1:
  data1 += str(site) + " " 
data1 += "\n"
outfile.write(data1)

data2 = ""
for site in branch2:
  data2 += str(site) + " " 
data2 += "\n"
outfile.write(data2)
outfile.close()

#
#counts = [0 for x in xrange(15)]
#for i in results:
#  counts[i] += 1
#
#print "**wild type"
#print "AB",counts[4]
#print "**single delete"
#print "A",counts[1]
#print "B",counts[2]
#print "**double delete"
#print "Null",counts[0]
#print "**LOH"
#print "AA",counts[3]
#print "BB",counts[5]
#print "**single gain"
#print "AAB",counts[7]
#print "ABB",counts[8]
#print "**gain with LOH"
#print "AAA",counts[6]
#print "BBB",counts[9]
#print "**double gain"
#print "AAAB",counts[11]
#print "ABBB",counts[13]
#print "**double gain with LOH"
#print "AAAA",counts[10]
#print "BBBB",counts[14]
#print "**balanced gain"
#print "AABB",counts[12]
