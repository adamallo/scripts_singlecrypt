# simulate evolution down a tree of two tips according to
# our SCA model.  Uses horserace logic to provide a quasi
# independent check of the matrix math in the evaluator.
# Considers states up to 6 copies

import random
import math
import dendropy

#null = 0
#a = 7
#b = 1
#aa = 13
#ab = 8
#bb = 2
#aaa = 18
#aab = 14
#abb = 9
#bbb = 3
#aaaa = 22
#aaab = 19
#aabb = 15
#abbb = 10
#bbbb = 4
#aaaaa = 25
#aaaab = 23
#aaabb = 20
#aabbb = 16
#abbbb = 11
#bbbbb = 5
#aaaaaa = 27
#aaaaab = 26
#aaaabb = 24
#aaabbb = 21
#aabbbb = 17
#abbbbb = 12
#bbbbbb = 6

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

def simbranchinseq(endtime,seq):
  answers = []
  for curstate in seq:
    curtime = 0.0
    while curtime < endtime:
      (curtime,curstate) = nextstate(curstate,endtime)
    answers.append(curstate)
  return answers


def simseqtree(numloci,tree):
  root_states=[ab]*numloci
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
      seq_list.append(simbranchinseq(length,par_seq))
    else:
       # no tail node: root
      n_prev_seq = len(seq_list)
      if root_states is not None:
        seq_list.append(root_states)
      else:
        assert n_prev_seq > 0
        n_prev_seq -= 1
    
    

# main program

numloci = 500
#t1 = 0.2
#t2 = 0.3
#branch1 = simbranch(t1,numloci)
#branch2 = simbranch(t2,numloci)

#Sequence simulation using our model along a rooted input newick tree
tree=dendropy.Tree.get(path="tree.trees",rooting="force-rooted",schema="newick")
simseqtree(numloci,tree)

#Matrix which will obtain the data still pending from the tree
char_matrix= dendropy.StandardCharacterMatrix(taxon_namespace=tree.taxon_namespace,default_state_alphabet=None)
char_matrix.taxon_namespace = tree.taxon_namespace

#Object of class with methods to harvest the sequences from the tree
seq_evolver= dendropy.model.discrete.DiscreteCharacterEvolver(seq_attr="sequence")
seq_evolver.extend_char_matrix_with_characters_on_tree(char_matrix=char_matrix,tree=tree)

#Resulting dataset
dataset = dendropy.DataSet()
dataset.add_char_matrix(char_matrix=char_matrix)
dataset.write(path="alignment.out",schema="nexus")

#outfile = open("data","w")
#
#data1 = ""
#for site in branch1:
#  data1 += str(site) + " " 
#data1 += "\n"
#outfile.write(data1)
#
#data2 = ""
#for site in branch2:
#  data2 += str(site) + " " 
#data2 += "\n"
#outfile.write(data2)
#outfile.close()

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
