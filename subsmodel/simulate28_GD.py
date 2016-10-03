# simulate evolution down a tree of two tips according to
# our SCA model.  Uses horserace logic to provide a quasi
# independent check of the matrix math in the evaluator.
# Considers states up to 6 copies

import random
import math
import dendropy

c = 0.3
d = 0.1
g = 0.2

# probability of genome doubling per unit time
pGD = 0.6

failnum = -1.0
bigtime = 10000000.0

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
def simbranchinseq(brlen,seq,alreadydoubled):
  wasdoubled = alreadydoubled
  if not alreadydoubled:
    gdtime = drawtime(pGD)
  else:
    gdtime = bigtime
  if gdtime < brlen:   # genome doubling happens
    print "Genome Doubling!"
    wasdoubled = True
    answers = simbranchinseq_segment(gdtime,seq)
    for i in xrange(len(answers)):
      answers[i] = gendouble(answers[i])
    answers = simbranchinseq_segment(brlen - gdtime,answers)
  else:
    answers = simbranchinseq_segment(brlen,seq)
  return (answers, wasdoubled)
    
def simbranchinseq_segment(endtime,seq):
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
      assert hasattr(par, "doubled")
      wasdoubled = getattr(par, "doubled")
      (newseq,wasdoubled) = simbranchinseq(length,par_seq,wasdoubled)
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


# main program

numloci = 50
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

# HACK:  attempts to change the taxon names to indicate which were
# genome-doubled encounteed very mysterious behavior from dendropy, so
# instead we write our own NEXUS file.

doublednames = []
for tip in tree.leaf_node_iter():
  if getattr(tip,"doubled"):
    doublednames.append(tip.taxon.label)

infile = open("alignment.out","r")
input = infile.readlines()
infile.close()
(names,seqs,header,footer) = readnexus(input)
outfile = open("alignment.out","w")
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
outline += str(len(seqs[0])) + ";\n"
outfile.write(outline)
outfile.write('\tFORMAT DATATYPE=STANDARD SYMBOLS="";\n')
outfile.write("\tMATRIX\n")
for name,seq in zip(names,seqs):
  if name in doublednames:
    name += "*"
  outline = "\t\t" + name + "\t" + seq + "\n"
  outfile.write(outline)
outfile.write("\t;\n")
outfile.write("END;\n")
outfile.close()
