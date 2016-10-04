import dendropy
from dendropy import TreeList,Taxon,Node
import sys
import argparse

xml_sufix=".xml"
parser = argparse.ArgumentParser(description="Parses a Newick tree file, modifying the branch lengths from number of generations to years and adding an outgroup")
parser.add_argument("-gt",type=float,default=0,required=False,help="Generation time")
parser.add_argument("-od",type=float,default=0,required=False,help="Outgroup branch length")
parser.add_argument("-i",type=str,default="infile.tree",required=True,help="Input Newick tree file")
parser.add_argument("-o",type=str,default="outtree.tree",required=False,help="Output Newick tree file")
parser.add_argument("-ox",type=str,default="tree",required=True,help="Output truncated XML file for BEAST")
parser.add_argument("-a",type=float,default=40,required=True,help="Time from the DOB to the Cenancestor")
args = parser.parse_args()

trees=TreeList.get_from_path(args.i,schema="newick",rooting="force-rooted")
if args.gt != 0:
	print "Scaling branch lengths to time with generation time %d\n" % args.gt
	for tree in trees:
		for edge in tree.preorder_edge_iter():
			#print "DEBUG: %s" % edge.length
			if edge.length != None:
				edge.length=edge.length/args.gt

if args.od != 0:
	print "Adding outgroup with branch length %d\n" % args.od
	namespace=trees.taxon_namespace
	outgroup= Taxon("outgroup")
	namespace.add_taxon(outgroup)
	ntree=0
	labels=namespace.labels()
	labels.remove("outgroup")
	for tree in trees:
		outgroup_node=Node(taxon=outgroup,edge_length=args.od)
		new_root_node=Node()
		tree.seed_node.edge_length=args.od-tree.seed_node.distance_from_tip()
		new_root_node.add_child(tree.seed_node)
		new_root_node.add_child(outgroup_node)
		tree.seed_node=new_root_node	
		##Truncated XML output
		tree.calc_node_root_distances()
		outfile = open(args.ox + str(ntree) + xml_sufix,"w") #Truncated XML
		outfile.write("<?xml version=\"1.0\" standalone=\"yes\"?>\n<beast>\n<taxa id=\"taxa\">\n")
		#Get labels and discard the outgroup label
		for label in labels:
			date= args.a + tree.find_node_for_taxon(namespace.get_taxon(label)).distance_from_root()#Date equals the distance from the root plus a constant
			outfile.write("\t<taxon id=\"" + label + "\">\n\t\t<date value=\"" + str(round(date,0)) + "\" direction=\"forwards\" units=\"years\"/>\n\t</taxon>\n")
		outfile.write("<taxa/>")
		ntree+=1
		outfile.write('\n<generalDataType id="cnv">\n\t<state code="@"/> <!-- Genotype: 0,0 ; Beast State: 0 -->\n\t<state code="A"/> <!-- Genotype: 0,1 ; Beast State: 1 -->\n\t<state code="B"/> <!-- Genotype: 0,2 ; Beast State: 2 -->\n\t<state code="C"/> <!-- Genotype: 0,3 ; Beast State: 3 -->\n\t<state code="D"/> <!-- Genotype: 0,4 ; Beast State: 4 -->\n\t<state code="E"/> <!-- Genotype: 0,5 ; Beast State: 5 -->\n\t<state code="F"/> <!-- Genotype: 0,6 ; Beast State: 6 -->\n\t<state code="G"/> <!-- Genotype: 1,0 ; Beast State: 7 -->\n\t<state code="H"/> <!-- Genotype: 1,1 ; Beast State: 8 -->\n\t<state code="I"/> <!-- Genotype: 1,2 ; Beast State: 9 -->\n\t<state code="J"/> <!-- Genotype: 1,3 ; Beast State: 10 -->\n\t<state code="K"/> <!-- Genotype: 1,4 ; Beast State: 11 -->\n\t<state code="L"/> <!-- Genotype: 1,5 ; Beast State: 12 -->\n\t<state code="M"/> <!-- Genotype: 2,0 ; Beast State: 13 -->\n\t<state code="N"/> <!-- Genotype: 2,1 ; Beast State: 14 -->\n\t<state code="O"/> <!-- Genotype: 2,2 ; Beast State: 15 -->\n\t<state code="P"/> <!-- Genotype: 2,3 ; Beast State: 16 -->\n\t<state code="Q"/> <!-- Genotype: 2,4 ; Beast State: 17 -->\n\t<state code="R"/> <!-- Genotype: 3,0 ; Beast State: 18 -->\n\t<state code="S"/> <!-- Genotype: 3,1 ; Beast State: 19 -->\n\t<state code="T"/> <!-- Genotype: 3,2 ; Beast State: 20 -->\n\t<state code="U"/> <!-- Genotype: 3,3 ; Beast State: 21 -->\n\t<state code="V"/> <!-- Genotype: 4,0 ; Beast State: 22 -->\n\t<state code="W"/> <!-- Genotype: 4,1 ; Beast State: 23 -->\n\t<state code="X"/> <!-- Genotype: 4,2 ; Beast State: 24 -->\n\t<state code="Y"/> <!-- Genotype: 5,0 ; Beast State: 25 -->\n\t<state code="Z"/> <!-- Genotype: 5,1 ; Beast State: 26 -->\n\t<state code="["/> <!-- Genotype: 6,0 ; Beast State: 27 -->\n\t<ambiguity code="-" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>\n\t<ambiguity code="?" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>\n</generalDataType\>\n')
trees.write(path=args.o,schema="newick",suppress_rooting=True)
		
