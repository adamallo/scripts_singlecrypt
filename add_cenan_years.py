import dendropy
from dendropy import TreeList,Taxon,Node
import sys
import argparse

parser = argparse.ArgumentParser(description="Parses a Newick tree file, modifying the branch lengths from number of generations to years and adding an outgroup")
parser.add_argument("-gt",type=float,default=0,required=False,help="Generation time")
parser.add_argument("-od",type=float,default=0,required=False,help="Outgroup branch length")
parser.add_argument("-i",type=str,default="infile.tree",required=True,help="Input Newick tree file")
parser.add_argument("-o",type=str,default="outtree.tree",required=False,help="Output Newick tree file")
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
trees.write(path=args.o,schema="newick",suppress_rooting=True)		
