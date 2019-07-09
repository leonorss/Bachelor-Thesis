from graphviz import Digraph
import numpy as np
import os

# load the tree in matrix format build by generateBinaryTree.py
tree = np.loadtxt(fname = snakemake.input[0], ndmin = 2)

# start a new graph
dot = Digraph(comment='Tree')

# create a node for every column (cell) and give it either the name -1 if it is not a child
# or the given child name, which will be used from this point forward
for i in range(0, int(tree.size/3)):
    if(tree[2][i] != (-1)):
        nodeName = str(int(tree[2][i])) + ", " + str(int(tree[2][i]) + int(snakemake.params[0]))
    else:
        nodeName = str(int(tree[2][i]))
    dot.node(str(i), nodeName)

# create the edges between the nodes, according the loaded tree
for i in range(0, int(tree.size/3)):
    if tree[0][i] != (-1):
        dot.edge(str(i), str(int(tree[0][i])))
    if tree[1][i] != (-1):
        dot.edge(str(i), str(int(tree[1][i])))

name = "Data/Trees/Graph" + snakemake.wildcards.treename

# render the graph
dot.render(name, view=False)

# delete unnecessary files created by the tree rendering
shellCommand = "rm " + name
os.system(shellCommand)
