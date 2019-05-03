from graphviz import Digraph
import numpy as np

tree = np.loadtxt(fname = "binaryTree10.txt")

dot = Digraph(comment='Tree')

for i in range(0, int(tree.size/3)):
    dot.node(str(i))

for i in range(0, int(tree.size/3)):
    if tree[0][i] != (-1):
        dot.edge(str(i), str(int(tree[0][i])))
    if tree[1][i] != (-1):
        dot.edge(str(i), str(int(tree[1][i])))

dot.render('Graph.pdf', view=True)
