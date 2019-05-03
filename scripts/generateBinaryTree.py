import random
import numpy as np
import sys

# set seed for reproducibilaty of random tree
configSeed = sys.argv[2]
random.seed(configSeed)

# class Node, which has a value, a left and a right child (represented by their values)
# if value = -1, then no Node (child) exists
class Node:
    counter = 0
    child = 0

    def __init__(self, val=(-1), right=(-1), left=(-1)):
        self._val = val
        self._right = right
        self._left = left

# generates a random tree with n leaves (2*n -1 total nodes) and root with value=0
# returns the generated tree as a 2dim vector where at position [x][0] the value of the node x's
# right child and at [x][1] the value of node x's left child stands
def generateTree(n):
#    tree = np.ones([2, (2*n-1)])
    tree = np.ones((3, (2*n-1)))
    tree = tree * (-1)

    root = Node(0)
    rootOfTree = arbitraryTree(tree, root, 2*n-1)

    return tree

# recursivaly generates the tree
def arbitraryTree(tree, root, n):

    if n > 3:

        leftTreeSize = random.randrange(1, n-2, 2)
        rightTreeSize = (n-1) - leftTreeSize

        rootRight = Node()
        rootLeft = Node()

        root._right = arbitraryTree(tree, rootRight, rightTreeSize)
        root._left = arbitraryTree(tree, rootLeft, leftTreeSize)

        Node.counter = Node.counter+1

        # we don't need to do it for the root, because we've already given it the value 0
        if (Node.counter < (tree.size/3)):
            root._val = Node.counter

        tree[0][root._val] = root._right._val
        tree[1][root._val] = root._left._val

        return root

    elif n == 3:
        rootRight = Node()
        rootLeft = Node()

        root._right = arbitraryTree(tree, rootRight, 1)
        root._left = arbitraryTree(tree, rootLeft, 1)

        Node.counter = Node.counter+1
        root._val = Node.counter

        tree[0][root._val] = root._right._val
        tree[1][root._val] = root._left._val

        return root

    else:
        Node.counter = Node.counter+1
        tree[2][Node.counter] = Node.child
        Node.child = Node.child + 1
        root._val = Node.counter
        return root

# using the number of Leaves specified in the config file
N = int(sys.argv[1])

# call to function generateTree with N leaves
testTree = generateTree(N)

# saving tree to the specified output file
np.savetxt(sys.argv[2], testTree, fmt='%i')
