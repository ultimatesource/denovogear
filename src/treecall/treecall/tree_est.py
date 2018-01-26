#!/usr/bin/env python2

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>
# Author: Rachel Schwartz <Rachel dot Schwartz at asu dot edu>
# Author: Kael Dai <Kael dot Dai at asu dot edu>

from __future__ import print_function
import warnings
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import sys
import itertools
import numpy as np
from scipy.stats import sem
import vcf

from utils import *

with warnings.catch_warnings(ImportWarning):
    from ete2 import Tree

warnings.filterwarnings('error')

DELTA=0.0001  #move this so it's not global

def neighbor_main(args):
    """generate neighbor-joining tree then do recursive NNI and recursive reroot

    Args:
        vcf(str): input vcf/vcf.gz file, "-" for stdin
        output(str): output basename
        mu (int): mutation rate in Phred scale, default 80
            WHY IS THE DEFAULT 80????? IE 10^-8!!!!!!!!!!!!!!!!
        het (int): heterozygous rate in Phred scale, default 30
        min_ev(int): minimum evidence in Phred scale for a site to be considered
            default 60
    
    Output:
        newick trees
    
    """
   
    print(args, file=sys.stderr)
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    #variants =  np.array (tuple): variant info (chrom, pos, ref)  for each variant
    #DPRs = np.array (int): Number of high-quality bases observed for each of the 2 most common alleles for each variant
    #PLs = np.array (int): List of Phred-scaled genotype likelihoods for each of the 2 most common alleles (3 genotypes) for each variant
    
    GTYPE3 = np.array(('RR','RA','AA'))
    base_prior = make_base_prior(args.het, GTYPE3) # base genotype prior; heterozygous rate in Phred scale, default 30; e.g. for het=30 [ 3.0124709,  33.012471,  3.0124709]
    mm,mm0,mm1 = make_mut_matrix_gtype3(args.mu) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)  # pairwise differences between samples based only on PLs (should include mutation, but also shouldn't matter)
    allscores = []
    
    fo = open(args.output+'.scores.txt','w')
    
    for i in range(n_smpl+2):  #10 different starting trees
        print('Tree '+str(i+1)+' of '+str(n_smpl+2))
        tree = init_star_tree(n_smpl)
        internals = np.arange(n_smpl)
        
        #2nd to last tree is nj tree (tho with raw scores not adjusted for saturation)
        if i == n_smpl:
            D,tree = neighbor_joining(D.copy(), tree.copy(), internals) #haven't checked this; make nj tree and update D given internal nodes; pass copy
            
        #last tree is partition tree
        elif i == n_smpl+1:
            tree = Tree()
            #sem calculates the standard error of the mean
            #check if sem for col 1
            if sem(PLs[...,1],axis=1).mean() > sem(PLs[...,2],axis=1).mean():  
                partition(PLs[...,0:2], tree, np.arange(n_smpl), args.min_ev)
            else:
                partition(PLs, tree, np.arange(n_smpl), args.min_ev)
            
        #all other trees are semi-random
        else:
            tree.set_outgroup(str(i))
            tree.resolve_polytomy()
    
        tree = init_tree(tree.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
        tree = populate_tree_PL(tree.copy(), PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
        tree = calc_mut_likelihoods(tree.copy(), mm0, mm1)  #add PLs w mutation
        
#        tree.write(outfile=args.output+'.nj0.tre', format=5)
        
        rerooted = 1
        while rerooted > 0:
            best_tree,best_PL = recursive_NNI(tree.copy(), PLs, mm0, mm1, base_prior,DELTA)
            #print(best_tree)
            best_tree,best_PL,rerooted = recursive_reroot(best_tree.copy(), PLs,mm0, mm1, base_prior,DELTA)  #why are brlens negative?
            #print(best_tree)
            print('PL_per_site = %.4f' % (best_PL/n_site))
            best_tree.write(outfile=args.output+'.'+str(i)+'.tre', format=5)  #write best tree
            #replace sample numbers with actual names
            for node in best_tree.traverse("postorder"):
                if node.is_leaf():
                    node.name=vcffile.samples[int(node.name)]
                
            best_tree.write(outfile=args.output+'.'+str(i)+'names.tre', format=5)  #write best tree
            fo.write(str(i) + ' ' + str(best_PL) + "\n")
            allscores.append(best_PL)
        i+=1
    
    print(allscores)
    
    fo.close
    
def init_star_tree(n):
    """Creates a tree, adds n children in star with numbers as names

    Args:
        n (int): Number of children in tree

    Returns:
        Tree: 
    """
    
    tree = Tree()
    for i in xrange(n):
        tree.add_child(name=str(i))
    return tree


def pairwise_diff(PLs, i, j):
    """
    Calculates difference between pairs of samples - sum across variants - where dif for any given var depends on PL for each possible genotype
    For each var: Dif geno will be close to 1; Same geno will be close to 0
    But depends slightly on PL
    
    Returns:
        numpy.float128 (one value): sum of diffs across vars
    """
    
    pli = normalize2d_PL(PLs[:,i])  #adjust PLs for sample i slightly
    plj = normalize2d_PL(PLs[:,j])  #adjust PLs for sample j slightly
    p = phred2p(pli+plj) # n x g
    return (1-p.sum(axis=1)).sum()  


def make_D(PLs):
    """
    Get pairwise differences between samples based on PLs (e.g. for generating nj tree)
    
    Args:
        PLs (np.array (longdouble)): List of Phred-scaled genotype likelihoods
            for each of the 2 most common alleles for each variant
        
    Returns:
        np.array (longdouble): matrix of n x n where n = number nodes in bifurcating tree
            includes diff btwn tips but not internal nodes - computed during tree building
    """
    n,m,g = PLs.shape   #n_site,n_smpl,n_gtype
    D = np.zeros(shape=(2*m-2,2*m-2), dtype=np.longdouble) #2*m-2 = the total number of nodes of unrooted tree (incl leaves and internal)
    for i,j in itertools.combinations(xrange(m),2):
        D[i,j] = pairwise_diff(PLs, i, j)
        D[j,i] = D[i,j]
    return D


def neighbor_joining(D, tree, internals):
    #fsum will have better precision when adding distances across sites
    #based on PLs not mutation
    """
    
    Args:
        D (np.array): pairwise differences between samples based on PLs (passing copy)
        tree (Tree): tree of class Tree with num tips = num samples
        internals (np.array): array of sample numbers
        
    Returns:
        Tree
        D (np.array): update pairwise differences now there are internal nodes to compare
    
    """
    print('neighbor_joining() begin', end=' ', file=sys.stderr)
    m = len(internals)
    while m > 2:  #if m is 2 then only two connected to root
        d = D[internals[:,None],internals]  #initially D matrix w/o 0 distance btwn internal nodes; then add in nodes as they have distances
        u = d.sum(axis=1)/(m-2)

        Q = np.zeros(shape=(m,m), dtype=np.longdouble)
        for i,j in itertools.combinations(xrange(m),2):  #std Q matrix calc
            Q[i,j] = d[i,j]-u[i]-u[j]
            Q[j,i] = Q[i,j]
        #print(Q.astype(int))
        np.fill_diagonal(Q, np.inf)
        #print(np.unique(Q, return_counts=True))
        i,j = np.unravel_index(Q.argmin(), (m,m))  #location in matrix of smallest Q value (ie closest nodes/tips)
        l = len(D)+2-m

        for k in xrange(m):
            D[l,internals[k]] = D[internals[k],l] = d[i,k]+d[j,k]-d[i,j]
        D[l,internals[i]] = D[internals[i],l] = vi = (d[i,j]+u[i]-u[j])/2
        D[l,internals[j]] = D[internals[j],l] = vj = (d[i,j]+u[j]-u[i])/2

        ci = tree&str(internals[i])
        cj = tree&str(internals[j])
        ci.detach()
        cj.detach()
        node = Tree(name=str(l))
        node.add_child(ci,dist=int(vi))
        node.add_child(cj,dist=int(vj))
        tree.add_child(node)
        #print(tree)

        internals = np.delete(internals, [i,j])
        internals = np.append(internals, l)
        m = len(internals)
        print('.', end='', file=sys.stderr)

    print(' done', file=sys.stderr)
    return D,tree

def update_PL(node, mm0, mm1):
    """
    PL for nodes depend on children so must be updated if node children change due to nni/reroot
    
    node has
        PL0: PLs given no mutation
        PLm: w mutation
        children
        
    changes PLm and PL0 for node and all its children (recursive); also sid for children
    
    returns:
        Tree: could be subtree; includes all recursively updated PL0 and PLm and node labels (sid)
    """

    n,g = node.PL0.shape
    l = 2*len(node)-2
    #node.PL0 = np.zeros((n,g), dtype=np.longdouble)
    node.PL0.fill(0.0)
    node.PLm = np.zeros((l,n,g), dtype=np.longdouble)
    for child in node.children:
        sid = sorted(map(int,child.get_leaf_names()))
        if child.sid != sid:  #sid is supposed to be names of leaves - could be dif due to swapping in nni
            newchild = update_PL(child, mm0, mm1)
            newchild.sid = sid
            child.detach()
            node.add_child(newchild)
        node.PL0 += p2phred(np.dot(phred2p(child.PL0), mm0)) 
    i = 0
    for child in node.children:
        sister = child.get_sisters()[0]
        if not child.is_leaf():
            l = child.PLm.shape[0]
            node.PLm[i:(i+l)] = p2phred(np.dot(phred2p(child.PLm), mm0)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
            i += l
        node.PLm[i] = p2phred(np.dot(phred2p(child.PL0), mm1)) + p2phred(np.dot(phred2p(sister.PL0), mm0)) 
        i += 1

    return node

def score(tree, base_prior):
    """
    used to compare trees
    
    Args:
        Tree
            PL0: PLs for all vars for leaf or dot product of child's matrix and mutation matrix
                where mm0: substitution rate matrix with non-diagonal set to 0
            PLm:
                    mm1: substitution rate matrix with diagonal set to 0
    
    """
    Pm = phred2p(tree.PLm+base_prior).sum(axis=(0,2))       #why add baseprior
    P0 = phred2p(tree.PL0+base_prior).sum(axis=1)
    return p2phred(Pm+P0).sum()


def partition(PLs, tree, sidx, min_ev):
    """requires
            make_selection_matrix2
            calc_minimum_pt_cost
            
    """
    
    if tree.is_root():
        print('partition() begin', end=' ', file=sys.stderr)
    m = len(sidx) # number of samples under current node
    print(m, end='.', file=sys.stderr)
    if m == 2:
        child1 = tree.add_child(name=str(sidx[0]))
        child1.add_features(samples=np.atleast_1d(sidx[0]))
        child2 = tree.add_child(name=str(sidx[1]))
        child2.add_features(samples=np.atleast_1d(sidx[1]))
    elif m > 2:
        smat = make_selection_matrix2(m)
        pt, cost = calc_minimum_pt_cost(PLs, smat, min_ev)
        k0 = pt==0
        sidx0 = np.atleast_1d(sidx[k0])
        child = tree.add_child(name=','.join(sidx0.astype(str)))
        child.add_features(samples=sidx0)
        if len(sidx0) > 1:
            partition(PLs[:,k0,], child, sidx0, min_ev)
        k1 = pt==1
        sidx1 = np.atleast_1d(sidx[k1])
        child = tree.add_child(name=','.join(sidx1.astype(str)))
        child.add_features(samples=sidx1)
        if len(sidx1) > 1:
            partition(PLs[:,k1,], child, sidx1, min_ev)
    else:
        print('m<=1: shouldn\'t reach here', file=sys.stderr)
        sys.exit(1)
    if tree.is_root():
        print(' done', file=sys.stderr)


def calc_minimum_pt_cost(PLs, smat, min_ev):
    n,m,g = PLs.shape
    pt_cost = np.inf
    for k in smat:
        x0 = PLs[:,k==0,].sum(axis=1) # dim = n_site x 2
        x0min = x0.min(axis=1) # dim = n_site x 1
        x0max = x0.max(axis=1) # dim = n_site x 1
        x1 = PLs[:,k==1,].sum(axis=1) # dim = n_site x 2
        x1min = x1.min(axis=1) # dim = n_site x 1
        x1max = x1.max(axis=1) # dim = n_site x 1
        # take everything
        #c = (x0 + x1).sum()
        # cap the penalty by mu
        #c = (x0>mu).sum()*mu + x0[x0<=mu].sum() + (x1>mu).sum()*mu + x1[x1<=mu].sum()
        # ignore sites where signal from either partitions is weak
        #c = (x0min+x1min)[(x0max>min_ev) & (x1max>min_ev)].sum()
        # ignore sites where signals from both partitions are weak
        c = (x0min+x1min)[(x0max>min_ev) | (x1max>min_ev)].sum()
        # some weird cost function that broadly penalize partition of similar samples
        #k0 = x0.argmin(axis=1)
        #k1 = x1.argmin(axis=1)
        #c = np.minimum(x0[k0],x1[k1]).sum() + (k0==k1).sum()*mu
        if c < pt_cost:
            pt_cost = c
            pt = k
    return pt, pt_cost


def make_selection_matrix(m, t=20):
    n = 2**(m-1)
    if m>3 and m<=t: # special treatment for intermediate size
        l = (m,)*n
        x = np.array(map(tuple, map(str.zfill, [b[2:] for b in map(bin, xrange(4))], (3,)*4)), dtype=np.byte)
        y = np.zeros((n,m),dtype=np.byte)
        for i in xrange(m-3):
            a,b = x.shape
            y[0:a,-b:] = x
            y[a:(2*a),-b:] = x
            y[a:(2*a),-b] = 1
            x = y[0:(2*a),-(b+1):]
        for s in y:
            yield s
    else:
        for i in xrange(n):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)


def make_selection_matrix2(m, t=20):
    n = 2**(m-1)
    if m>3 and m<=t: # special treatment for intermediate size
        l = (m,)*n
        x = np.array(map(tuple, map(str.zfill, [b[2:] for b in map(bin, xrange(4))], (3,)*4)), dtype=np.byte)
        y = np.zeros((n,m),dtype=np.byte)
        for i in xrange(m-3):
            a,b = x.shape
            y[0:a,-b:] = x
            y[a:(2*a),-b:] = x
            y[a:(2*a),-b] = 1
            x = y[0:(2*a),-(b+1):]
        for s in y:
            yield s
    elif m<=3:
        for i in xrange(n):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)
    else:
        r1 = np.random.randint(1,m-1,2**t)
        r2 = np.random.rand(2**t)
        x = ((1+r2)*2**r1).astype(int)
        for i in iter(x):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)


def reroot(tree, PLs, mm0, mm1, base_prior,DELTA):
    """
    
    return:
        Tree
        np.array (PLs)
        int: flag if rerooted (1) or not (0)
    """
    '''
              /-A              /-A              /-B
           /-|              /-|              /-|
    -root-|   \-B => -root-|   \-C => -root-|   \-C
          |                |                |
           \-C              \-B              \-A
    '''

    best_tree = tree
    best_PL = score(tree, base_prior)
    flag = 0

    for node in tree.iter_descendants('postorder'):  #go through all nodes including tips but not root
        tree_reroot = tree.copy()
        new_root = tree_reroot.search_nodes(sid=node.sid)[0]  #gets node of interest
        tree_reroot.set_outgroup(new_root)  #sets node of interest to outgroup
        
        tree_reroot = init_tree(tree_reroot)  #tree has nid's (node id) and sid's (list of tip names - sorted)
        tree_reroot = populate_tree_PL(tree_reroot, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
        tree_reroot = calc_mut_likelihoods(tree_reroot, mm0, mm1)  #add PLs w mutation
        
        #tree_reroot = update_PL(tree_reroot, mm0, mm1)  #new PL given decendants
        PL_reroot = score(tree_reroot, base_prior) 
        #print(tree_reroot)
        #print(PL_reroot)
        if PL_reroot < best_PL * (1-DELTA): #new best tree only if significantly better ie trees could be similar but status quo wins
            best_tree = tree_reroot
            best_PL = PL_reroot
            flag = 1
            
    return best_tree,best_PL,flag


def recursive_reroot(tree, PLs,mm0, mm1, base_prior,DELTA):
    """
    starting at tips, work up tree, get best way of rooting subtree 
    """

    print('recursive_reroot() begin', file=sys.stderr)
    PL = score(tree, base_prior)
    for node in tree.iter_descendants('postorder'):  #go through all nodes including tips but not root
        rerooted = 0
        new_tree = tree.copy()
        node_leaves = node.get_leaf_names()
        new_node = new_tree.get_common_ancestor(node_leaves)  #get corresponding node in new tree
        try:
            #print('Node:')
            #print(node)
            new_tree.set_outgroup(new_node) #reroot
    
            #recalculate
            new_tree = init_tree(new_tree.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
            new_tree = populate_tree_PL(new_tree.copy(), PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
            new_tree = calc_mut_likelihoods(new_tree.copy(), mm0, mm1)  #add PLs w mutation
                   
            PL_new = score(new_tree, base_prior)
            
            #print('This tree:')
            #print(new_tree)
            #print('Has this score:')
            #print(PL_new)
    
            if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                best_tree = new_tree.copy()
                PL = PL_new
                rerooted = 1
        except:
            #print('Can\'t set node as outgroup:')
            #print(node)
            pass
            
    if rerooted == 1:  #there was a better tree
        tree = best_tree.copy()
        #print('Best tree:')
        #print(tree)
    #else:                        
        #print('No change to tree:')
        #print(tree)
        #print(PL)
        
    print(' done', end='', file=sys.stderr)
    print(tree)
    print(PL)
    return tree,PL,rerooted


def nearest_neighbor_interchange(node, PLs,mm0, mm1, base_prior,DELTA):
    '''
    Args:
    
    Return:
        node: Tree (could be a subtree of the full tree)
        np.array(PL)
        int: flag to indicate nni happened
        
    Process:
    
              /-A              /-A              /-A
           /-|              /-|              /-|
          |   \-B          |   \-C          |   \-D
    -node-|       => -node-|       => -node-|
          |   /-C          |   /-B          |   /-B
           \-|              \-|              \-|
              \-D              \-D              \-C
           ||               ||               ||
           \/               \/               \/
        reroot()         reroot()         reroot()
    '''
    
    c1,c2 = node.children  #children of root node
    possible_rearrangements = []
    
    #children are leaves - don't need to swap anything
    if c1.is_leaf() and c2.is_leaf():
        return [None]
    
    #one child is a leaf - rerooting will provide all possible combinations - flagged if rerooted
    elif c1.is_leaf():
        c21,c22 = c2.children
        node.set_outgroup(c22)
        possible_rearrangements.append(node.copy())
        node.set_outgroup(c21)
        possible_rearrangements.append(node.copy())
        return possible_rearrangements
        
    elif c2.is_leaf():
        c12,c11 = c1.children
        node.set_outgroup(c12)
        possible_rearrangements.append(node.copy())
        node.set_outgroup(c11)
        possible_rearrangements.append(node.copy())
        return possible_rearrangements

    else:
        #rerootings of original tree
        node_copy1 = node.copy()
        c1,c2 = node_copy1.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        for n in [c11,c12,c21,c22]:
            node_copy1.set_outgroup(n)
            possible_rearrangements.append(node_copy1.copy())

        #2nd tree - swap relationships and reroot
        node_copy2 = node.copy()
        c1,c2 = node_copy2.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        c12 = c12.detach()
        c22 = c22.detach()
        c1.add_child(c22)
        c2.add_child(c12)
        
        possible_rearrangements.append(node_copy2.copy())
        for n in [c11,c12,c21,c22]:
            node_copy2.set_outgroup(n)
            possible_rearrangements.append(node_copy2.copy())
            
        #3rd tree - swap relationships and reroot
        node_copy3 = node.copy()
        c1,c2 = node_copy3.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        c12 = c12.detach()
        c21 = c21.detach()
        c1.add_child(c21)
        c2.add_child(c12)
        
        possible_rearrangements.append(node_copy3.copy())
        for n in [c11,c12,c21,c22]:
            node_copy3.set_outgroup(n)
            possible_rearrangements.append(node_copy3.copy())

        return possible_rearrangements


def recursive_NNI(tree, PLs, mm0, mm1, base_prior,DELTA):
    #recursive just means traverse the tree 
    """
    
    Args:
        tree(Tree)
        mm0: mutation matrix (np array of float) (non-diagonal set to 0)
        mm1: mutation matrix (np array of float) (diagonal set to 0)
        base_prior (np.array): Base prior probs depending on het pl

    Returns:
        Tree (tree)
        np.array (PL): phred-scaled likelihooods
        
    tree resulting from each round of nni (working from tips to root) printed to trees_tried.txt
    
    """
    print('recursive_NNI() begin', end=' ', file=sys.stderr)
    PL = score(tree, base_prior)
    #goes until can get through tree w/o nni at any node
    #a la phylip
    num_nnis=1
    #print(tree)
    #print(PL)
    while(num_nnis>0):
        num_nnis=0
        print('Start nni round')
#        for node in tree.traverse('postorder'):
#            print(node)
        for node in tree.traverse('postorder'):
            #goes through each node, does nni if better
            if node.is_leaf():
                continue
            print('.', end='', file=sys.stderr)
            #print(node)
            possible_rearrangements = nearest_neighbor_interchange(node.copy(),PLs, mm0, mm1, base_prior,DELTA)
#            print('Original original tree:')
#            print(tree)
            if possible_rearrangements[0] is not None:
#                for r in possible_rearrangements:
#                    print(r)
                if node.is_root():  #because can't get parent as for below
                    for r in possible_rearrangements:
                        #print('Rearranged node:')
                        #print(r)
                        
                        new_tree = init_tree(r.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
                        new_tree = populate_tree_PL(new_tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
                        new_tree = calc_mut_likelihoods(new_tree, mm0, mm1)  #add PLs w mutation

                        PL_new = score(new_tree, base_prior)
                        
                        #print('This tree:')
                        #print(new_tree)
                        #print('Has this score:')
                        #print(PL_new)

                        if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                            best_tree = new_tree.copy()
                            #print(best_tree)
                            PL = PL_new
                            num_nnis = 1
                    if num_nnis == 1:  #there was a better tree
                        tree = best_tree.copy()
                        #print('Best tree so far:')
                        #print(tree)
                    #else:                        
                     #   print('No change to tree:')
                     #   print(tree)
                     #   print(PL)

                else:
                    for r in possible_rearrangements:
                        #print('Rearranged node:')
                        #print(r)
                        new_tree = tree.copy()
                        node_leaves = node.get_leaf_names()
                        new_node = new_tree.get_common_ancestor(node_leaves)  #get corresponding node in new tree
                        
                        parent = new_node.up
                        new_node.detach()
                        parent.add_child(r)
                        
                        #modify copy of tree
                        new_tree = init_tree(new_tree)  #tree has nid's (node id) and sid's (list of tip names - sorted)
                        new_tree = populate_tree_PL(new_tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
                        new_tree = calc_mut_likelihoods(new_tree, mm0, mm1)  #add PLs w mutation
                        
                        PL_new = score(new_tree, base_prior)
                                                
                        #print('This tree:')
                        #print(new_tree)
                        #print('Has this score:')
                        #print(PL_new)
                        
                        if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                            PL = PL_new
                            num_nnis = 1
                            best_tree = new_tree.copy()
#                        print('Original tree:')
#                        print(tree)

                    if num_nnis == 1:  #there was a better tree
                        tree = best_tree.copy()
                        #print('Best tree so far:')
                        #print(tree)
                        break  #take best tree and start over because now nni's will be all different
                    #else:                        
                        #print('No change to tree:')
                        #print(tree)
                        #print(PL)
                    
        #print(str(num_nnis)+' nnis', end='', file=sys.stderr)
        #print(PL)
        
    print(' done', file=sys.stderr)
    print(tree)
    print(PL)
    return tree,PL_new