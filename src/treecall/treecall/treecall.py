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
#from editdistance import eval as strdist
import vcf

with warnings.catch_warnings(ImportWarning):
    from ete2 import Tree

warnings.filterwarnings('error')

from tree_est import *
from utils import *
from geno import *

def compat_main(args):
    """calculate pairwise compatibility between all pairs of sites

    Args:
        args.vcf (str): input vcf/vcf.gz file or - stdin
        args.output (str): file to output compatibility matrix
        args.min_ev (int): minimum evidence in Phred scale for a site to be considered, default 60

    Output to file:
        np.array: compatibility matrix

    """
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    #n_site, n_smpl = PLs.shape[0:2]
    #sidx = np.arange(n_smpl)

    compats = calc_compat(PLs) #matrix of num_var x num_var containing 0 for compatible sites (ie same pattern) or 1 if not
    c = compats.sum(axis=-1)  #not used?
    k = (c == 0) & ~find_singleton(PLs)  #not used?

    gzout = args.output + '.gz'
    np.savetxt(gzout, compats, fmt='%d', delimiter='\t')


def calc_compat(PLs):
    """Create a pairwise compatibility matrix (numpy array) for variants
    
    0 if sites are compatible, 1 if not

    Args:
        PLs (np.array (int)): List of Phred-scaled genotype likelihoods for each of the 2 most common alleles for each variant

    Returns:
        np.array: matrix of num_var x num_var containing 0 for compatible sites (ie same pattern) or 1 if not

    """
    print('calc_compat() begin', end=' ', file=sys.stderr)
    n,m,g = PLs.shape       #get array dimensions - ie n=num_variants, m=num samples, g=num genotypes
    nidx = np.arange(n)     #ints from 0 to num var
    midx = np.arange(m)     #ints from 0 to num samples
    kn = np.tile(nidx,m).reshape(m,n)   #repeat nidx m times
    km = np.repeat(midx,n).reshape(m,n)

    PLs = PLs.astype(int)
    gt = (PLs[...,0]==0).astype(np.byte) # num var x num samples; genotype as 1 or 0 (is or not homozygous ref)
    #WHY ARE NON ZEROS ONE PL V THE OTHER????
    non_zeros = PLs[kn,km,gt.T].T # n x m; cell is PL of homo ref if it's a het or PL of het if homo ref
    #still not sure what below is doing to get matrix or whether it's correct
    groups = (2*gt[i]+gt for i in xrange(n)) # each n x m
    cost = (np.minimum(non_zeros[i], non_zeros) for i in xrange(n)) # each n x m

    compats = np.zeros(shape=(n,n), dtype=np.int32)  #0 array of num_var x num_var to compare each var 
    for i in xrange(n):
        grp = groups.next()
        cst = cost.next()
        compats[i,i:] = map(min, map(np.bincount, grp[i:], cst[i:]))  #only need to fill in half of matrix (symmetrical)
        #
    compats = compats + compats.T - np.diag(compats.diagonal())  #make symmetrical

    print(' done', file=sys.stderr)
    return compats


def find_singleton(PLs):
    #check if there is only one sample w/o homozygous ref genotype 
    is_singleton = (PLs[...,0]>0).sum(axis=1)==1  #get PL for 0/0 for each sample, check if each >0, count these, check if only one
    return is_singleton

# TODO
def split_main(args):
    print(args, file=sys.stderr)
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    n_site, n_smpl = PLs.shape[0:2]
    maf = (DPRs[...,1]>0).sum(axis=1).astype(np.float)/((DPRs.sum(axis=2)>0).sum(axis=1))
    maf[maf>0.5] = 1 - maf[maf>0.5]
    maf_order = maf.argsort()


# TODO
def rsplit_main(args):
    print(args, file=sys.stderr)
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    n,m,g = PLs.shape

    maf = (DPRs[...,1]>0).sum(axis=1).astype(np.float)/((DPRs.sum(axis=2)>0).sum(axis=1))
    maf[maf>0.5] = 1 - maf[maf>0.5]
    k = np.random.choice(n, 100)
    maf_order = maf[k].argsort()[::-1]

    tree = Tree()
    subdiv(PLs[k][maf_order], tree)


# TODO
def subdiv(PLs, tree):
    n,m,g = PLs.shape
    tree.sid = range(m)
    for i in xrange(n):
        for leaf in tree.get_leaves():
            PL = PLs[i,leaf.sid,0]
            k0 = PL==0
            k1 = ~k0
            p0 = PL[k0].sum()
            p1 = PL[k1].sum()

            c0 = Tree()
            c0.sid = leaf.sid[k0]
            leaf.add_child(c0)
            c1 = Tree()
            c1.sid = leaf.sid[k1]
            leaf.add_child(c1)


def annotate_nodes(tree, attr, values):
    #fix this so it returns something and doesn't try to use a global variable
    for node in tree.iter_descendants('postorder'):
        setattr(node, attr, values[node.nid])


def tview_main(args):
    """Takes input newick tree, prints tree with labels

    Args:
        args.tree (str): string containing newick to be converted to Tree
        args.attrs (str): node attributes given by a comma separated list
        args.label (str): leaves label
        
    Prints tree

    """
    tree = Tree(args.tree)
    if args.attrs:
        attrs = args.attrs.split(',')
        if 'label' in attrs and args.label:
            label = read_label(args.label)
            for leaf in tree.iter_leaves():
                leaf.add_feature('label', label.get(leaf.name))
        print(tree.get_ascii(attributes=attrs, show_internal=False))
    else:
        print(tree)


def read_label(filename):
    """from tab delim file: dict of
        key: first col in file / index
        value: second col in file (or first col if only one)
    """
    label = {}
    with open(filename) as f:
        i = 0
        for line in f:
            c = line.rstrip().split('\t')
            if len(c) > 1:
                label[c[0]] = c[1]
            else:
                label[str(i)] = c[0]
            i += 1
    return label


def annotate_main(args):
    print(args, file=sys.stderr)

    tree = Tree(args.tree)
    tree = init_tree(tree)

    gtcall = read_gtcall(args.gtcall)
    for node in tree.iter_descendants('postorder'):
        k = gtcall['mut_smpl'] == ','.join(map(str,node.sid))
        node.dist = k.sum()+1
    tree.write(outfile=args.output, format=5)


def read_gtcall(filename):
    dtype=[('chrom','a10'),('pos','i4'),('ref','a1'),
           ('null_p','f8'),('mut_p','f8'),
           ('null_base','a2'),('null_base_p','f8'),
           ('mut_base','a2'),('mut_alt','a2'),('mut_conf_p','f8'),
           ('mut_loc','i4'),('mut_smpl','a128')]
    if filename == '-':
        filename = sys.stdin
    return np.loadtxt(filename, dtype=dtype)

def compare_main(args):
    	
    """compare tree topologies

    Args:
        args.tree (str): input tree(s), in Newick format
        args.ref (str): reference tree, in Newick format
        
    Prints:
        tree
        result['norm_rf']: normalized robinson-foulds distance (from 0 to 1)
        result['ref_edges_in_source']: compatibility score of the target tree with respect to the source tree (how many edges in reference are found in the source)
        result['source_edges_in_ref']: compatibility score of the source tree with respect to the reference tree (how many edges in source are found in the reference)
        dstat: sum of differences between two distance matrices / sum of ref matrix
        rstat: avg ratio between corresponding pairwise distances

    """
    
    print(args, file=sys.stderr)
    ref_tree = Tree(args.ref)
    ref_tree_leafnames = [l.name for l in ref_tree.get_leaves()]
    leaf_idx = {l:i for i,l in enumerate(ref_tree_leafnames)}  #how to get int for leaf name consistent btwn trees
    
    ref_am = tree2adjacency(ref_tree,leaf_idx)   #matrix of "distances" for ref (node counts)
    for f in args.tree:
        tree = Tree(f)
        tree_leafnames = [l.name for l in tree.get_leaves()]
        if set(tree_leafnames) != set(ref_tree_leafnames):
            print('leaf names are not the same', file=sys.stderr)
        am = tree2adjacency(tree,leaf_idx)   #matrix of "distances" for comparison
        if ref_am.shape != am.shape:
            print('%s incompatible with %s' % (f, args.ref), file=sys.stderr)
        else:
            k = ref_am > 0

            diff = np.abs(ref_am - am)
            dstat = diff[k].sum()/k.sum()

            ratio = am[k]/ref_am[k]
            ratio[ratio>1] = 1.0/ratio[ratio>1]
            rstat = np.power(ratio.prod(), 1.0/k.sum())

            result = ref_tree.compare(tree, unrooted=True)  #comparison calculated by ete2

            # <tree>,<norm_rf>,<ref_edge_in_tree>,<tree_edge_in_ref>,<diff_adj>,<ratio_adj>
            print('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (f, result['norm_rf'], result['ref_edges_in_source'], result['source_edges_in_ref'], dstat, rstat))


def tree2adjacency(tree,leaf_idx):
    """Get the number of nodes between leaves formatted as np array

    Args:
        tree (Tree): 

    Returns:
        np array (float): matrix of size num_leaves x num_leaves containing the number of nodes between pairs of leaves
            num nodes comes from the get_distance function in ete2 using topology_only

    """
    leaves = tree.get_leaves()
    m = len(leaves)
    adjmat = np.zeros(shape=(m,m), dtype=int)
    for l1 in leaves:  #leaf names are not in numeric order
        i = leaf_idx[l1.name]
        for l2 in leaves:
            j = leaf_idx[l2.name]
            adjmat[i,j] = l1.get_distance(l2, topology_only=True)
    return adjmat.astype(float)


def make_gt2sub():
    """dict of
        key: pairs of (different) genotypes in tuple
        value: what (single) substitution happened (None if >1 sub)
    
    """
    NT4 = np.array(('A','C','G','T'))
    GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
    
    base_code = {nt:int(10**(i-1)) for i,nt in enumerate(NT4)}
    gt_code = {gt:base_code[gt[0]]+base_code[gt[1]] for gt in GTYPE10}
    sub_decode = {base_code[nt1]-base_code[nt2]:(nt1,nt2) for nt1,nt2 in itertools.permutations(NT4,2)}
    gt2sub = {(gt1,gt2):sub_decode.get(gt_code[gt1]-gt_code[gt2]) for gt1,gt2 in itertools.permutations(GTYPE10,2)}
    return gt2sub


def make_sub2tstv():
    """dict of
        key: base/base tuple
        value: tv or ts
    """
    NT4 = np.array(('A','C','G','T'))
    base_code = {'A':0,'C':1,'G':0,'T':1}
    
    tstv = ['ts','tv']
    sub2tstv = {(nt1,nt2):tstv[abs(base_code[nt1]-base_code[nt2])] for nt1,nt2 in itertools.permutations(NT4,2)}
    sub2tstv[None] = None
    return sub2tstv


if __name__ == '__main__':
    import argparse
    
    GTYPE3 = np.array(('RR','RA','AA'))
    GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
    DELTA = 1e-4

    parser = argparse.ArgumentParser()
    subp = parser.add_subparsers(metavar='<command>', help='sub-commands')

    #tview uses tview_main, read_label
    parser_tview = subp.add_parser('tview', help='view tree')
    parser_tview.add_argument('tree', metavar='<nwk>', type=str, help='input tree in Newick format')
    parser_tview.add_argument('-a', metavar='STR', dest='attrs', type=str, help='node attributes to print given by a comma separated list')
    parser_tview.add_argument('-l', metavar='FILE', dest='label', type=str, help='leaves label')
    parser_tview.set_defaults(func=tview_main)

    #compare uses compare_main, tree2adjacency
    parser_compare = subp.add_parser('compare', help='compare tree topology')
    parser_compare.add_argument('-t', metavar='FILE', dest='tree', type=str, nargs='+', required=True, help='input tree(s), in Newick format')
    parser_compare.add_argument('-r', metavar='FILE', dest='ref', type=str, required=True, help='reference tree, in Newick format')
    parser_compare.set_defaults(func=compare_main)

    #compat uses compat_main, read_vcf, calc_compat, find_singleton
    parser_compat = subp.add_parser('compat', help='calculate pairwise compatibility between all pairs of sites')
    parser_compat.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_compat.add_argument('output', metavar='<output>', type=str, help='output compatibility matrix')
    parser_compat.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_compat.set_defaults(func=compat_main)

    #nbjoin uses neighbor_main, read_vcf, make_base_prior (normalize_PL), make_mut_matrix (phred2p, gtype_distance), make_D (pairwise_diff, normalize2d_PL, phred2p), init_star_tree, neighbor_joining
    #init_tree, populate_tree_PL, calc_mut_likelihoods (p2phred), recursive_NNI (nearest_neighbor_interchange, update_PL, score), recursive_reroot (reroot, update_PL)
    parser_nbjoin = subp.add_parser('nbjoin', help='neighbor-joining')
    parser_nbjoin.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_nbjoin.add_argument('output', metavar='output', type=str, help='output basename')
    parser_nbjoin.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_nbjoin.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30')
    parser_nbjoin.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_nbjoin.set_defaults(func=neighbor_main)

    #gtype uses genotype_main
    parser_gtype = subp.add_parser('gtype', help='genotype samples with help of a lineage tree')
    parser_gtype.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_gtype.add_argument('output', metavar='<output>', type=str, help='output basename')
    parser_gtype.add_argument('-t', metavar='FILE', dest='tree', type=str, required=True, help='lineage tree')
    parser_gtype.add_argument('-n', metavar='INT', dest='nsite', type=int, default=1000, help='number of sites processed once, default 1000')
    parser_gtype.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_gtype.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30, 0 for uninformative')
    parser_gtype.set_defaults(func=genotype_main)

    #annot uses annotate_main
    parser_annot = subp.add_parser('annot', help='annotate lineage tree with genotype calls')
    parser_annot.add_argument('gtcall', metavar='<gtcall>', type=str, help='input gtype calls, "-" for stdin')
    parser_annot.add_argument('output', metavar='<outnwk>', type=str, help='output tree in Newick format')
    parser_annot.add_argument('-t', metavar='FILE', dest='tree', type=str, required=True, help='lineage tree')
    parser_annot.set_defaults(func=annotate_main)

    #parser_split = subp.add_parser('split', help='a top-down method that partition samples at a sequence of variants ordered by decreasing MAF')
    #parser_split.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    #parser_split.add_argument('output', metavar='<output>', type=str, help='output basename')
    #parser_split.set_defaults(func=split_main)

    #parser_rsplit = subp.add_parser('rsplit', help='similar to "split" but involves random shuffling of variants instead of ordering by MAF')
    #parser_rsplit.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    #parser_rsplit.add_argument('output', metavar='<output>', type=str, help='output basename')
    #parser_rsplit.add_argument('-n', metavar='INT', dest='n_rep', type=int, default=1000, help='number of random shuffles, default 1000')
    #parser_rsplit.add_argument('-b', metavar='INT', dest='n_bin', type=int, default=1, help='number of MAF bins within each random shuffles occur, default 1')
    #parser_rsplit.set_defaults(func=rsplit_main)

    try:
        args = parser.parse_args()
        args.func(args)
    except KeyboardInterrupt:
        sys.exit(1)
