#!/usr/bin/env python2

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>

from __future__ import print_function
import warnings
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import sys
import itertools
import numpy as np
from scipy.stats import sem
from editdist import distance as strdist
from treecall.memoize import memoized
from treecall.pyvcf import Vcf,VcfFile,warning

with warnings.catch_warnings(ImportWarning):
    from ete2 import Tree

warnings.filterwarnings('error')

NT4 = np.array(('A','C','G','T'))
GTYPE3 = np.array(('RR','RA','AA'))
GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
I3,I10 = len(GTYPE3),len(GTYPE10)
F3,F10 = float(I3),float(I10)
DELTA = 1e-4

def iter_vcf(vcffile):
    vcffile.open()
    if vcffile.seekable:
        line = '#'
    else:
        line = vcffile.fmt_line
    while True:
        if line[0] == '#':
            pass
        else:
            yield Vcf(line)
        try:
            line = vcffile.next()
        except:
            vcffile.close()
            break


def read_vcf(filename, evidence=60):
    print('read_vcf() begin', end=' ', file=sys.stderr)
    vcffile = VcfFile(filename)
    fmt = vcffile.fmt
    strsplit = str.split

    a2g = np.array((
        ((0,0,0), (0,1,2), (0,3,5), (0,6,9)),
        ((2,1,0), (2,2,2), (2,4,5), (2,7,9)),
        ((5,3,0), (5,4,2), (5,5,5), (5,8,9)),
        ((9,6,0), (9,7,2), (9,8,5), (9,9,9)),
    ))

    variants,DPRs,PLs = [],[],[]
    for v in vcffile:
        try:
            variants.append((v.CHROM,v.POS,v.REF))
            dpr = np.array(v.extract_gtype('DPR', fmt, strsplit, ','), dtype=np.uint16)
            ak = dpr.sum(axis=0).argsort(kind='mergesort')[-2:][::-1] # allele ordered by decreasing depth, take only the two most common alleles
            DPRs.append(dpr[...,ak])
            pl = np.array(v.extract_gtype('PL', fmt, strsplit, ','), dtype=np.uint16)
            gk = a2g[ak[0],ak[1]] # take only gtypes formed by the two most common alleles, ordered by increasing PL (decreasing GL)
            PLs.append(pl[...,gk])
        except Exception as e:
            print(e, file=sys.stderr)
            print(v, file=sys.stderr)

    variants = np.array(variants)
    DPRs = np.array(DPRs, dtype=np.uint16)
    PLs = np.array(PLs, dtype=np.uint16)

    k_ev = (PLs.sum(axis=1)>=evidence).sum(axis=1)==3
    variants,DPRs,PLs = variants[k_ev],DPRs[k_ev],PLs[k_ev]

    print(' done', file=sys.stderr)
    return vcffile, variants, DPRs, PLs


def read_vcf_records(vcffile, fmt, maxn=1000):
    print('read next %d sites' % maxn, end=' ', file=sys.stderr)
    variants,DPRs,PLs = [],[],[]
    i = 0
    for v in vcffile:
        i += 1
        try:
            variants.append((v.CHROM,v.POS,v.REF))
            dpr = np.array(v.extract_gtype('DPR', fmt, v.get_DPR4), dtype=np.uint16)
            DPRs.append(dpr)
            pl = np.array(v.extract_gtype('PL', fmt, v.get_PL10), dtype=np.longdouble)
            PLs.append(pl)
        except Exception as e:
            print(e, file=sys.stderr)
            print(v, file=sys.stderr)
        if i == maxn:
            print('... %s:%s ...' % (v.CHROM, v.POS), end=' ', file=sys.stderr)
            break

    variants = np.array(variants)
    DPRs = np.array(DPRs, dtype=np.uint16)
    PLs = np.array(PLs, dtype=np.longdouble)

    print(' done', file=sys.stderr)
    return variants, DPRs, PLs


def compat_main(args):
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    #n_site, n_smpl = PLs.shape[0:2]
    #sidx = np.arange(n_smpl)

    compats = calc_compat(PLs)
    c = compats.sum(axis=-1)
    k = (c == 0) & ~find_singleton(PLs)

    gzout = args.output + '.gz'
    np.savetxt(gzout, compats, fmt='%d', delimiter='\t')


def calc_compat(PLs):
    print('calc_compat() begin', end=' ', file=sys.stderr)
    n,m,g = PLs.shape
    nidx = np.arange(n)
    midx = np.arange(m)
    kn = np.tile(nidx,m).reshape(m,n)
    km = np.repeat(midx,n).reshape(m,n)

    PLs = PLs.astype(int)
    gt = (PLs[...,0]==0).astype(np.byte) # n x m
    non_zeros = PLs[kn,km,gt.T].T # n x m
    groups = (2*gt[i]+gt for i in xrange(n)) # each n x m
    cost = (np.minimum(non_zeros[i], non_zeros) for i in xrange(n)) # each n x m

    compats = np.zeros(shape=(n,n), dtype=np.int32)
    for i in xrange(n):
        grp = groups.next()
        cst = cost.next()
        compats[i,i:] = map(min, map(np.bincount, grp[i:], cst[i:]))
    compats = compats + compats.T - np.diag(compats.diagonal())

    print(' done', file=sys.stderr)
    return compats


def find_singleton(PLs):
    n,m,b = PLs.shape
    is_singleton = (PLs[...,0]>0).sum(axis=1)==1
    return is_singleton


def neighbor_main(args):
    print(args, file=sys.stderr)
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    base_prior = make_base_prior(args.het, GTYPE3) # base genotype prior
    mm,mm0,mm1 = make_mut_matrix(args.mu, GTYPE3) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)
    tree = init_star_tree(n_smpl)
    internals = np.arange(n_smpl)
    neighbor_joining(D, tree, internals)

    init_tree(tree)
    populate_tree_PL(tree, PLs, mm0, 'PL0')
    calc_mut_likelihoods(tree, mm0, mm1)

    print(tree)
    tree.write(outfile=args.output+'.nj0.nwk', format=5)
    best_tree,best_PL = recursive_NNI(tree, mm0, mm1, base_prior)
    print(best_tree)
    best_tree,best_PL = recursive_reroot(best_tree, mm0, mm1, base_prior)
    print(best_tree)
    print('PL_per_site = %.4f' % (best_PL/n_site))
    best_tree.write(outfile=args.output+'.nj.nwk', format=5)


def init_star_tree(n):
    tree = Tree()
    for i in xrange(n):
        tree.add_child(name=str(i))
    return tree


def pairwise_diff(PLs, i, j):
    pli = normalize2d_PL(PLs[:,i])
    plj = normalize2d_PL(PLs[:,j])
    p = phred2p(pli+plj) # n x g
    return (1-p.sum(axis=1)).sum()


def make_D(PLs):
    n,m,g = PLs.shape
    D = np.zeros(shape=(2*m-2,2*m-2), dtype=np.longdouble)
    for i,j in itertools.combinations(xrange(m),2):
        D[i,j] = pairwise_diff(PLs, i, j)
        D[j,i] = D[i,j]
    return D


def neighbor_joining(D, tree, internals):
    print('neighbor_joining() begin', end=' ', file=sys.stderr)
    m = len(internals)
    while m > 2:
        d = D[internals[:,None],internals]
        u = d.sum(axis=1)/(m-2)

        Q = np.zeros(shape=(m,m), dtype=np.longdouble)
        for i,j in itertools.combinations(xrange(m),2):
            Q[i,j] = d[i,j]-u[i]-u[j]
            Q[j,i] = Q[i,j]
        #print(Q.astype(int))
        np.fill_diagonal(Q, np.inf)
        #print(np.unique(Q, return_counts=True))
        i,j = np.unravel_index(Q.argmin(), (m,m))
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
    return tree


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


def partition_main(args):
    print(args, file=sys.stderr)
    base_prior = make_base_prior(args.het, GTYPE3) # base genotype prior
    mm,mm0,mm1 = make_mut_matrix(args.mu, GTYPE3) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    n_site,n_smpl = PLs.shape[0:2]

    tree = Tree()
    if sem(PLs[...,1],axis=1).mean() > sem(PLs[...,2],axis=1).mean():
        partition(PLs[...,0:2], tree, np.arange(n_smpl), args.min_ev)
    else:
        partition(PLs, tree, np.arange(n_smpl), args.min_ev)

    init_tree(tree)
    PLs = PLs.astype(np.longdouble)
    populate_tree_PL(tree, PLs, mm, 'PL')
    populate_tree_PL(tree, PLs, mm0, 'PL0')
    calc_mut_likelihoods(tree, mm0, mm1)

    print(tree)
    tree.write(outfile=args.output+'.pt0.nwk', format=5)
    best_tree,best_PL = recursive_NNI(tree, mm0, mm1, base_prior)
    best_tree,best_PL = recursive_reroot(best_tree, mm0, mm1, base_prior)
    print(best_tree)
    print('PL_per_site = %.4f' % (best_PL/n_site))
    best_tree.write(outfile=args.output+'.pt.nwk', format=5)


def genotype_main(args):
    print(args, file=sys.stderr)

    tree = Tree(args.tree)
    init_tree(tree)

    base_prior = make_base_prior(args.het, GTYPE10) # base genotype prior
    mm,mm0,mm1 = make_mut_matrix(args.mu, GTYPE10) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    vcffile = VcfFile(args.vcf)
    fmt = vcffile.fmt
    fout = open(args.output, 'w')
    fout.close()
    fout = open(args.output, 'a')
    score = 0
    while True:
        variants, DPRs, PLs = read_vcf_records(vcffile, vcffile.fmt, args.nsite)
        records,s = genotype(PLs, tree, variants, mm, mm0, mm1, base_prior)
        np.savetxt(fout, records, fmt=['%s','%d','%s','%.2e','%.2e','%s','%.2e','%s','%s','%.2e','%d','%s'], delimiter='\t')
        score += s
        if len(PLs) < args.nsite:
            break
    print('sum(PL) = %.2f' % score)
    fout.close()


def genotype(PLs, tree, variants, mm, mm0, mm1, base_prior):
    # calculate total likelihoods for each genotypes
    populate_tree_PL(tree, PLs, mm, 'PL') # dim(tree.PL) = site x gtype
    tree_PL = tree.PL + base_prior
    # calculate no-mutation likelihoods for each genotypes
    #try:
    populate_tree_PL(tree, PLs, mm0, 'PL0') # dim(tree.PL0) = site x gtype
    #except Exception as e:
    #    print('populate_tree_PL():', e, file=sys.stderr)
    #    sys.exit(1)
    tree_PL0 = tree.PL0 + base_prior
    # calculate mutation likelihoods for each genotypes and mutation locations
    calc_mut_likelihoods(tree, mm0, mm1)
    mut_PLs = np.swapaxes(tree.PLm,0,1) # site x location x gtype
    mut_PLs += base_prior
    n,l,g = mut_PLs.shape # n sites, l locations, g gtypes
    nn = np.arange(n)

    k = tree_PL.argmin(axis=1) # most likely base genotype for each site
    tree_P_per_site = phred2p(tree_PL).sum(axis=1) # total tree likelihood

    k0 = tree_PL0.argmin(axis=1) # most likely non-mutation base genotype for each site
    null_PL = tree_PL0[nn,k0] # best non-mutation likelihood (across genotypes) for each site
    null_P_per_site = phred2p(tree_PL0).sum(axis=1) # total non-mutation likelihood

    k1 = np.array([np.unravel_index(s.argmin(), (l,g)) for s in mut_PLs]) # site x 2, most likely mutation event for each site
    k1l = k1[:,0] # most likely location
    k1g = k1[:,1] # most likely base genotype
    mut_PL = mut_PLs[nn,k1l,k1g] # best mutation likelihood (across location and genotypes) for each site
    mut_P_per_site = phred2p(mut_PLs).sum(axis=(1,2)) # total mutation likelihood

    null_PLs = np.array([node.PL0 for node in tree.iter_descendants(strategy='postorder')])
    k2 = null_PLs[k1l,nn,].argmin(axis=-1) # get most likely mutation mutant genotype

    node_sids = np.array([','.join(map(str,node.sid)) for node in tree.iter_descendants(strategy='postorder')])
    records = np.array(zip(
            variants[nn,0],                         # chrom
            variants[nn,1],                         # pos
            variants[nn,2],                         # ref
            null_P_per_site/tree_P_per_site,        # null_P
            mut_P_per_site/tree_P_per_site,         # mut_P
           #GTYPE10[k],                             # MLE_base_gtype
           #phred2p(tree_PL[nn,k])/tree_P_per_site, # MLE_base_gtype_P
            GTYPE10[k0],                            # MLE_null_base_gtype
            phred2p(null_PL)/tree_P_per_site,       # MLE_null_base_gtype_P
            GTYPE10[k1g],                           # MLE_mut_base_gtype
            GTYPE10[k2],                            # MLE_mut_alt_gtype
            phred2p(mut_PL)/tree_P_per_site,        # MLE_mut_base_gtype_P
            k1l,                                    # MLE_mut_location
            node_sids[k1l]),                        # MLE_mut_samples
        dtype=[
            ('chrom','a10'),('pos','i4'),('ref','a1'),
            ('null_p','f8'),('mut_p','f8'),
            ('null_base','a2'),('null_base_p','f8'),
            ('mut_base','a2'),('mut_alt','a2'),('mut_conf_p','f8'),
            ('mut_loc','i4'),('mut_smpl','a128')])
    score = p2phred(records['mut_p']+records['null_p']).sum()
    return records,score


def init_tree(tree):
    tree.leaf_order = map(int, tree.get_leaf_names())

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.sid = [int(node.name)]
        else:
            node.name = ''
            node.sid = []
            for child in node.children:
                node.sid.extend(child.sid)

    m = len(tree)
    for i,node in zip(xrange(2*m-1), tree.traverse(strategy='postorder')):
        node.nid = i
        node.sid = sorted(node.sid)


def p2phred(x):
    return -10.0*np.log10(x)

def phred2p(x):
    return 10.0**(-x/10.0)

def sum_PL(x, axis=None):
    return p2phred(phred2p(x).sum(axis=axis))

def normalize_PL(x):
    p = 10.0**(-x/10.0)
    return -10.0*np.log10(p/p.sum())

def normalize2d_PL(x):
    p = 10.0**(-x/10.0)
    return -10.0*np.log10(p/p.sum(axis=1)[:,None])


def gtype_distance(gt):
    n = len(gt)
    gt_dist = np.zeros((n,n), dtype=int)
    for i,gi in enumerate(gt):
        for j,gj in enumerate(gt):
            gt_dist[i,j] = min(strdist(gi,gj),strdist(gi,gj[::-1]))
    return gt_dist


def make_mut_matrix(mu, gtypes):
    pmu = phred2p(mu)
    gt_dist = gtype_distance(gtypes)
    mm = pmu**gt_dist
    np.fill_diagonal(mm, 2.0-mm.sum(axis=0))
    mm0 = np.diagflat(mm.diagonal()) # substitution rate matrix with non-diagonal set to 0
    mm1 = mm - mm0 # substitution rate matrix with diagonal set to 0
    return mm,mm0,mm1


def make_base_prior(het, gtypes):
    return normalize_PL(np.array([g[0]!=g[1] for g in gtypes], dtype=np.longdouble)*het)


def calc_mut_likelihoods(tree, mm0, mm1):
    n,g = tree.PL0.shape
    for node in tree.traverse(strategy='postorder'):
        if not node.is_leaf():
            node.PLm = np.zeros((2*len(node)-2,n,g), dtype=np.longdouble)

    for node in tree.traverse(strategy='postorder'):
        i = 0
        for child in node.children:
            sister = child.get_sisters()[0]
            if not child.is_leaf():
                l = child.PLm.shape[0]
                node.PLm[i:(i+l)] = p2phred(np.dot(phred2p(child.PLm), mm0)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
                i += l
            node.PLm[i] = p2phred(np.dot(phred2p(child.PL0), mm1)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
            i += 1


def update_PL(node, mm0, mm1):
    n,g = node.PL0.shape
    l = 2*len(node)-2
    #node.PL0 = np.zeros((n,g), dtype=np.longdouble)
    node.PL0.fill(0.0)
    node.PLm = np.zeros((l,n,g), dtype=np.longdouble)
    for child in node.children:
        sid = sorted(map(int,child.get_leaf_names()))
        if child.sid != sid:
            update_PL(child, mm0, mm1)
            child.sid = sid
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


def populate_tree_PL(tree, PLs, mm, attr):
    n,m,g = PLs.shape # n sites, m samples, g gtypes
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            setattr(node, attr, PLs[:,node.sid[0],])
        else:
            setattr(node, attr, np.zeros((n,g), dtype=np.longdouble))
            for child in node.children:
                setattr(node, attr, getattr(node, attr) + p2phred(np.dot(phred2p(getattr(child, attr)), mm)))

def score(tree, base_prior):
    Pm = phred2p(tree.PLm+base_prior).sum(axis=(0,2))
    P0 = phred2p(tree.PL0+base_prior).sum(axis=1)
    return p2phred(Pm+P0).sum()


def annotate_nodes(tree, attr, values):
    for node in tree.iter_descendants('postorder'):
        setattr(node, attr, values[node.nid])


def tview_main(args):
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
    init_tree(tree)

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


def partition(PLs, tree, sidx, min_ev):
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


def reroot(tree, mm0, mm1, base_prior):
    '''
              /-A              /-A              /-B
           /-|              /-|              /-|
    -root-|   \-B => -root-|   \-C => -root-|   \-C
          |                |                |
           \-C              \-B              \-A
    '''
    best_tree = tree
    best_PL = score(tree, base_prior)

    for node in tree.iter_descendants('postorder'):
        tree_reroot = tree.copy()
        new_root = tree_reroot.search_nodes(sid=node.sid)[0]
        tree_reroot.set_outgroup(new_root)
        update_PL(tree_reroot, mm0, mm1)
        PL_reroot = score(tree_reroot, base_prior)
        #print(tree_reroot)
        #print(PL_reroot)
        if PL_reroot < best_PL * (1-DELTA):
            best_tree = tree_reroot
            best_PL = PL_reroot
    return best_tree,best_PL


def recursive_reroot(tree, mm0, mm1, base_prior):
    print('recursive_reroot() begin', end=' ', file=sys.stderr)
    for node in tree.iter_descendants('postorder'):
        if node.is_leaf():
            continue
        print('.', end='', file=sys.stderr)
        new_node,new_PL = reroot(node,mm0,mm1,base_prior)
        parent = node.up
        parent.remove_child(node)
        parent.add_child(new_node)
        update_PL(tree, mm0, mm1)
    new_tree,new_PL = reroot(tree, mm0, mm1, base_prior)
    print(' done', end='', file=sys.stderr)
    #print(new_tree)
    #print(new_PL)
    return new_tree,new_PL


def nearest_neighbor_interchange(node, mm0, mm1, base_prior):
    '''
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
    c1,c2 = node.children
    if c1.is_leaf() and c2.is_leaf():
        return None,None
    if c1.is_leaf() or c2.is_leaf():
        return reroot(node, mm0, mm1, base_prior)

    #conf0
    node_copy0 = node.copy()
    node0,PL0 = reroot(node_copy0, mm0, mm1, base_prior)
    #conf1
    node_copy1 = node.copy()
    c1,c2 = node_copy1.children
    c11,c12 = c1.children
    c21,c22 = c2.children
    c12 = c12.detach()
    c22 = c22.detach()
    c1.add_child(c22)
    c2.add_child(c12)
    update_PL(node_copy1, mm0, mm1)
    node1,PL1 = reroot(node_copy1, mm0, mm1, base_prior)
    #conf2
    node_copy2 = node.copy()
    c1,c2 = node_copy2.children
    c11,c12 = c1.children
    c21,c22 = c2.children
    c12 = c12.detach()
    c21 = c21.detach()
    c1.add_child(c21)
    c2.add_child(c12)
    update_PL(node_copy2, mm0, mm1)
    node2,PL2 = reroot(node_copy2, mm0, mm1, base_prior)

    if PL1 < PL0 * (1-DELTA):
        if PL1 < PL2:
            return node1,PL1
        else:
            return node2,PL2
    if PL2 < PL0 * (1-DELTA):
        return node2,PL2
    else:
        return node0,PL0


def recursive_NNI(tree, mm0, mm1, base_prior):
    print('recursive_NNI() begin', end=' ', file=sys.stderr)
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            continue
        print('.', end='', file=sys.stderr)
        node_nni,PL_nni = nearest_neighbor_interchange(node, mm0, mm1, base_prior)
        if node_nni is None:
            continue
        if node.is_root():
            tree = node_nni
            PL = PL_nni
        else:
            parent = node.up
            node.detach()
            parent.add_child(node_nni)
            update_PL(tree, mm0, mm1)
            PL = score(tree, base_prior)
    print(' done', file=sys.stderr)
    #print(tree)
    #print(PL)
    return tree,PL


def compare_main(args):
    print(args, file=sys.stderr)
    ref_tree = Tree(args.ref)
    ref_am = tree2adjacency(ref_tree)
    for f in args.tree:
        tree = Tree(f)
        am = tree2adjacency(tree)
        if ref_am.shape != am.shape:
            print('%s incompatible with %s' % (f, args.ref), file=sys.stderr)
        else:
            k = ref_am > 0

            diff = np.abs(ref_am - am)
            dstat = diff[k].sum()/k.sum()

            ratio = am[k]/ref_am[k]
            ratio[ratio>1] = 1.0/ratio[ratio>1]
            rstat = np.power(ratio.prod(), 1.0/k.sum())

            result = ref_tree.compare(tree, unrooted=True)

            # <tree>,<norm_rf>,<ref_edge_in_tree>,<tree_edge_in_ref>,<diff_adj>,<ratio_adj>
            print('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (f, result['norm_rf'], result['ref_edges_in_source'], result['source_edges_in_ref'], dstat, rstat))


def tree2adjacency(tree):
    leaves = tree.get_leaves()
    m = len(leaves)
    adjmat = np.zeros(shape=(m,m), dtype=int)
    for l1 in leaves:
        i = int(l1.name)
        for l2 in leaves:
            j = int(l2.name)
            adjmat[i,j] = l1.get_distance(l2, topology_only=True)
    return adjmat.astype(float)


def make_gt2sub():
    base_code = {nt:int(10**(i-1)) for i,nt in enumerate(NT4)}
    gt_code = {gt:base_code[gt[0]]+base_code[gt[1]] for gt in GTYPE10}
    sub_decode = {base_code[nt1]-base_code[nt2]:(nt1,nt2) for nt1,nt2 in itertools.permutations(NT4,2)}
    gt2sub = {(gt1,gt2):sub_decode.get(gt_code[gt1]-gt_code[gt2]) for gt1,gt2 in itertools.permutations(GTYPE10,2)}
    return gt2sub


def make_sub2tstv():
    base_code = {'A':0,'C':1,'G':0,'T':1}
    tstv = ['ts','tv']
    sub2tstv = {(nt1,nt2):tstv[abs(base_code[nt1]-base_code[nt2])] for nt1,nt2 in itertools.permutations(NT4,2)}
    sub2tstv[None] = None
    return sub2tstv


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    subp = parser.add_subparsers(metavar='<command>', help='sub-commands')

    parser_tview = subp.add_parser('tview', help='view tree')
    parser_tview.add_argument('tree', metavar='<nwk>', type=str, help='input tree in Newick format')
    parser_tview.add_argument('-a', metavar='STR', dest='attrs', type=str, help='node attributes to print given by a comma separated list')
    parser_tview.add_argument('-l', metavar='FILE', dest='label', type=str, help='leaves label')
    parser_tview.set_defaults(func=tview_main)

    parser_compare = subp.add_parser('compare', help='compare tree topology')
    parser_compare.add_argument('-t', metavar='FILE', dest='tree', type=str, nargs='+', required=True, help='input tree(s), in Newick format')
    parser_compare.add_argument('-r', metavar='FILE', dest='ref', type=str, required=True, help='reference tree, in Newick format')
    parser_compare.set_defaults(func=compare_main)

    parser_compat = subp.add_parser('compat', help='calculate pairwise compatibility between all pairs of sites')
    parser_compat.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_compat.add_argument('output', metavar='<output>', type=str, help='output compatibility matrix')
    parser_compat.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_compat.set_defaults(func=compat_main)

    parser_nbjoin = subp.add_parser('nbjoin', help='neighbor-joining')
    parser_nbjoin.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_nbjoin.add_argument('output', metavar='output', type=str, help='output basename')
    parser_nbjoin.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_nbjoin.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30')
    parser_nbjoin.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_nbjoin.set_defaults(func=neighbor_main)

    parser_part = subp.add_parser('part', help='a top-down method that partition samples by sum of partition cost across all sites')
    parser_part.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_part.add_argument('output', metavar='<output>', type=str, help='output basename')
    parser_part.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_part.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30')
    parser_part.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_part.set_defaults(func=partition_main)

    parser_gtype = subp.add_parser('gtype', help='genotype samples with help of a lineage tree')
    parser_gtype.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_gtype.add_argument('output', metavar='<output>', type=str, help='output basename')
    parser_gtype.add_argument('-t', metavar='FILE', dest='tree', type=str, required=True, help='lineage tree')
    parser_gtype.add_argument('-n', metavar='INT', dest='nsite', type=int, default=1000, help='number of sites processed once, default 1000')
    parser_gtype.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_gtype.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30, 0 for uninformative')
    parser_gtype.set_defaults(func=genotype_main)

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
