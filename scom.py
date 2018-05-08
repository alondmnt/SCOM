# -*- coding: utf-8 -*-
"""
code for detecting Spatially Co-evolving Orthologous Modules (SCOMs)

Alon Diament, Tuller Lab.
"""

import os
import sys
import warnings
import scipy.io as sio
import numpy as np
from collections import namedtuple
from copy import deepcopy
import time


### for the sake of the example set, min_seed was set to 0 (20 originally) ###
DefModes = {'target': 'weight', 'total': 'sum', 'min_delta': 1e-3,
            'max_overlap': 0.2, 'min_size': 5, 'max_size': 30,
            'min_seed': 0, 'max_moves': 20e3, 'rand_seed': 0,
            'alpha': 0.9, 'cleanup': True, 'add': True, 'remove': True}

#DefModes['rand_seed'] = int(time.time())
np.random.seed(DefModes['rand_seed'])


def calc_modules(G, seeds, R=None, params=deepcopy(DefModes)):
    """G - green graph, R (optional) - red graph"""
    params = {k: v for k, v in params.items() if type(v) is not str}
    params['ignore_R'] = R is None  # for computing regular G-based modules
    params = namedtuple('Mode', params.keys())(**params)
    print(params)

    if R is None:
        V, I, U, init_score = init_monoseeds(seeds, G, params)
        print('{} single-module seeds (init_score = {:.2f})'.format(V.shape[1],
              init_score.mean()))
    else:
        V, I, U, init_score = init_biseeds(seeds, G, R, params)
        print('{} bi-module seeds (init_score = {:.2f})'.format(V.shape[1],
              init_score.mean()))

    move_count = np.zeros((4), dtype=int)
    move_success = np.ones((3), dtype=int)
    i_add = 0
    i_rep = 1
    i_del = 2
    i_iter = 3

    while np.sum(move_success):
        move_success[:] = 0
        if np.any(move_count > params.max_moves):
            print('exceeded optimization steps')
            break

        modperm = np.arange(V.shape[1])
        np.random.shuffle(modperm)  # no bias
        for m in modperm:
            if params.add:
                move_success[i_add] += add_bimodule_node(V, m, I, U, G, R, params)

            if params.remove:
                move_success[i_del] += del_bimodule_node(V, m, I, U, G, R, params)

        move_count[:3] += move_success
        move_count[i_iter] += 1

    if np.any(I / U > params.max_overlap):
        warnings.warn('overlap exceeded')

    print('[add: {}, replace: {}, remove: {}, iters: {}]'.format(*move_count))

    return tidy_bimodules(V, G, R, init_score, params)


def build_weighted_graph(G, alpha=0.9, p_null=None):
    """using log-odds ratio probability weighting to model graph edges.
    when [p_null] is missing, using a default null model for the weighting:
    random degree-preserving rewiring approximation."""
    W = G > 0
    deg = W.sum(axis=1).reshape((1, -1))
    M = W.sum()
    if p_null is None:
        print('build_weighted_graph: using degree-preserving theoretical approx.')
        p_null = np.minimum(1-1e-4, deg * deg.T / M)  # degree preserving rewiring
    W = np.log(alpha*W + (1-alpha)*(1-W)) - np.log(p_null*W + (1-p_null)*(1-W))
    # where we have an edge W_ij = log(alpha / p_null)
    # where we don't have an edge W_ij = log((1-alpha)/(1-p_null))
    if np.isinf(W).any():
        raise Exception('ValueError: Inf in weighted graph')
    np.fill_diagonal(W, 0)
    return W


def calc_charikar_all_subgraphs(G, modes=deepcopy(DefModes)):
    """we implement Charikar's 2-approximation for finding the densest subgraph
    in graph [G] and run it repeatedly to get all such subgraphs.
    possible modes: target=weight/degree, total=sum/average."""
    sets = []
    scores = []
    while G.any():
        this, score = calc_charikar_densest(G, modes)[:2]
        if len(this) == 0:
            break
        G[this, :] = 0
        G[:, this] = 0
        sets.append(this)
        scores.append(score)
    return sets, scores


def calc_charikar_densest(G, modes):
    """we implement Charikar's 2-approximation for a single densest subgraph in
    [G]."""
    if modes['target'] == 'degree':
        G = G > 0
    elif modes['target'] == 'weight':
        G = G.copy()
    target = G.sum(axis=1, dtype=np.float)
    nodes_left = list(target.nonzero()[0])
    nodes_out = []
    nV = len(nodes_left)
    total = G.sum() / 2

    if nV >= modes['min_size'] and nV <= modes['max_size']:
        max_dense = total
        if modes['total'] == 'average':
            max_dense /= nV
        maxV = list(nodes_left)
    elif nV > modes['max_size']:
        max_dense = -np.inf
        maxV = []
    elif nV < modes['min_size']:
        return [], None

    while len(nodes_left) > 0:
        i = target[nodes_left].argmin()
        v = nodes_left[i]
        del nodes_left[i]
        nodes_out.append(v)

        total -= target[v]
        neis = [nodes_left[n] for n in G[v, nodes_left].nonzero()[0]]
        target[neis] -= G[v, neis]
        G[v, neis] = 0
        G[neis, v] = 0
        nV -= 1

        cur_dense = total
        if modes['total'] == 'average':
            cur_dense /= nV
        if cur_dense >= max_dense and nV >= modes['min_size'] and \
                                      nV <= modes['max_size']:
            # for equal value we prefer minimal size
            maxV = list(nodes_left)
            max_dense = cur_dense
    nodes_out.reverse()

    return maxV, max_dense, nodes_out


def bimodule_objective(v1, v2, G, R, ignore_R=False):

    # intra GREEN score
    obj = G[v1, :][:, v1].sum()
    obj += G[v2, :][:, v2].sum()

    if not ignore_R:
        # inter RED score
        obj += 2*R[v1, :][:, v2].sum()  # symmetry (we count twice above)

    return obj


def update_overlap(I, U, V, m, v, action=1):
    """ I: interscection matrix
        U: union matrix
        V: all modules
        m: module index
        v: vertex index
        action: add/remove"""
    inmod = V[0][:, v] + V[1][:, v]  # v is in module
    if np.sum(inmod):
        I[m][inmod] = I[m][inmod] + action  # intersection increased
        I[:, m][inmod] = I[:, m][inmod] + action
        I[m, m] = 0

    outmod = ~inmod
    if np.sum(outmod):
        U[m][outmod] = U[m][outmod] + action  # union increased
        U[:, m][outmod] = U[:, m][outmod] + action
        U[m, m] = 1


def init_monoseeds(seeds, G, params):
    """ adapted from init_bimodules. """
    # CONVERT SEEDS
    N = len(seeds)
    # we init V to have the same dimensions as bimodules, but 
    V = np.zeros((2, N, G.shape[0]), dtype=bool)  # submodules x seeds x nodes
    init_score = np.zeros((N))
    for i, s in enumerate(seeds):
        V[0][i][s] = True
        init_score[i] = bimodule_objective(V[0][i], V[1][i], G, None, ignore_R=True)

    # OVERLAP
    inmod = V[0].astype(np.int16)  # [i,j] = 1 >> node j in module i
    I = inmod.dot(inmod.T)  # intersect
    inmod = np.reshape(inmod.sum(axis=1), (1, -1))
    U = inmod + inmod.T - I  # union: A + B - (A^B)
    np.fill_diagonal(I, 0)
    np.fill_diagonal(U, 1)

    if params.cleanup:
        # FILTER OVERLAP (discard)
        while (I / U > params.max_overlap).any():
            # remove module with largest total overlap with others
            ind = np.argmax(np.sum(I / U, axis=1))
            V = np.delete(V, ind, 1)
            I = np.delete(np.delete(I, ind, 0), ind, 1)
            U = np.delete(np.delete(U, ind, 0), ind, 1)
            init_score = np.delete(init_score, ind, 0)

    return V, I, U, init_score


def init_biseeds(seeds, G, R, params):
    # CONVERT SEEDS
    N = len(seeds)
    tmp = np.zeros((N, G.shape[0]), dtype=bool)  # seeds x nodes
    for i, s in enumerate(seeds):
        tmp[i][s] = True
    seeds = tmp

    # ALL SEED-PAIRS SCORES
    cand_score = np.full((N, N), -np.inf)
    for i in range(N):
        sz = np.sum(seeds[i])
        if sz > params.max_size:
            print('seed {} too large'.format(i))
            continue
        if sz < params.min_size:
            print('seed {} too small'.format(i))
            continue
        for j in range(i):
            sz = np.sum(seeds[j])
            if sz > params.max_size:
                continue
            if sz < params.min_size:
                continue
            if np.sum(np.logical_and(seeds[i], seeds[j])):
                continue
            cand_score[i, j] = bimodule_objective(seeds[i], seeds[j], G, R)
            cand_score[j, i] = cand_score[i, j]

    mod_per_seed = 10  # take k-most optimized pairs, later filtered for overlaps
    init_score = np.sort(cand_score, axis=1)[:, -mod_per_seed:].flatten(order='F')
    v1 = np.arange(N)
    v2 = np.argsort(cand_score, axis=1)[v1, -mod_per_seed:].flatten(order='F')
    v1 = np.concatenate(mod_per_seed*(v1,))

    # MODULES FROM SEED PAIRS
    # [s, i, j] = True >> node j in submodule s of module i
    V = np.zeros((2, len(v1), G.shape[0]), dtype=bool)  # submodules x modules x nodes
    for i in range(len(v1)):
        if cand_score[v1[i], v2[i]] < params.min_seed:
            continue
        V[0][i][seeds[v1[i]]] = True
        V[1][i][seeds[v2[i]]] = True
    V = V[:, np.any(V[0], axis=1), :]

    # OVERLAP
    inmod = (V[0] + V[1]).astype(np.int16)  # [i,j] = True >> node j in module i
    I = inmod.dot(inmod.T)  # intersect
    inmod = np.reshape(inmod.sum(axis=1), (1, -1))
    U = inmod + inmod.T - I  # union: A + B - (A^B)
    np.fill_diagonal(I, 0)
    np.fill_diagonal(U, 1)

    # FILTER OVERLAP (discard)
    while (I / U > params.max_overlap).any():
        # remove module with largest total overlap with others
        ind = np.argmax(np.sum(I / U, axis=1))
        V = np.delete(V, ind, 1)
        I = np.delete(np.delete(I, ind, 0), ind, 1)
        U = np.delete(np.delete(U, ind, 0), ind, 1)
        init_score = np.delete(init_score, ind, 0)

    return V, I, U, init_score


def tidy_bimodules(V, G, R, init_score, params):
    """ packing the modules in two formats:
        [modules] dict that can be exported to a cell array of matlab structs.
        as well [sets1], [sets2] arrays (for the 2 submodules) and [scores]. """
    modules = []
    sets1 = []
    sets2 = []
    scores = []
    for i in range(V.shape[1]):
        S1 = np.nonzero(V[0][i])[0] + 1  # NOTE!: matlab indexing
        S2 = np.nonzero(V[1][i])[0] + 1
        sets1.append(S1)
        sets2.append(S2)
        scores.append(bimodule_objective(V[0][i], V[1][i], G, R, params.ignore_R))
        modules.append({'v1': S1,
                        'v2': S2,
                        'init_score': init_score[i],
                        'score': scores[-1],
                        'G1_score': (G[np.ix_(V[0][i], V[0][i])] > 0).sum(),
                        'G2_score': (G[np.ix_(V[1][i], V[1][i])] > 0).sum()})
        if R is not None:
            modules[-1]['R_score'] = (R[np.ix_(V[0][i], V[1][i])] > 0).sum()

    return modules, sets1, sets2, scores


def add_bimodule_node(V, m, I, U, G, R, params):
    """ V: modules 3-dim array (submodules x modules x graph nodes)
        m: current module index
        I: size of intersection between all modules
        U: size of union between all modules
        G: green graph edges
        R: red graph edges
        params: optimization prarmeters"""
    res = 0

    if np.sum(V[:, m]) >= 2*params.max_size:
        # reacehd max size
        return res

    N = G.shape[0]
    refscore = bimodule_objective(V[0][m], V[1][m], G, R, params.ignore_R)

    subperm = np.arange(V.shape[0])
    np.random.shuffle(subperm)  # no bias
    for i in subperm:  # for each submodule
        if np.sum(V[i][m]) >= params.max_size or np.sum(V[i][m]) == 0:
            continue
        cand_score = np.full((N), -np.inf)
        for cand in range(N):  # for each node
            if np.sum(V[:, m, cand]):
                # already in module
                continue
            inmod = V[0] + V[1]
            if (params.max_overlap) < 1 and np.sum((I[m][inmod[:, cand]] + 1) /
                                                   U[m][inmod[:, cand]]
                                                   > params.max_overlap):
                # addition will increase overlap too much
                # (for the observed modules overlap might increase due to
                # increased intersection)
                continue

            V[i, m, cand] = True  # add momentarily submodule i
            cand_score[cand] = bimodule_objective(V[0][m], V[1][m], G, R, params.ignore_R)
            V[i, m, cand] = False

        if np.max(cand_score) - refscore < params.min_delta:
            continue
        cand = np.argmax(cand_score)
        update_overlap(I, U, V, m, cand, 1)
        V[i, m, cand] = True
        res += 1

    return res


def del_bimodule_node(V, m, I, U, G, R, params):
    """ V: modules 3-dim array (submodules x modules x graph nodes)
        m: current module index
        I: size of intersection between all modules
        U: size of union between all modules
        G: green graph edges
        R: red graph edges
        params: optimization prarmeters"""
    res = 0

    N = G.shape[0]
    refscore = bimodule_objective(V[0][m], V[1][m], G, R, params.ignore_R)

    subperm = np.arange(V.shape[0])
    np.random.shuffle(subperm)  # no bias
    for i in subperm:  # for each submodule
        if np.sum(V[i][m]) <= params.min_size:
            continue
        cand_score = np.full((N), -np.inf)
        inmod = V[0] + V[1]
        for cand in np.nonzero(V[i][m])[0]:  # for each contained vertex
            if (params.max_overlap < 1) and np.sum(I[m][~inmod[:, cand]] /
                                                   (U[m][~inmod[:, cand]] - 1)
                                                   > params.max_overlap):
                # deletion will increase overlap too much
                # (for the observed modules overlap might increase due to
                # decreased union)
                continue
            V[i, m, cand] = False  # mask momentarily
            cand_score[cand] = bimodule_objective(V[0][m], V[1][m], G, R, params.ignore_R)
            V[i, m, cand] = True

        if np.max(cand_score) - refscore < params.min_delta:
            continue
        best = np.argmax(cand_score)
        update_overlap(I, U, V, m, best, -1)
        V[i, m, best] = False
        res += 1

    return res


def run_bimodules(infile, modes=deepcopy(DefModes), empiric_null=None):
    """ detect dense SEPARATED subgraphs.
        INPUT:
            [infile]: matlab datafile containing a graphs [G] of green edges
                    and [R] of red edges
            [modes]: optional, a dictionary of program parameters
            [empiric_null]: optional, list of datafiles, each containing
                    matrix [P] of probabilties for observing an edge between
                    nodes (i, j) in [G] and [R], respectively
        OUTPUT: (written to file)
            [submodules1]: cell array containing subsets of indices
                (NOTE!: matlab indices starting from 1)
            [submodules2]: cell array containing subsets of indices
            [scores]: log odds ratio scores for the modules
            [modes]: program parameters used
            [datafiles]: data input to program
    """
    print('detecting sepSCOMs')
    outfile = '%s sepSCOM' % infile

    G = sio.loadmat(infile)
    R = sio.loadmat(infile)['R']
    G = G['G']

    if not empiric_null:
        Pg = None
        Pr = None
        empiric_null = np.nan
    else:
        print('\nloading {}'.format(empiric_null[0]))
        Pg = sio.loadmat(empiric_null[0])['P']
        print('\nloading {}'.format(empiric_null[1]))
        Pr = sio.loadmat(empiric_null[1])['P']
        epsilon = 1e-4
        Pg = np.maximum(epsilon, np.minimum(1-epsilon, Pg))
        Pr = np.maximum(epsilon, np.minimum(1-epsilon, Pr))

    Wg = build_weighted_graph(G, modes['alpha'], Pg)
    Wr = build_weighted_graph(R, modes['alpha'], Pr)

    itime = time.time()
    seeds, scores = calc_charikar_all_subgraphs(Wg, modes)
    print('seeds took {:.2f} secs, {} seeds'.format(time.time() - itime,
          len(seeds)))

    itime = time.time()
    modules, sets1, sets2, scores = calc_modules(
            build_weighted_graph(G, modes['alpha'], Pg), seeds, Wr, params=modes)
    print('bi-modules took {:.2f} secs'.format(time.time() - itime))

    print('saving to: {}'.format(outfile))
    datafiles = {'infile': infile, 'outfile': outfile,
                 'empiric_null': empiric_null}
    sio.savemat(outfile,
                {'submodules1': sets1, 'submodules2': sets2,
                 'modes': modes, 'datafiles': datafiles})


def run_mono(infile, modes=deepcopy(DefModes), empiric_null=None, seedfile=None):
    """ detect dense subgraphs.
        INPUT:
            [infile]: matlab datafile containing a graph [G]
            [modes]: optional, a dictionary of program parameters
            [empiric_null]: optional, matlab datafile containing matrix [P] of
                    probabilties for observing an edge between nodes (i, j)
            [seedfile]: for expanding existing SCOMs
        OUTPUT: (written to file)
            [modules]: cell array containing subsets of indices
                (NOTE!: matlab indices starting from 1)
            [scores]: log odds ratio scores for the modules
            [modes]: program parameters used
            [datafiles]: data input to program
    """
    print('detecting mono-SCOMs')
    outfile = '%s SCOM' % infile

    G = sio.loadmat(infile)['G']

    if not empiric_null:
        P = None
        empiric_null = np.nan
    else:
        print('\nloading {}'.format(empiric_null))
        P = sio.loadmat(empiric_null)['P']
        epsilon = 1e-4
        P = np.maximum(epsilon, np.minimum(1-epsilon, P))

    if seedfile is None:
        # generate seeds using Charikar
        Wg = build_weighted_graph(G, modes['alpha'], P)
        itime = time.time()
        seeds, scores = calc_charikar_all_subgraphs(Wg, modes)
        seeds = [s for s, c in zip(seeds, scores) if c >= modes['min_seed']]
        seedfile = np.nan
        print('seeds took {:.2f} secs, {} seeds'.format(time.time() - itime,
              len(seeds)))
    else:
        # use existing seeds for expSCOM
        seeds = sio.loadmat(seedfile)['modules'][0].tolist()
        seeds = [(s[0] - 1).tolist() for s in seeds]  # to 0-indexing
        # override modes
        modes['remove'] = False
        modes['cleanup'] = False
        modes['max_size'] = 100
        modes['max_overlap'] = 1
        outfile = '%s expSCOM' % infile

    itime = time.time()
    modules, sets1, sets2, scores = calc_modules(
            build_weighted_graph(G, modes['alpha'], P), seeds, params=modes)
    print('mono-modules took {:.2f} secs'.format(time.time() - itime))

    print('saving to: {}'.format(outfile))
    datafiles = {'infile': infile, 'outfile': outfile,
                 'empiric_null': empiric_null, 'seedfile': seedfile}
    sio.savemat(outfile,
                {'modules': sets1, 'scores': scores,
                 'modes': modes, 'datafiles': datafiles})


def main(argv):
    if argv[1] == 'div':
        run_mono(argv[2])
    elif argv[1] == 'con':
        run_mono(argv[2])
    elif argv[1] == 'sep':
        run_bimodules(sys.argv[2])
    elif argv[1] == 'exp':
        run_mono(argv[2], seedfile=argv[3])
    else:
        print('unknown option')


if __name__ == '__main__':
#    run_bimodules('example/conserv net 25')
    main(sys.argv)
