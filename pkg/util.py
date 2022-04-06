"""
Module provides methods for linguistic reconstruction.
"""
from lingpy.align.sca import Alignments, get_consensus
from lingpy.sequence.sound_classes import prosodic_string, class2tokens
from lingpy.align.multiple import Multiple
from lingpy.align.pairwise import edit_dist, nw_align
from lingpy.evaluate.acd import _get_bcubed_score as get_bcubed_score
from collections import defaultdict
from lingpy.align.sca import normalize_alignment
import random
import networkx as nx
from networkx.algorithms.clique import find_cliques
from itertools import combinations
from itertools import permutations


def ungap(alignment, languages, proto):
    cols = []
    pidxs = []
    for i, taxon in enumerate(languages):
        if taxon == proto:
            pidxs += [i]
    merges = []
    for i in range(len(alignment[0])):
        col = [row[i] for row in alignment]
        col_rest = [site for j, site in enumerate(col) if j not in pidxs]
        if "-" in col_rest and len(set(col_rest)) == 1:
            merges += [i]
    if merges:
        new_alms = []
        for i, row in enumerate(alignment):
            new_alm = []
            mergeit = False
            started = True
            for j, cell in enumerate(row):
                if j in merges or mergeit:
                    mergeit = False
                    if not started: #j != 0:
                        if cell == "-":
                            pass
                        else:
                            if not new_alm[-1]:
                                new_alm[-1] += cell
                            else:
                                new_alm[-1] += '.'+cell
                    else:
                        mergeit = True
                        if cell == "-":
                            new_alm += [""]
                        else:
                            new_alm += [cell]
                else:
                    started = False
                    new_alm += [cell]
            for k, cell in enumerate(new_alm):
                if not cell:
                    new_alm[k] = "-"
            new_alms += [new_alm]
        return new_alms
    return alignment


def clean_sound(sound):
    itms = sound.split('.')
    return ".".join([s.split('/')[1] if "/" in s else s for s in itms])


def alm2tok(seq, gap="-"):
    """Turn an alignent into a sequence."""
    return [clean_sound(x) for x in unsegment(seq) if x != gap]


def unsegment(seq):
    out = []
    for itm in seq:
        out += itm.split('.')
    return out


def missing_items(length):
    out = []
    for i in range(1, length):
        tmp = i * [1] + (length-i) * [0]
        for tup in permutations(tmp):
            if tup not in out:
                out += [tup]
    return out


class CorPaRClassifier(object):

    def __init__(self, missing=0, threshold=1):
        self.G = nx.Graph()
        self.missing = 0
        self.threshold = threshold

    def compatible(self, ptA, ptB):
        match, mismatch = 0, 0
        for a, b in zip(ptA, ptB):
            if not a or not b:
                pass
            elif a == b:
                match += 1
            else:
                mismatch += 1
        return match, mismatch

    def consensus(self, nodes):
        
        cons = []
        for i in range(len(nodes[0])):
            nocons = True
            for node in nodes:
                if node[i] != self.missing:
                    cons += [node[i]]
                    nocons = False
                    break
            if nocons:
                cons += [self.missing]
        return tuple(cons)

    def fit(self, X, y):
        """
        TODO: handle gaps sufficiently!
        TODO: think about post-processing procedure for candidates!
        """
        # get identical patterns
        P = defaultdict(list)

        for i, row in enumerate(X):
            P[tuple(row + [y[i]])] += [i]
        for k, v in P.items():
            self.G.add_node(k, freq=len(v))
        
        # make graph
        for (pA, vA), (pB, vB) in combinations(P.items(), r=2):
            match, mismatch = self.compatible(pA, pB)
            if not mismatch and match >= self.threshold:
                self.G.add_edge(pA, pB, weight=match)
                
        self.patterns = defaultdict(lambda : defaultdict(list))
        self.lookup = defaultdict(lambda : defaultdict(int))
        # get cliques
        for nodes in find_cliques(self.G):
            cons = self.consensus(list(nodes))
            self.patterns[cons[:-1]][cons[-1]] = len(nodes)
            for node in nodes:
                self.lookup[node[:-1]][cons[:-1]] += len(nodes)
        self.candidates = {}
        self.predictions = {}
        for ptn in self.patterns:
            self.predictions[ptn] = [x for x, y in sorted(
                self.patterns[ptn].items(),
                key=lambda p: p[1],
                reverse=True)][0]
        for ptn in self.lookup:
            ptnB = [x for x, y in sorted(self.lookup[ptn].items(),
                key=lambda p: p[1],
                reverse=True)][0]
            self.predictions[ptn] = self.predictions[ptnB]

        # make index of data points for quick search based on attested data
        self.ptnlkp = defaultdict(list)
        for ptn in self.patterns:
            for i in range(len(ptn)):
                if ptn[i] != self.missing:
                    self.ptnlkp[i, ptn[i]] += [ptn]

        
    def predict(self, matrix):
        out = []
        for row in matrix:
            ptn = tuple(row)
            try:
                out += [self.predictions[ptn]]
            except KeyError:
                candidates = []
                visited = set()
                for i in range(len(ptn)-1):
                    if ptn[i] != self.missing:
                        for ptnB in self.ptnlkp[i, ptn[i]]:
                            if ptnB not in visited:
                                visited.add(ptnB)
                                match, mismatch = self.compatible(ptn, ptnB)
                                if match and not mismatch:
                                    candidates += [(ptnB, match+len(ptn))]
                                elif match-mismatch:
                                    candidates += [(ptnB, match-mismatch)]
                if candidates:
                    ptn = [x for x, y in sorted(
                        candidates,
                        key=lambda p: p[1],
                        reverse=True)][0]
                    self.predictions[tuple(row)] = self.predictions[ptn]
                    out += [self.predictions[tuple(row)]]
                else:
                    out += [self.missing]
        return out
        


class ReconstructionBase(Alignments):
    """
    """
    def __init__(
            self, infile, target=None, ref="cogids", fuzzy=True, 
            transcription="form", missing="Ø", gap="-"):
        Alignments.__init__(self, infile, fuzzy=fuzzy, ref=ref,
                transcription=transcription)
        self.target = target
        self.missing = missing
        self.gap = gap
        self.languages = [t for t in self.cols if t != target]
        self.target = target
        self.tgtidx = self.cols.index(target)
        self.lngidx = {t: self.cols.index(t) for t in self.languages}

    def iter_sequences(self, valid_target=True, aligned=False):
        """
        Iterate over aligned or unaligned sequences with or without the target \
                sequence.
        """
        seq_ref = self._alignments if aligned else self._segments
        for cogid, idxs in self.etd[self._ref].items():
            if valid_target and idxs[self.tgtidx]:
                if self._mode == "fuzzy":
                    target = self[
                            idxs[self.tgtidx][0],
                            seq_ref
                            ].n[
                                    self[
                                        idxs[self.tgtidx][0], self._ref
                                        ].index(cogid)]
                else:
                    target = self[idxs[self.tgtidx][0], seq_ref]
                alignment, languages = [], []
                for j, lng in enumerate(self.languages):
                    lidx = self.lngidx[lng]
                    if idxs[lidx]:
                        languages += [lng]
                        idx = idxs[lidx][0]
                        if self._mode == "fuzzy":
                            alm = self[idx, seq_ref].n[
                                    self[idx, self._ref].index(cogid)]
                        else:
                            alm = self[idx, seq_ref]
                        alignment.append([clean_sound(x) for x in alm])
                alignment.append([clean_sound(x) for x in target])
                if aligned:
                    alignment = normalize_alignment(alignment)
                languages.append(self.target)
                yield cogid, alignment, languages
            elif not valid_target:
                alignment, languages = [], []
                for j, lng in enumerate(self.cols):
                    lidx = self.lngidx[lng]
                    if idxs[lidx]:
                        languages += [lng]
                        idx = idxs[lidx][0]
                        if self._mode == "fuzzy":
                            alm = self[idx, seq_ref].n[
                                    self[idx, self._ref].index(cogid)]
                        else:
                            alm = self[idx, seq_ref]
                        alignment.append(alm)
                if aligned: 
                    alignment = normalize_alignment(alignment)
                yield cogid, alignment, languages


class Trainer(ReconstructionBase):
    """
    Base class to split data into test and training sets.
    """

    def split(self, proportions=None, seed=None, minrefs=2):
        """
        Split into training, test, and evaluation set.

        @param minrefs: minimal number of reflexes when splitting the data. 
        """
        if seed:
            random.seed(seed)
        proportions = proportions or (80, 10, 10)
        assert sum(proportions) == 100
        
        sample = {}
        for cogid, alignment, languages in self.iter_sequences(
                valid_target=True):
            if alignment[:-1] and len(languages) > minrefs:
                sample[cogid] = (alignment[:-1], alignment[-1], languages[:-1])
        cogids = list(sample)

        trn_set = random.sample(
                cogids, int(
                    proportions[0]/100 * len(cogids)+0.5))
        rest = [cogid for cogid in cogids if cogid not in trn_set]
        tst_set = random.sample(
                rest, 
                int(proportions[1]/(proportions[1]+proportions[2])*len(rest)+0.5)
                )
        dev_set = [cogid for cogid in rest if cogid not in tst_set]

        # make new word lists
        D = {0: ["doculect", "concept", "form", "tokens", "alignment", "cogid"]}
        idx = 1
        for cogid in trn_set:
            alms, trg, lngs = sample[cogid]
            for alm, lng in zip(alms, lngs):
                tks = [x for x in alm if x != "-"]
                D[idx] = [lng, str(cogid), "".join(tks), tks, alm, cogid]
                idx += 1
            tks = [x for x in trg if x != "-"]
            D[idx] = [self.target, str(cogid), "".join(tks), tks, trg, cogid]
            idx += 1
        wl_trn = PatternReconstructor(D, target=self.target, ref="cogid", fuzzy=False)

        lst_tst = []
        for cogid in tst_set:
            alms, trg, lngs = sample[cogid]
            lst_tst += [[cogid, trg, alms, lngs]]
        lst_dev = []
        for cogid in dev_set:
            alms, trg, lngs = sample[cogid]
            lst_dev += [[cogid, trg, alms, lngs]]

        return wl_trn, lst_tst, lst_dev


class OneHot(object):

    def __init__(self, matrix):
        self.vals = []
        for i in range(len(matrix[0])):
            cols = [row[i] for row in matrix]
            self.vals += [sorted(set(cols))+["?"]]
    
    def __call__(self, matrix):
        out = [[] for row in matrix]
        for i, vals in enumerate(self.vals):
            for j in range(len(matrix)):
                template = [0 for k in vals]
                try:
                    template[matrix[j][i]] = 1
                except IndexError:
                    template[-1] = 1
                out[j] += template
        return out



def transform_alignment(
        seqs, 
        languages, 
        all_languages,
        align=True,
        training=True,
        missing="Ø", 
        gap="-",
        startend=True,
        prosody=False,
        position=False,
        firstlast=False,
        ):
    if align:
        seqs = [[s for s in seq if s != gap] for seq in seqs]
        msa = Multiple([[s for s in seq if s != gap] for seq in seqs])
        msa.prog_align()
        alms = [alm for alm in msa.alm_matrix]
    else:
        seqs = [[s for s in seq if s != gap] for seq in seqs]
        alms = normalize_alignment([s for s in seqs])
    if training:
        alms = ungap(alms, languages, languages[-1])
        these_seqs = seqs[:-1]
    else:
        these_seqs = seqs
    matrix = [[missing for x in all_languages] for y in alms[0]]
    for i in range(len(alms[0])):
        for j, lng in enumerate(languages):
            lidx = all_languages.index(lng)
            snd = clean_sound(alms[j][i])
            matrix[i][lidx] = snd
    if position:
        for i in range(len(matrix)):
            matrix[i] += [i+1]
    if startend:
        matrix[0] += [1]
        for i in range(1, len(matrix)-1):
            matrix[i] += [2]
        if len(matrix) > 1:
            matrix[-1] += [3]
    if prosody:
        for i, c in enumerate(
                get_consensus(
                    [class2tokens(
                        prosodic_string(
                            seqs[j], 
                            _output="CcV"
                            ), 
                        alms[j]
                        ) for j in range(len(these_seqs))],
                    gaps=True)):
            matrix[i] += [c]
    if firstlast:
        if training:
            all_seqs = len(all_languages)-1
        else:
            all_seqs = len(all_languages)
        for i, row in enumerate(matrix):
            for j in range(all_seqs):
                matrix[i] += [matrix[0][j], matrix[-1][j]]
    
    # for debugging
    for row in matrix:
        assert len(row) == len(matrix[0])
    return matrix


class PatternReconstructor(ReconstructionBase):
    """
    Automatic reconstruction with correspondence patterns.
    """

    def fit(self, clf=None, onehot=False, func=None, aligned=False):
        """
        Fit a classifier to the data.

        @param clf: a classifier with a predict function.
        """
        self.patterns = defaultdict(lambda : defaultdict(list))
        self.occurrences = defaultdict(list)
        self.func = func or transform_alignment
        
        for cogid, alignment, languages in self.iter_sequences(
                valid_target=True, aligned=aligned):
            if len(alignment) >= 2:
                matrix = self.func(
                        alignment,
                        languages,
                        self.languages+[self.target],
                        training=True)
                for i, row in enumerate(matrix):
                    ptn = tuple(row[:len(self.languages)]+row[len(self.languages)+1:])
                    self.patterns[ptn][row[len(self.languages)]] += [
                        (cogid, i)]
                    for j, lng in enumerate(self.languages):
                        if row[j] not in [self.missing]:
                            self.occurrences[lng, j, row[j]] += [(cogid, i)]
                    for j in range(len(self.languages)+1, len(row)):
                        self.occurrences["feature-{0}".format(j-1), j-1, row[j]] += [(cogid, i)]
        
        self.snd2idx = {(i, self.missing): 0 for i in
                range(len(matrix[0]))}
        for i in range(len(matrix[0])):
            self.snd2idx[i, self.gap] = 1

        idxtracker = {i: 2 for i in range(len(matrix[0]))}
        for lng, lidx, sound in self.occurrences:
            last_idx = idxtracker[lidx]
            if (lidx, sound) not in self.snd2idx:
                self.snd2idx[lidx, sound] = last_idx
                idxtracker[lidx] += 1

        self.tgt2idx = {}
        idx = 1
        for pattern in self.patterns:
            for sound in self.patterns[pattern]:
                if sound not in self.tgt2idx:
                    self.tgt2idx[sound] = idx
                    idx += 1

        self.matrix = []
        self.solutions = []
        for pattern, sounds in self.patterns.items():
            for sound, vals in sounds.items():
                freq = len(vals)
                tidx = self.tgt2idx[sound]
                row = []
                for i in range(len(pattern)):
                    sidx = self.snd2idx[i, pattern[i]]
                    row += [sidx]
                for cogid, idx in vals:
                    self.matrix += [row]
                    self.solutions += [tidx]
        self.dim = len(self.matrix[0])
        if clf is not None:
            self.clf = clf 
        else:
            self.clf = CorPaRClassifier() 
        if onehot:
            self.onehot = OneHot(self.matrix)
            self.clf.fit(self.onehot(self.matrix), self.solutions)
        else:
            self.clf.fit(self.matrix, self.solutions)
        self.idx2tgt = {v: k for k, v in self.tgt2idx.items()}

    def predict(
            self, alignment, languages, unknown="?", onehot=False,
            desegment=True):
        """
        Predict a target word from an alignment.

        @param desegment: Return the form without gaps and ungapped tokens.
        """
        matrix = self.func(alignment, languages, self.languages, training=False)
        for row in matrix:
            if len(row) != self.dim:
                print(alignment, languages, matrix)
                input("predict")
        #alignment = normalize_alignment(alignment)
        new_matrix = [[0 for char in row] for row in matrix]
        for i, row in enumerate(matrix):
            for j, char in enumerate(row):
                new_matrix[i][j] = self.snd2idx.get((j, char), 0)
        if hasattr(self, "onehot"):
            new_matrix = self.onehot(new_matrix)
        out = [self.idx2tgt.get(idx, unknown) for idx in self.clf.predict(new_matrix)]
        return alm2tok(out) if desegment else out


def eval_by_dist(data, func=None, **kw):
    """
    Evaluate by measuring distances between sequences.
    
    @note: Defaults to the unnormalized edit distance.
    """
    func = func or edit_dist
    scores = []
    for seqA, seqB in data:
        if not seqA:
            seqA = ["?"]
        if not seqB:
            seqB = ["?"]
        scores += [func(seqA, seqB, **kw)]
    return sum(scores)/len(scores)


def eval_by_bcubes(data, func=None, **kw):
    numsA, numsB = {"": 0}, {"": 0}
    func = nw_align
    almsA, almsB = [], []
    for seqA, seqB in data:
        if not seqA:
            seqA = ["?"]
        if not seqB:
            seqB = ["?"]
        almA, almB, score = func(seqA, seqB, **kw)
        for a, b in zip(almA, almB):
            if not a in numsA:
                numsA[a] = max(numsA.values())+1
            if not b in numsB:
                numsB[b] = max(numsB.values())+1
            almsA += [numsA[a]]
            almsB += [numsB[b]]
    p, r = get_bcubed_score(almsA, almsB), get_bcubed_score(almsB, almsA)
    return 2*(p*r)/(p+r)
                


