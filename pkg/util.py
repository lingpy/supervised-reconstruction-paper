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
from lingrex.reconstruct import ReconstructionBase, PatternReconstructor


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
        for cogid, alignment, languages in self.iter_sequences():
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

