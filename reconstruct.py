"""
Test the reconstruction method.
"""
from lingrexpi.reconstruct import (
        Trainer, PatternReconstructor, CorPaRClassifier, simple_align,
        clean_sound, eval_by_dist, eval_by_bcubes, alm2tok)
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import CategoricalNB
from pathlib import Path
from itertools import combinations
from lingpy import log
from statistics import mean
import json
from lingpy.align.pairwise import edit_dist


def align_ppf(x, y, z, **kw):
    return simple_align(x, y, z, align=True, position=True, prosody=True,
            firstlast=True, **kw)


def align_pp(x, y, z, **kw):
    return simple_align(x, y, z, align=True, position=True, prosody=True,
            firstlast=False, **kw)


def align_p(x, y, z, **kw):
    return simple_align(x, y, z, align=True, position=True, prosody=False,
            firstlast=False, **kw)


def align(x, y, z, **kw):
    return simple_align(x, y, z, align=True, position=False, prosody=False,
            firstlast=False, **kw)



datasets = [
        ("ltkkaren", "ProtoKaren"),
        ("yanglalo", "ProtoLalo"),
        ("wangbai", "ProtoBai"),
        ("carvalhopurus", "ProtoPurus")
        ]


methods = [
        ("APosProsFl", align_ppf),
        ("APosPros", align_pp),
        ("APos", align_p),
        ("A", align)
        ]


classifiers = [
        ("Categorical", CategoricalNB(), False),
        ("ExtremeTrees", ExtraTreesClassifier(), False),
        ("SVM", SVC(), True),
        ("CorPaR", CorPaRClassifier(), False)
        ]


BC, ED, NED = {}, {}, {}
RUNS = 10

for ds, proto in datasets:
    for clf_name, clf, onehot in classifiers:
        for meth_name, meth in methods:
            trn = Trainer(
                    str(Path("data", ds+".tsv")), 
                    ref="cogids",
                    fuzzy=True,
                    target=proto
                    )
            scores_bc, scores_ed, scores_ned = [], [], []
            for i in range(RUNS):
                wl, test_set, _ = trn.split(proportions=(80, 20, 0))
                wl.fit(clf=clf, func=meth, onehot=onehot)
                results = []
                for cogid, target, alignment, languages in test_set:
                    if len(alignment) > 1:
                        results += [[
                            alm2tok(
                                wl.predict(
                                    [alm2tok(a) for a in alignment],
                                    languages
                                    )
                                ),
                            alm2tok(target)]]
                scores_bc += [eval_by_bcubes(results)]
                scores_ed += [
                        eval_by_dist(results) 
                            ]

                scores_ned += [
                        eval_by_dist(
                            results, 
                            func=lambda x, y: edit_dist(x, y, normalized=True)
                            )
                            ]
                print("[i] {0} {1} {2} {3} {4:.2f} {5:.2f} {6:.2f}".format(
                    i+1,
                    ds,
                    clf_name,
                    meth_name,
                    scores_bc[-1],
                    scores_ed[-1],
                    scores_ned[-1]))

            BC[ds, clf_name, meth_name] = mean(scores_bc)
            ED[ds, clf_name, meth_name] = mean(scores_ed)
            NED[ds, clf_name, meth_name] = mean(scores_ned)
            print("[ii] {0} {1} {2} {3:.2f} {4:.2f} {5:.2f}".format(
                ds, clf_name, meth_name, mean(scores_bc), mean(scores_ed),
                mean(scores_ned)
                ))
with open("results.json", "w") as f:
    f.write(json.dumps(
        {
            "bcubes": BC, "ed": ED, "ned": NED,
            "datasets": data,
            }, indent=2))

