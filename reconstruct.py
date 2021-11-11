"""
Test the reconstruction method.
"""
from lingrexpi.reconstruct import (
        Trainer, PatternReconstructor, CorPaRClassifier, transform_alignment,
        clean_sound, eval_by_dist, eval_by_bcubes)
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import CategoricalNB
from pathlib import Path
from itertools import combinations
from lingpy import log
from statistics import mean
import json
from lingpy.align.pairwise import edit_dist
from tqdm import tqdm as progressbar
from tabulate import tabulate
from functools import partial


align_psf = partial(
        transform_alignment, align=True, position=True, prosody=True,
        firstlast=True)
align_ps = partial(
        transform_alignment, align=True, position=True, prosody=True,
        firstlast=False)
align_sf = partial(
        transform_alignment, align=True, position=False, prosody=True,
        firstlast=True)
align_pf = partial(
        transform_alignment, align=True, position=True, prosody=False,
        firstlast=True)
align_p = partial(
        transform_alignment, align=True, position=True, prosody=False,
        firstlast=False)
align_s = partial(
        transform_alignment, align=True, position=False, prosody=False,
        firstlast=True)
align_f = partial(
        transform_alignment, align=True, position=True, prosody=True,
        firstlast=True)
align = partial(
        transform_alignment, align=True, position=False, prosody=False,
        firstlast=False)

datasets = [
        ("pharaocoracholaztecan", "ProtoCoracholNahua"),
        ("davletshinaztecan", "ProtoAztecan"),
        ("carvalhopurus", "ProtoPurus"),
        ("ltkkaren", "ProtoKaren"),
        ("wangbai", "ProtoBai"),
        ("yanglalo", "ProtoLalo"),
        ("meloniromance", "Latin"),
        ]


methods = [
        ("APSF", align_psf),
        ("APS", align_ps),
        #("APF", align_pf),
        #("ASF", align_sf),
        #("AP", align_p),
        #("AS", align_s),
        #("AF", align_f),
        ("A", align)
        ]


classifiers = [
        #("Categorical", CategoricalNB(), False),
        #("ExtremeTrees", ExtraTreesClassifier(), False),
        #("RandomForest", RandomForestClassifier(), False),
        ("SVM", lambda : SVC(kernel="linear"), True),
        #("MLP1", MLPClassifier(solver="lbfgs", max_iter=1000), True),
        ("CorPaR", CorPaRClassifier, False),
        #("MLP2", MLPClassifier(
        #    solver="lbfgs", hidden_layer_sizes=20,
        #    max_iter=1000), True),
        #("MLP3", MLPClassifier(
        #    solver="lbfgs", hidden_layer_sizes=(20, 10),
        #    max_iter=1000), True)
        ]


BC, ED, NED, S =  {}, {}, {}, {}
RUNS = 100
DEBUG = False
table = []
for ds, proto in datasets:
    for clf_name, clf, onehot in classifiers:
        for meth_name, meth in methods:
            trn = Trainer(
                    str(Path("data", ds+".tsv")), 
                    ref="cogids",
                    fuzzy=True,
                    target=proto
                    )
            scores_bc, scores_ed, scores_ned, sizes = [], [], [], []
            for i in progressbar(range(RUNS)):
                wl, test_set, _ = trn.split(proportions=(80, 20, 0))
                sizes.append([len(wl.etd["cogid"]), len(test_set)])
                wl.fit(clf=clf(), func=meth, onehot=onehot, aligned=False)
                results = []
                for cogid, target, alignment, languages in test_set:
                    assert len(alignment) > 1
                    results += [[
                        wl.predict(
                            alignment,
                            languages
                            ),
                        target]]
                    if DEBUG:
                        print(" ".join(results[-1][0]))
                        print(" ".join(results[-1][1]))
                        print("---")
                with open(
                        Path(
                            "results", 
                            "-".join([ds, clf_name, meth_name, str(i+1)])+"-res.tsv"
                            ), 
                        "w") as f:
                    for row in results:
                        f.write(" ".join(row[0])+"\t"+" ".join(row[1])+"\n")
                with open(
                        Path(
                            "results", 
                            "-".join([ds, clf_name, meth_name, str(i+1)])+"-tst.tsv"
                            ), 
                        "w") as f:
                    for cogid, target, alignment, languages in test_set:
                        for alm, lng in zip(alignment, languages):
                            f.write(str(cogid)+"\t"+lng+"\t"+" ".join(alm)+"\n")
                wl.output("tsv", filename=str(
                    Path("results", "-".join([ds, clf_name, meth_name,
                        str(i+1)])+"-wl")), ignore="all", prettify=False)
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
                if DEBUG:
                    print("[i] {0} {1} {2} {3} {4:.2f} {5:.2f} {6:.2f}".format(
                        i+1,
                        ds,
                        clf_name,
                        meth_name,
                        scores_bc[-1],
                        scores_ed[-1],
                        scores_ned[-1]))

            BC[ds+"-"+clf_name+"-"+meth_name] = scores_bc
            ED[ds+"-"+clf_name+"-"+meth_name] = scores_ed
            NED[ds+"-"+clf_name+"-"+meth_name] = scores_ned
            S[ds+"-"+clf_name+"-"+meth_name] = sizes 
            print("[i] {0} {1} {2} {3:.2f} {4:.2f} {5:.2f}".format(
                ds, clf_name, meth_name, mean(scores_bc), mean(scores_ed),
                mean(scores_ned)
                ))
            table += [[
                ds, clf_name, meth_name, 
                mean(scores_bc),
                mean(scores_ed),
                mean(scores_ned)
                ]]

with open(Path("results", "results.json"), "w") as f:
    f.write(json.dumps(
        {
            "bc": BC, "ed": ED, "ned": NED,
            "datasets": datasets,
            "table": table,
            "methods": [m[0] for m in methods],
            "classifiers": [c[0] for c in classifiers]
            }, indent=2))
print(tabulate(table, tablefmt="pipe"))
