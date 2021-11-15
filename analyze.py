"""
Analyze the data.
"""
from lingrexpi.reconstruct import (
        PatternReconstructor, CorPaRClassifier,
        transform_alignment)
from sklearn.svm import SVC
from sys import argv
from pathlib import Path
import json
from functools import partial

align_psf = partial(
        transform_alignment, align=True, position=True, prosody=True, startend=True)
align_ps = partial(
        transform_alignment, align=True, position=True, prosody=True, startend=False)
align_sf = partial(
        transform_alignment, align=True, position=False, prosody=True, startend=True)
align_pf = partial(
        transform_alignment, align=True, position=True, prosody=False, startend=True)
align_p = partial(
        transform_alignment, align=True, position=True, prosody=False, startend=False)
align_s = partial(
        transform_alignment, align=True, position=False, prosody=True, startend=False)
align_f = partial(
        transform_alignment, align=True, position=False, prosody=False, startend=True)
align = partial(
        transform_alignment, align=True, position=False, prosody=False, startend=False)


if len(argv) >= 2:
    if argv[1] == "svm":
        classifier = "svm"
        clf = lambda : SVC(kernel="linear")
        onehot = True
    elif argv[1] == "corpar":
        clf = CorPaRClassifier
        onehot = False
        classifier = "corpar"
RUNS = 100

datasets = [
        ("Burmish", "ProtoBurmish"),
        ("Purus", "ProtoPurus"),
        ("Karen", "ProtoKaren"),
        ("Bai", "ProtoBai"),
        ("Lalo", "ProtoLalo"),
        ("Romance", "Latin"),
        ]

methods = [

        ("PosStr", align_ps),
        ("PosIni", align_pf),
        ("StrIni", align_sf),
        ("Pos", align_p),
        ("Str", align_s),
        ("Ini", align_f),
        ("none", align),
        ("PosStrIni", align_psf),
        ]


for ds, proto in datasets:
    print("[i] analyzing {0}".format(ds))
    for i in range(RUNS):
        print("[i] analyzing {0} test {1}".format(ds, i+1))
        wlpath = str(Path(
            "results", "testlists", ds, "test-{0}.tsv".format(i+1)))
        with open(Path(
            "results", "testitems", ds, "test-{0}.json".format(i+1))) as f:
            tests = json.load(f)
        for meth_name, meth in methods:
            res_path = Path(
                "results", classifier, ds+"-"+meth_name+"-"+str(i+1)+".tsv")
            if not res_path.exists():
                pt = PatternReconstructor(
                        wlpath, ref="cogid", fuzzy=False, target=proto)
                pt.fit(clf=clf(), func=meth, onehot=onehot, aligned=False)
                results = []
                for cogid, target, alignment, languages in tests:
                    results += [[
                        pt.predict(
                            alignment,
                            languages),
                        target]]
                with open(res_path, "w") as f:
                    for a, b in results:
                        f.write(" ".join(a)+"\t"+" ".join(b)+"\n")

