"""
Prepare test and training data for the experiments.
"""
from lingrexpi.reconstruct import Trainer
from pathlib import Path
import json
from tabulate import tabulate
from tqdm import tqdm as progressbar

datasets = {
        "wangbai": ("ProtoBai", "Bai"),
        "hillburmish": ("ProtoBurmish", "Burmish"),
        "ltkkaren": ("ProtoKaren", "Karen"),
        "yanglalo": ("ProtoLalo", "Lalo"),
        "carvalhopurus": ("ProtoPurus", "Purus"),
        "meloniromance": ("Latin", "Romance"),
        }
table = []
for ds, (proto, name) in progressbar(datasets.items()):
    trn = Trainer(
            str(Path("data", ds+".tsv")),
            ref="cogids",
            fuzzy=True,
            target=proto
            )
    cognates, words = 0, 0
    etd = trn.get_etymdict(ref="cogids")
    for cogid, idxs in etd.items():
        lngs = [trn[idx[0], "doculect"] for idx in idxs if idx]
        if proto in lngs and len(lngs) > 2:
            cognates += 1
            words += len(lngs)
    table += [[name, trn.width, cognates, words]]

    for i in range(100):
        wl, test_set, _ = trn.split(proportions=(90, 10, 0))
        wl.output(
                "tsv", 
                filename=str(Path(
                    "results", "testlists", name, "test-{0}".format(i+1))), 
                ignore="all", 
                prettify=False
                )
        with open(Path(
                "results", "testitems", name, "test-{0}.json".format(i+1)), 
                "w") as f:
            json.dump(test_set, f)
print(tabulate(table))
