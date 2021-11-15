from lingpy.compare.partial import Partial

def run(wordlist):
    part = Partial(wordlist)
    part.get_partial_scorer(runs=10000)
    part.partial_cluster(
            cluster_method="infomap", ref="cogids",
            method="lexstat", threshold=0.55)
    D = {0: ["doculect", "concept", "value", "form", "tokens", "cogids"]}
    for idx in part:
        D[idx] = [part[idx, h] for h in D[0]]
    return D
