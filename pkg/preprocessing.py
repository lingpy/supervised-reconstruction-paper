from lingpy.compare.partial import Partial

def run(wordlist):
    part = Partial(wordlist)
    part.partial_cluster(cluster_method="infomap", ref="cogids", method="sca", threshold=0.45)
    D = {0: ["doculect", "concept", "value", "form", "tokens", "cogids"]}
    for idx in part:
        D[idx] = [part[idx, h] for h in D[0]]
    return D
