"""
Download the data in CLDF format and convert the datasts.
"""

from git import Repo
import json
from pathlib import Path

def download(datasets, pth):
    """
    Download all datasets as indicated with GIT.
    """
    for dataset, conditions in datasets.items():
        if pth.joinpath(dataset, "cldf", "cldf-metadata.json").exists():
            print("[i] skipping existing dataset {0}".format(dataset))
        else:
            repo = Repo.clone_from(
                    "https://github.com/"+conditions["path"]+".git",
                    pth / dataset)
            repo.git.checkout(conditions["version"])
            print("[i] downloaded {0}".format(dataset))


if __name__ == "__main__":
    with open("data.json") as f:
        data = json.load(f)

    download(data, Path("cldf-datasets"))

