# Code and Data to Replicate the Study "A New Framework for Fast Automated Phonological Reconstruction Using Trimmed Alignments and Sound Correspondence Patterns"

## Data

The data is originally coded in CLDF (Forkel et al. 2018) and then converted to TSV formats needed as input for the framework. The conversion is done automatically with the help of the [pyedictor](https://github.com/lingpy/pyedictor) tool. Two datasets are for now only offered in the form of a TSV file, since they were not published in this form before, but repositories in the form of CLDF are planned to be published later. 

If you are interested in how the original data was converted to the TSV files from the CLDF datasets, please check the shell-scripts in the folder `repos`.

## Code

The code consists of three scripts. To install all dependencies needed for replication, we recommend to use `pip`:

```
$ pip install -r requirements.txt
```

To prepare the split of test and training data, the script `prepare.py` loads all six datasets and then creates 100 datesets in which only 90 percent of the original data (counted by the number of valid cognate sets) are retained.

```
$ python prepare.py
```

To run the actual analysis for SVM classifiers, type:

```
$ python analyze.py svm
```

For the CorPaR classifier, type:

```
$ python analyze.py corpar
```

The results in the form of preditions are written to TSV files in the folders `results/svm` and `results/corpar` respectively. To analyze them by computing edit distance and B-Cubed F-scores, we use the script `evaluate.py`:

```
$ python evaluate.py
```

This will create three plots, two showing individual results, one showing combined results for individual datasets, as well as a list of general and individual results in the form of TSV files in the folder `results`. In addition, the major results will be printed to the terminal.

