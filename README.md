# Code and Data to Replicate the Study "A New Framework for Fast Automated Phonological Reconstruction Using Trimmed Alignments and Sound Correspondence Patterns"

## Data

The data is originally coded in CLDF (Forkel et al. 2018) and then converted to TSV formats needed as input for the framework. The conversion is done automatically with the help of the [pyedictor](https://github.com/lingpy/pyedictor) tool.  

To download the data, you can use this Python script:

```
$ python download.py
```

To prepare the data, we offer a Make file that uses `pyedictor` to convert the wordlists into their desired format:

```
$ make all
```


## Code

The code consists of three scripts. To install all dependencies needed for replication, we recommend to use `pip`:

```
$ pip install -r requirements.txt
```

To prepare the split of test and training data, the script `prepare.py` loads all six datasets and then creates 100 datesets in which only 90 percent of the original data (counted by the number of valid cognate sets) are retained.

```
$ python prepare.py
```

This will provide as output a table in LaTeX format:

```
\begin{tabular}{lrrr}
\hline
 Bai     & 10 &  459 &  3866 \\
 Burmish &  9 &  269 &  1711 \\
 Karen   & 11 &  365 &  3231 \\
 Lalo    &  8 & 1251 &  7815 \\
 Purus   &  4 &  199 &   693 \\
 Romance &  6 & 4147 & 18806 \\
\hline
\end{tabular}
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

