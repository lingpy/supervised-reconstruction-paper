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

```
\begin{tabular}{llrrr}
\hline
 Classifier   & Analysis   &     ED &    NED &     BC \\
\hline
 svm          & PosStrIni  & 0.7491 & 0.1598 & 0.8110 \\
 svm          & PosStr     & 0.7478 & 0.1594 & 0.8115 \\
 svm          & PosIni     & 0.7701 & 0.1624 & 0.8077 \\
 svm          & StrIni     & 0.7578 & 0.1601 & 0.8110 \\
 svm          & Pos        & 0.7685 & 0.1618 & 0.8084 \\
 svm          & Str        & 0.7681 & 0.1614 & 0.8086 \\
 svm          & Ini        & 0.7895 & 0.1641 & 0.8061 \\
 svm          & none       & 0.8059 & 0.1673 & 0.8006 \\
 corpar       & PosStrIni  & 0.8503 & 0.1816 & 0.7862 \\
 corpar       & PosStr     & 0.8655 & 0.1826 & 0.7854 \\
 corpar       & PosIni     & 0.8425 & 0.1802 & 0.7882 \\
 corpar       & StrIni     & 0.8402 & 0.1771 & 0.7924 \\
 corpar       & Pos        & 0.8836 & 0.1847 & 0.7840 \\
 corpar       & Str        & 0.9048 & 0.1851 & 0.7848 \\
 corpar       & Ini        & 0.8342 & 0.1763 & 0.7946 \\
 corpar       & none       & 0.9379 & 0.1898 & 0.7821 \\
\hline
\end{tabular}
[i] svm
\begin{tabular}{lrrrrr}
\hline
 DATASET   &   PosStrIni &   StrIni &    Str &    Ini &   none \\
\hline
 Bai       &      0.7848 &   0.7870 & 0.7832 & 0.7846 & 0.7770 \\
 Burmish   &      0.8388 &   0.8418 & 0.8420 & 0.8405 & 0.8226 \\
 Karen     &      0.8696 &   0.8736 & 0.8734 & 0.8731 & 0.8723 \\
 Lalo      &      0.7232 &   0.7214 & 0.7204 & 0.7202 & 0.7191 \\
 Purus     &      0.9011 &   0.9021 & 0.9016 & 0.9013 & 0.9022 \\
 Romance   &      0.7487 &   0.7401 & 0.7310 & 0.7171 & 0.7103 \\
\hline
\end{tabular}
[i] corpar
\begin{tabular}{lrrrrr}
\hline
 DATASET   &   PosStrIni &   StrIni &    Str &    Ini &   none \\
\hline
 Bai       &      0.7485 &   0.7581 & 0.7560 & 0.7572 & 0.7560 \\
 Burmish   &      0.8319 &   0.8449 & 0.8422 & 0.8458 & 0.8331 \\
 Karen     &      0.8564 &   0.8581 & 0.8614 & 0.8604 & 0.8581 \\
 Lalo      &      0.6852 &   0.6874 & 0.6890 & 0.6893 & 0.6871 \\
 Purus     &      0.8688 &   0.8865 & 0.8730 & 0.8897 & 0.8880 \\
 Romance   &      0.7262 &   0.7192 & 0.6871 & 0.7253 & 0.6705 \\
\hline
\end{tabular}

```
