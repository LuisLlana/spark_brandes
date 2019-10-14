This project contains a parallelization of the Brandes Algorithm using Spark.
There are 3 implementations of the Brandes Algorithm.
- betweenness_centrality_ori.py: this file contains the original
  algorithm from the networtx python module.

- betweenness_centrality_mr.py: This module contains a sequential
  algorithm based on the original one but using the map/reduce
  functions of python.

- betweenness_centrality_spark.py: This file contains the
  implementation using Spark.


There also an implementation of the Girvan Newman algorithm using the
previous implementation of the brandes algorithm. This implementations
is based on the one provided by Kazem Jahanbakhsh
(https://github.com/kjahan).


The datasets have been taken from the Stanford Large Network Dataset
Collection (https://snap.stanford.edu/data/)

The program to meassure the Brandes Algorithm is
benchmark_brandes.py
The constant DATASET_DIR should be rewritten to adjust the directory
containing the datasets.
