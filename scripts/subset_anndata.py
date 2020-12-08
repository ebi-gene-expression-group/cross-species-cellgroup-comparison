#!/usr/bin/env python

import sys
import scanpy as sc
import pandas as pd

ad = sc.read(sys.argv[1])
obs = pd.read_csv(sys.argv[2], sep = "\\t", index_col = 0)
ad[obs.index].write(sys.argv[3])
