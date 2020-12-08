#!/usr/bin/env python

import sys
import scanpy as sc

ad = sc.read(sys.argv[1])
ad.obs.to_csv(sys.argv[2], sep="\t")
