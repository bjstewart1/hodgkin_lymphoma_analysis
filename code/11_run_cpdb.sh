#!/bin/bash
cd /lustre/scratch117/cellgen/team297/bs16/current_projects/HLN_project/data
mkdir HLN_cpdb
cellphonedb method statistical_analysis cpdb_meta.txt cpdb_adata.h5ad --iterations=1000 --threads=10 --output-path=HLN_cpdb --result-precision=3 --threshold=0.1 --verbose