#!/bin/bash

sbatch 	-n 10 \
	-N 1 \
	--mem=100G \
	--wrap="cp -r /datacommons/wraylab/Alejo_Files/alejo/singlecell/atlas_metamorphosis/ /work/cjm124/scRNAanalysis/"


