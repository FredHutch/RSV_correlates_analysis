#!/bin/bash
ml fhR/4.0.2-foss-2019b
ml jbigkit
sbatch -c10 --mem 10G --time=7-0 --array=1-32 run_sl_assays.sh
