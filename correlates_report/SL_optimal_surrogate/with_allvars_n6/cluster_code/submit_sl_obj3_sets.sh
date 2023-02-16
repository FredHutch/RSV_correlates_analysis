#!/bin/bash
ml fhR/4.0.2-foss-2019b
ml jbigkit
# sbatch -c10 --mem 10G --time=7-0 --array=1-32 run_sl_assays.sh
sbatch -c10 --mem 10G --time=7-0 --array=1-29 run_sl_obj3_sets.sh
# sbatch -c10 --mem 33G --time=7-0 --array=2 run_sl_obj3_sets.sh