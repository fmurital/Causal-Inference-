#!/usr/bin/env bash
# Runs the whole project in order. Each step logs to logs/<script>.log
set -e
cd "$(dirname "$0")"
mkdir -p logs figures output
for s in 01_prepare_data 02_eda 03_main_analysis 04_sensitivity 04b_sensitivity_figures 06_draw_dag 07_paper_numbers; do
  echo "Running $s ..."; Rscript $s.R > logs/$s.log 2>&1
done
echo "Running 05_variance_check (slow) ..."; Rscript 05_variance_check.R > logs/05_variance_check.log 2>&1
echo ALL_DONE
