#!/bin/tcsh
# usage (venv active, project folder): tcsh tools/compare_all.csh 5hpt 9lqk 7d47
foreach name ($argv)
  echo "=============== $name"
  python tools/compare_molprobity.py data/$name.cif mp_$name/molprobity_coot.py > results/${name}_vs_molprobity.txt
  cat results/${name}_vs_molprobity.txt
end
