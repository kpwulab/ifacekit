#!/bin/tcsh
# Run Phenix MolProbity on each model, keep the outlier lists, delete the huge probe file,
# then compare with ifacekit. Run it from the project folder in a NEW terminal window
# (no CCP4 setup, no .venv), so that Phenix uses its own Python.
# usage: tcsh tools/run_molprobity.csh data/9lqk.cif data/7d47.cif
set PH = ~/Applications/phenix-2.0-5936/bin
set ROOT = `pwd`
foreach f ($argv)
  set name = `basename $f | sed 's/\.[^.]*$//'`
  mkdir -p $ROOT/mp_$name
  cd $ROOT/mp_$name
  echo "== MolProbity $name"
  $PH/phenix.molprobity $ROOT/$f > mp.log
  rm -f molprobity_probe.txt
  tail -12 mp.log
  cd $ROOT
end
rm -f mp_5hpt/molprobity_probe.txt
echo "done. Now open another terminal with the venv (source .venv/bin/activate.csh) and run tools/compare_all.csh"
