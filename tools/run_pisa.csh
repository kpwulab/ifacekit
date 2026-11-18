#!/bin/tcsh
# Run PISA on every cif/pdb given and save the interface XML next to it.
# usage: tcsh tools/run_pisa.csh data/7d47.cif data/7yqm.cif ...
mkdir -p /tmp/ifacekit_pisa
foreach f ($argv)
  set name = `basename $f | sed 's/\.[^.]*$//'`
  set dir = `dirname $f`
  echo "== $name"
  pisa $name -analyse $f
  pisa $name -xml interfaces > $dir/${name}_pisa_interfaces.xml
end
