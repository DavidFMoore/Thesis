#! /bin/bash

for chan in '20_50'; do #'50_80' '80_110' '110_140' '140_170'; do
  #./scripts/sim_pspec.py -c ${chan} -r 1000  -C psa_null scripts/6C.cf
  #./scripts/sim_pspec.py -c ${chan} -r 2000  -C psa_null scripts/6C.cf
  #./scripts/sim_pspec.py -c ${chan} -r 5000  -C psa_null scripts/6C.cf
  ./scripts/sim_pspec.py -c ${chan} -r 10000 -C psa_null scripts/VLSS.cf
done
