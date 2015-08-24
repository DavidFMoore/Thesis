#! /bin/bash

for chan in '20_50' '50_80' '80_110' '110_140' '140_170'; do
  ./scripts/sim_correlated_pspec.py -c ${chan} -C psa_null scripts/6C.cf
done
