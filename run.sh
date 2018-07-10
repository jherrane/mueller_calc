#!/bin/bash
base="./mueller_calc"
geom_fold="Geometries"
T_fold="T"

# Sizes (in um)
# 0.0500 0.0598 0.0715 0.0855 0.1022 0.1223 0.1462 0.17480.2090 0.2500
declare -a shapes=("ell" "ob" "pro" "sph")
declare -a compositions=("P-aC")

BEGIN=1
END=15

for is in $(seq 0 3); do
   for i in $(seq $BEGIN $END); do
      mesh=" --mesh "$geom_fold"/${shapes[$is]}-P-aC-$i.h5"
      let "j = i + is*END"

      T=" -T $T_fold/T-$j.h5"
      out=" --out -$j"
      
      comm="$base$mesh$T$out"
      repl="sed -i '9s|.*|"$comm"|' script.sh"
      
      eval "$comm"       
   done
done
