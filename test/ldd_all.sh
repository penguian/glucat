#!/bin/bash
for d in *; do
  f=$(basename $d);
  e="$d/$f";
  if [[ -f "$e" ]]; then
    ls "$e";
    ldd $e;
  fi;
done
