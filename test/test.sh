#!/bin/bash
dir=${0%/*}
if [ $# -gt 0 ]; then
 tests=$*;
else
 tests='01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16';
fi
for i in ${tests}
do
echo
echo "Test "$i":"
${dir}/../test$i/test$i
done
