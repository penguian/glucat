#!/bin/bash
dir=${0%/*}
tests=$(echo {00..19})
for i in ${tests}
do
echo
echo "Test "$i":"
${dir}/../test$i/test$i $*
done
