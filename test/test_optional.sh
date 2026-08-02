#!/bin/bash
dir=${0%/*}
tests=$(echo {00..19})
for i in ${tests}
do
echo
echo "Test "$i":"
if [ -x ./test$i ]; then
 ./test$i $*
elif [ -x ../test$i/test$i ]; then
 ../test$i/test$i $*
else
 ${dir}/../test$i/test$i $*
fi
done
