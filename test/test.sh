#!/bin/bash
dir=${0%/*}
if [ $# -gt 0 ]; then
 tests=$*;
else
 tests=$(echo {00..19});
fi
for i in ${tests}
do
echo
echo "Test "$i":"
${dir}/../test$i/test$i
done
