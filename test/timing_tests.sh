#!/bin/bash
dir=${0%/*}
if [ $# -gt 3 ]; then
 sp=$1
 ss=$2
 sg=$3
 st=$4
else
 sp=8
 ss=11
 sg=11
 st=8
fi

echo "products "${sp}":"
${dir}/../products/products ${sp} | tee products-${sp}.out
echo
echo "squaring "${ss}":"
${dir}/../squaring/squaring ${ss} | tee squaring-${ss}.out
echo
echo "gfft_test "${sg}":"
${dir}/../gfft_test/gfft_test ${sg} | tee gfft_test-${sg}.out
echo
echo "transforms "${st}":"
${dir}/../transforms/transforms ${st} | tee transforms-${st}.out
