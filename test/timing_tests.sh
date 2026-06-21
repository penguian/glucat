#!/bin/bash
dir=${0%/*}
sp=${1:-8}
ss=${2:-11}
sg=${3:-11}
st=${4:-8}
se=${5:-8}
sv=${6:-16}

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
echo
echo "expressions "${se}":"
${dir}/../expressions/expressions ${se} | tee expressions-${se}.out
echo
echo "versor "${sv}":"
${dir}/../versor/versor ${sv} | tee versor-${sv}.out
