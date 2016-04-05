#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   diff-one-config-output.sh : Find differences in test output for one configuration.
#
#   begin                : Sun 2016-04-03
#   copyright            : (C) 2016 by Paul C. Leopardi
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   For licence details, see ${package_dir}/COPYING

# set -x

here=$(cd $(dirname ${0})/;pwd)
. ${here}/define-config-options.sh

line=$1

abbrev=$(eval echo $(match_option_pattern ${line} 1))

pushd ${package_dir}/.. \
  > /dev/null

  pushd ${package_dir}.${line} \
    > /dev/null

    pushd test_runtime \
      > /dev/null

      diff -ub ${package_dir}/test_runtime/test.configure.${abbrev}.out test.configure.${abbrev}.out \
        > test.configure.diff

      cat test.configure.diff

    popd \
      > /dev/null
      
    pushd pyclical \
      > /dev/null

      diff -ub ${package_dir}/pyclical/test.out test.out \
        > test.out.diff

      cat test.out.diff
      
    popd \
      > /dev/null

  popd \
    > /dev/null

popd \
  > /dev/null
