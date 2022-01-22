#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   fast-diff-one-config-output.sh : Find differences in test output for one configuration.
#
#   begin                : Sun 2016-04-03
#   copyright            : (C) 2016-2022 by Paul C. Leopardi
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

      diff -ub ${package_dir}/test_runtime/fast-test.configure.${abbrev}.out fast-test.configure.${abbrev}.out \
        > fast-test.configure.diff

      cat fast-test.configure.diff

    popd \
      > /dev/null
      
    pushd pyclical \
      > /dev/null
      if [ -f fast-test-check.out ]
      then
        diff -ub ${package_dir}/pyclical/test.out fast-test.check.out \
          > fast-test.out.diff

        cat fast-test.out.diff
      fi
    popd \
      > /dev/null

  popd \
    > /dev/null

popd \
  > /dev/null
