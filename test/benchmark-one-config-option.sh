#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   benchmark-one-config-option.sh : Test one defined configuration.
#
#   begin                : Sun 2026-05-31
#   copyright            : (C) 2016-2026 by Paul C. Leopardi
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   For licence details, see ${package_dir}/COPYING

# set -x

here=$(cd $(dirname ${0})/;pwd)
config_options_file="${here}/benchmark-config-options.txt"
. ${here}/define-config-options.sh

line=$1
shift
args=$*

abbrev=$(eval echo $(match_option_pattern ${line} 1))
options=$(eval echo $(match_option_pattern ${line} 2))

pushd ${package_dir}/.. \
  > /dev/null

  rm -rf ${package_dir}.${line}
  cp -a ${package_dir} ${package_dir}.${line}

  pushd ${package_dir}.${line} \
    > /dev/null

    rm -rf benchmarks_runtime
    mkdir benchmarks_runtime
    ./configure ${options} \
      > /dev/null
    cp -a config.log benchmarks_runtime
    make clean \
      > /dev/null
    make ${args} \
      > benchmarks_runtime/make.log
    test/ldd_all.sh \
      > benchmarks_runtime/ldd.log

    # Run the actual benchmark in a subshell
    (
      source benchmarks/env-${abbrev}.sh
      cd benchmarks_runtime
      ../test/timing_tests.sh \
        > timing_tests.log
    )

    make clean \
      > /dev/null

  popd \
    > /dev/null

popd \
  > /dev/null
