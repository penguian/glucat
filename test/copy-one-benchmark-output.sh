#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   copy-one-config-output.sh : Copy test output for one defined configuration.
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
dest=$2
abbrev=$(eval echo $(match_option_pattern ${line} 1))
pushd ${package_dir}/.. \
  > /dev/null

  pushd ${package_dir}.${line} \
    > /dev/null

    rm -rf ${dest}/${abbrev}
    cp -a benchmarks_runtime ${dest}/${abbrev}

  popd \
    > /dev/null

popd \
  > /dev/null
