#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   copy-all-benchmark-outputs.sh : Copy benchmark output for all defined configurations.
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

dest=${1:-$(pwd)}
here=$(cd $(dirname ${0})/ && pwd)
config_options_file="${here}/benchmark-config-options.txt"
. ${here}/define-config-options.sh

for ((line=1;line<=${nbr_options_lines};line++))
do
  echo ${line}
  ${here}/copy-one-benchmark-output.sh ${line} ${dest}
done
echo "."
