#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   test-all-config-options.sh : Test all defined configurations.
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

here=$(cd $(dirname ${0})/ && pwd)
. ${here}/define-config-options.sh

args=$*

for ((line=1;line<=${nbr_options_lines};line++))
do
  echo ${line}
  ${here}/fast-test-one-config-option.sh ${line} ${args}
done
