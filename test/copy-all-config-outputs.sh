#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   copy-all-config-outputs.sh : Copy test output for all defined configurations.
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

here=$(cd $(dirname ${0})/ && pwd)
. ${here}/define-config-options.sh

for ((line=1;line<=${nbr_options_lines};line++))
do
  echo ${line}
  ${here}/copy-one-config-output.sh ${line}
done
