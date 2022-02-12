#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   check-config-options.sh : Check common definitions for systematic tests.
#
#   begin                : Wed 2022-01-19
#   copyright            : (C) 2022 by Paul C. Leopardi
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   For licence details, see ${package_dir}/COPYING

# set -x

here=$(cd $(dirname ${BASH_SOURCE[0]})/ && pwd)
package_dir=$(cd ${here}/.. && pwd)

config_options_file=${config_options_file:-"${package_dir}/test/config-options.txt"}

vars=$(grep '$[_A-Z]*/' $config_options_file | sed 's/^.*\$\([_A-Z]*\).*$/\1/' | sort -u)

for varname in $vars; do
  var=${!varname}
  if [[ "$var" == "" ]]; then
    echo "$varname is undefined or empty."
    exit 1
  elif [[ ! -d $var ]]; then
    echo "$var : no such directory."
    exit 1
  fi
done
