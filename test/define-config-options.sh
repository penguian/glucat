#!/usr/bin/env bash

#   GluCat : Generic library of universal Clifford algebra templates
#   define-config-options.sh : Common definitions for systematic tests.
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

here=$(cd $(dirname ${BASH_SOURCE[0]})/ && pwd)
package_dir=$(cd ${here}/.. && pwd)

config_options_file=${config_options_file:-"${package_dir}/test/config-options.txt"}
config_options_pattern='^\([a-z][a-z-]*\):\(.*\)$'
nbr_options_lines=$(wc -l ${config_options_file} | awk '{print $1}')

match_option_pattern(){
  line=$1
  which=$2
  echo "$(sed -n ${line}'s/'${config_options_pattern}'/\'${which}'/p' ${config_options_file})"
}
