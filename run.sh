#! /bin/bash

eta="$1"
en="$2"
nevents="$3"
idx="$4"
input_dir="$5"

cmsRun python/Config.py $eta $en $nevents $idx $step $input_dir
