#! /bin/bash

eta="$1"
en="$2"
nevents="$3"
idx="$4"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_13_2_0_pre3/KFv0.1/debugging/0_PU"

cmsRun dumper.py $eta $en $nevents $idx $root
