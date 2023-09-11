#! /bin/bash

eta="$1"
en="$2"
nevents="$3"
idx="$4"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_13_2_0_pre3/samples/200_PU"

cmsRun python/Config.py $eta $en $nevents $idx $root
