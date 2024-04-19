from Analyzers.Ntuplizer.Ntuplizer_cfi import *

def customiseTICLForKalmanFilterNtuplizer(process):
    process.kfNtuplizer = ntuplizer.clone()
    return process
