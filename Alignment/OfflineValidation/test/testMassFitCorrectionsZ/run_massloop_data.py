# Runs all the steps
# Once for Iter0: ./massscales_data ./massfit ./resolfit
# For niter times: ./massscales_data (command updated to use results of the previous iteration) ./massfit ./resolfit 
# Can be run in TOYS MODE to loop over many toys, where in ./massscales_data pseudodata with known AeM bias is generated from MC
# Authors: Cristina Alexe, Lorenzo Bianchini

import argparse
import os
import sys
import copy
import math
import ROOT
import time
   
parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--out_folder',   default='out' , help = 'name of output subdirectory')
parser.add_argument('--tag',   default='PostVFP' , help = 'type of data used')
parser.add_argument('--pathToDataFiles',   default='./inoutfiles/Run2022_E_F_G-reshaped/' , help = 'path to data')
parser.add_argument('--lumiData',  type=float, default=-1. , help = 'for massscales') # 27.4970
parser.add_argument('--pathToMCFiles',   default='./inoutfiles/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_CVH_reshaped/' , help = 'path to MC')
parser.add_argument('--lumiMC',  type=float, default=-1. , help = 'for massscales')# 14.241
parser.add_argument('--nominalResolutionFile',   default='./inoutfiles/NominalResolution/Run3Summer22EE_nominal_resolution_coefficients.root' , help = 'nominal resolution')
parser.add_argument('--niter', dest = 'niter'  , type = int,  default=1, help='number of iterations after the 0th')
parser.add_argument('--forceIter', dest = 'forceIter'  , type = int,  default=-1, help='will only do a specific iteration and skip the rest')
parser.add_argument('--ntoys', dest = 'ntoys' , type = int, default=0, help='number of toys, default is 0, set to >0 for TOYS MODE')

args = parser.parse_args()

def loop_one(seed, toy_number):    

    assert args.forceIter <= args.niter 

    out_folder=args.out_folder
    tag = args.tag
    pathToDataFiles = args.pathToDataFiles
    pathToMCFiles = args.pathToMCFiles
    lumiData = args.lumiData
    lumiMC = args.lumiMC
    nominalResolutionFile = args.nominalResolutionFile

    cmd_histo_iter0 = './massscales_data --firstIter=-1 --lastIter=2 '+\
        ' --out_folder='+out_folder+' '+\
        ' --tag='+tag+' '+\
        ' --run=Iter0 '+\
        ' --pathToDataFiles='+pathToDataFiles+' '+\
        ' --pathToMCFiles='+pathToMCFiles+' '+\
        ' --nRMSforGausFit=-1 '+\
        ' --minNumEvents=100 --minNumEventsPerBin=10 '+\
        ' --minNumMassBins=4 '+\
        ' --rebin=2 '+\
        ' --fitNorm --fitWidth '+\
        f" --lumiData={lumiData} --lumiMC={lumiMC} --scaleToData "

    if args.ntoys>0 : # in TOYS MODE overwrite the tag and the Iter 0 massscales command
        tag = args.tag+'_toy'+str(toy_number)
        cmd_histo_iter0 = './massscales_data --firstIter=-1 --lastIter=2 '+\
        ' --out_folder='+out_folder+' '+\
        ' --tag='+tag+' '+\
        ' --run=Iter0 '+\
        ' --pathToDataFiles='+pathToDataFiles+' '+\
        ' --pathToMCFiles='+pathToMCFiles+' '+\
        ' --nRMSforGausFit=-1 '+\
        ' --minNumEvents=100 --minNumEventsPerBin=10 '+\
        ' --minNumMassBins=4 '+\
        ' --rebin=2 '+\
        ' --fitNorm --fitWidth '+\
        ' --toysMode '+\
        ' --nominalResolutionFile='+nominalResolutionFile+' '+\
        ' --biasResolutionRange=0.1 '+\
        ' --scaleToData ' +\
        ' --seed='+str(seed)
        
    if not args.forceIter>0:
        print(cmd_histo_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_histo_iter0)
        print("\n")
    cmd_fit_iter0 = './massfit --ntoys=1 --bias=-1 '+\
        ' --out_folder='+out_folder+' '+\
        '--tag='+tag+' '+\
        '--run=Iter0 '
    if not args.forceIter>0:
        print(cmd_fit_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_fit_iter0)
        print("\n")
    cmd_resol_iter0 = './resolfit --ntoys=1 --bias=-1 '+\
        ' --out_folder='+out_folder+' '+\
        ' --tag='+tag+' '+\
        ' --run=Iter0 '+\
        ' --nominalResolutionFile='+nominalResolutionFile+' '+\
        ' --maxSigmaErr=0.1 '
    if not args.forceIter>0:
        print(cmd_resol_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_resol_iter0)
        print("\n")

    for iter in range(1, args.niter+1):
        if (args.forceIter>0 and iter!=args.forceIter) or args.forceIter==0 :
            continue
        cmd_histo_iteri = cmd_histo_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        cmd_histo_iteri += ' --usePrevMassFit '+\
            ' --tagPrevMassFit='+tag+' '+\
            ' --runPrevMassFit=Iter'+str(iter-1)+' '
        cmd_histo_iteri += ' --usePrevResolFit '+\
            ' --tagPrevResolFit='+tag+' '+\
            ' --runPrevResolFit=Iter'+str(iter-1)+' '
        print(cmd_histo_iteri)
        if not args.dryrun:
            os.system(cmd_histo_iteri)
            print("\n")
        cmd_fit_iteri = cmd_fit_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_fit_iteri)
        if not args.dryrun:
            os.system(cmd_fit_iteri)
            print("\n")
        cmd_resol_iteri = cmd_resol_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_resol_iteri)
        if not args.dryrun:
            os.system(cmd_resol_iteri)
            print("\n")
    return


if __name__ == '__main__':
    start = time.time()
    if args.ntoys==0 :
        print('Running on data and MC')
        loop_one(seed=0,toy_number=0)
        end = time.time()
        print('Done', args.niter, 'iterations in', (end - start)/60., 'min.')
    else :
        iseed = 4357
        # Run over many toys
        for itoy in range(0, args.ntoys):
            # The seed must increase each time by the total_number_of_threads*10 + 2 + 1
            iseed += itoy*3843
            print('Running toy with seed '+str(iseed))
            loop_one(seed=iseed,toy_number=itoy)
            '''
            if args.ntoys>1 :
                for iter in range(0, args.niter+1):
                    if args.forceIter>0 and iter!=args.forceIter:
                        continue
                    cmd_hadd = 'hadd -f massfit_'+args.tag+'_merged_Iter'+str(iter)+'.root massfit_'+args.tag+'_toy*_Iter'+str(iter)+'.root'
                    print(cmd_hadd)
                    if not args.dryrun:
                        os.system(cmd_hadd)
            '''
        end = time.time()
        print(args.ntoys, 'toys run in', (end - start)/60., 'min.')
