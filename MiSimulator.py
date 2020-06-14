#import PySimpleGUI as sg

from Mouse import *
from tqdm import tqdm
import random

import sys

import bisect
import os
import numpy as np

def runSim2wayRegBreedingScheme(N, MaxGenerations):
    """ This method simulates N 2-way breedings and outputs 3 arrays
    Each output array is a histogram of how many times a value occurs
    For example: histo2way[5] = The number of times it took exactly 5 generations for the animals to become fully inbred
    histo2way - # of generations to fully inbred
    histo2way99 - # of generations to reach 99% inbred (the combined heterozygous fraction of the male and female is >= 99%)
    inter2way - # of intervals in the inbred animals, with a max of 199 intervals
    """
    histo2way = [0 for i in xrange(MaxGenerations)]
    histo2way99 = [0 for i in xrange(MaxGenerations)]
    inter2way = [0 for i in xrange(200)]
    for rep in xrange(N):
        newF = Mouse('A', 'Female')
        newM = Mouse('B', 'Male')
        reached99 = False
        for g in xrange(1,MaxGenerations):
            F, M = newF, newM
            newF = mate(F, M, "Female")
            newM = mate(F, M, "Male")
            if(reached99 == False and hetFractionBetweenAnimals(newF,newM) <= 0.01):
                histo2way99[g]+=1
                reached99 = True
            if inbred(newF, newM):
                histo2way[g] += 1
                inter2way[min(newF.intervals(), 199)] += 1
                break
    return (histo2way,inter2way, histo2way99)



def runMisimulator(N, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, maxmates, allatonce, cage_setup):
    ''''This simulates the desired breeding scenario N times and produces summary statistics, the worst case breeding scenario and the best case breeding scenario'''
 
    Tmaxhist = [0 for i in xrange(N)]
    NumBreedhist = [0 for i in xrange(N)]
    PedigreeList = [[] for i in xrange(N)]
    BreederList = [[] for i in xrange(N)] 
    G1PedigreeList = [[] for i in xrange(N)]
    
    with tqdm(total=N, file=sys.stdout) as pbar:
        for rep in xrange(N):
            A = Mouse("A", "Female")
            B = Breeder(A.chromosome, A.generation, A, maxmates, maxage)

            init_Ped, init_Breed, init_Numbreed, init_tmax, init_G1 = B.MakePedigree(desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup)
            pbar.update(0.25)
        #print "Initial Pedigree:\n"

        #print init_Ped
        
        #print "G1 Litter and Times: "

        #for litter in init_G1:
        #    print litter
       
        #print "\nG2 --> Gn Breeders Set up: \n", PrintGraph(init_Ped, init_Breed, 3, 6, init_G1, tBetween, tBirthtoMate) 

            opt_Ped, opt_Breed, opt_Numbreed, opt_tmax, opt_G1 = B.OptimiseBreeding(init_Breed, init_Ped, init_G1, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup)
            pbar.update(0.5)
        
        #print "Opt Pedigree:\n"

        #print opt_Ped

        #print "G1 Litter and Times: "

        #for litter in opt_G1:
        #    print litter

        #print "\nG2 --> Gn Breeders Set up: \n", PrintGraph(opt_Ped, opt_Breed, 3, 6, opt_G1, tBetween, tBirthtoMate)

            PedigreeList[rep], BreederList[rep], NumBreedhist[rep], Tmaxhist[rep], G1PedigreeList[rep] = B.RemoveExcess(opt_Breed, opt_Ped, opt_G1, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup)
            pbar.update(0.25)
        #print "Final Pedigree: \n"

        #print PedigreeList[rep]

        #print "G1 Litter and Times: "

        #for litter in G1PedigreeList[rep]:
        #    print litter

        #print "\nG2 --> Gn Breeders Set up: \n", PrintGraph(PedigreeList[rep], BreederList[rep], 3, 6, G1PedigreeList[rep], tBetween, tBirthtoMate)
        Tmaxhist = [tmax + 49 for tmax in Tmaxhist]
        

    return Tmaxhist, NumBreedhist, PedigreeList, BreederList, G1PedigreeList


""" Methods for collecting statistics """

def findAverageMiceTotal(hist):
    return float(np.mean(np.array(hist)))

def findStdMiceTotal(hist):
    return float(np.std(np.array(hist)))



if __name__ == "__main__":

    print "Welcome to MiSimulator!!! \n"
    
    N = 10000


    Genotype = float(raw_input('What is the rate limiting Genotype? '))

    goal_num = int(raw_input('How many breeders of the rate limiting genotype are required? '))

    cage_setup = int(raw_input('How are cages set up? Enter 1 for Breeder Pairs  -  2 for Trios: '))
 
    pairs =  int(raw_input('How many inital breeder pairs (BPrs) do you have? ')) if cage_setup == 1 else int(raw_input('How many Trios do you have? '))

    #parameters affecting days per generation

    tBirthtoMate = int(raw_input('How long does it take (in days) from the start of breeding until the new pups reach maturity? '))
    
    tBetween = int(raw_input('How long does it take for breeding between litters? '))

    #Parameters affecting mice breeding
    numOffspring = int(raw_input('What is the average litter size (to the nearest unit)? ') )

    maxmates = int(raw_input('How many times can a breeder mate before being disbanded? '))

    maxage = int(raw_input('What age are breeders disbanded? '))

    #Yes / No parameters which affect calculations

    allatonceprompt = raw_input('Are all experimental mice needed to be bred simultaneously (Y/n) ')
    
    allatonce = False if allatonceprompt == 'n' else True

    print "Number of Simulations: ", N

    Tmaxs, Numbreeds, Pedigrees, UsedBreeders, G1Ped = runMisimulator(N, Genotype, goal_num, numOffspring, pairs, tBirthtoMate, tBetween, maxmates, allatonce, cage_setup)


    Tmax_ave = findAverageMiceTotal(Tmaxs)
    Tmax_std = findStdMiceTotal(Tmaxs)

    Num_ave = findAverageMiceTotal(Numbreeds)
    Num_std = findStdMiceTotal(Numbreeds) 


    print "\nSummary of", N, "Simulations: \n"

    print "Time required to generate", goal_num, "breeders with Genotype Ratio of", Genotype, "\n"
    print "Average - ", Tmax_ave, "Standard Deviation - ", Tmax_std, "Shortest Time", np.amin(Tmaxs), "Longest Time", np.amax(Tmaxs), "\n" 


    print "Number of Breeders Reached at Goal:", "Average -", Num_ave, "Standard Deviation -", Num_std, "\n"
    
    real_breeders = [[] for i in xrange(N)]
    for index, breeders in enumerate(UsedBreeders): 
        real_breeders[index] = np.sum([pair[1] for pair in breeders if pair[1] > 0])

    

    #print "Number of G2 --> Gn Breeders Set up: ", "Average -",  np.mean(real_breeders), "Standard Deviation", float(np.std(real_breeders) )


    best_scen = np.where(Tmaxs == np.amin(Tmaxs))[0][0]

    med_scen = np.round(np.where(Tmaxs == np.median(Tmaxs))[0][0])


    print "\nBest Case Scenario - t =", Tmaxs[best_scen], ":", "\n", "G1 Litter and Times: "
    
    for litter in G1Ped[best_scen]:
        print litter

    
    #print "\nNumber of G2 --> Gn Breeders Set up:" , real_breeders[best_scen]
    print "\nG2 --> Gn Breeders Set up: \n", PrintGraph(Pedigrees[best_scen], UsedBreeders[best_scen], pairs, maxmates, G1Ped[best_scen], tBetween, tBirthtoMate)

    print "\nMedian Scenario - t =", Tmaxs[med_scen], ":", "\n", "G1 Litter and Times: "

    for litter in G1Ped[med_scen]:
        print litter
    

    #print "\nNumber of G2 --> Gn Breeders Set up:" , real_breeders[med_scen]
    print "\nG2 --> Gn Breeders Set up: \n", PrintGraph(Pedigrees[med_scen], UsedBreeders[med_scen], pairs, maxmates, G1Ped[med_scen], tBetween, tBirthtoMate)
    
