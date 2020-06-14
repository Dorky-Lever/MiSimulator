import numpy as np
import heapq as hq


ChromoLen = [
    197070000, 182000000, 160000000, 155030000, 152000000,
    150000000, 145000000, 132000000, 124000000, 130000000,
    122000000, 120000000, 121000000, 124000000, 103000000,
    98000000,   95000000,  91000000,  61000000,
    166000000,  16000000 ]

LinkageLen = [
    127.0, 114.0, 119.2, 84.0, 92.0,
     75.0,  74.0,  82.0, 79.0, 77.0,
     80.0,  66.0,  80.0, 69.0, 81.0,
     72.0,  81.6,  60.0, 55.7,
     96.5
]


class Mouse:
    def __init__(self, arg1, arg2):
        if (type(arg1) is str) and (type(arg2) is str):
            # constructor for Mouse(strain, sex)
            self.generation = 0
            self.chromosome = [[[(ChromoLen[i], arg1)],[(ChromoLen[i], arg1)]] for i in xrange(len(ChromoLen)-2)]
            if (arg2 == "Female"):
                self.chromosome.append([[(ChromoLen[19], arg1)],[(ChromoLen[19], arg1)]])
            else:
                self.chromosome.append([[(ChromoLen[19], arg1)],[(ChromoLen[20], arg1)]])
        elif ((type(arg1) is list) and (type(arg2) is int) and (len(arg1) == 20)):
            # constructor for Mouse(chromoList, generation)
            self.generation = arg2
            self.chromosome = arg1
        else:
            raise TypeError("Invalid arguments for Mouse()")
        self.age = 0
    
    def sex(self):
        # checks the lengths of the sex chromosome haploids
        if (self.chromosome[-1][0][-1][0] == self.chromosome[-1][1][-1][0]):
            return "Female"
        else:
            return "Male"
    
    def intervals(self):
        N = len(self.chromosome)
        total = 0
        for i in xrange(N):
            count = 0
            if (self.chromosome[i][0][-1][0] != self.chromosome[i][1][-1][0]):
                count += len(self.chromosome[i][0])
                count += len(self.chromosome[i][1])
            else:
                j, k = 0, 0
                while (j < len(self.chromosome[i][0])) or (k < len(self.chromosome[i][1])):
                    if (self.chromosome[i][0][j][0] < self.chromosome[i][1][k][0]):
                        j += 1
                    elif (self.chromosome[i][0][j][0] > self.chromosome[i][1][k][0]):
                        k += 1
                    else:
                        j += 1
                        k += 1
                    count += 1
            total += count
        return total        

    #assuming inbred animals for now
    def segmentSize(self):
	N = len(self.chromosome)
	sizes = []
	for i in xrange(N):
	    prevPos = 0
	    for j in xrange(len(self.chromosome[i][0])):
		pos = self.chromosome[i][0][j][0]
		size = pos - prevPos
		sizes.append(size)
		prevPos = pos	    
	
	return sizes
    
    def hetFraction(self):
        N = len(self.chromosome)
        hetro = 0.0
        total = 0.0
        for i in xrange(N):
            start = 0
            homoz = 0
            # test only chromosome pairs of the same length (not male sex XY) 
            if (self.chromosome[i][0][-1][0] == self.chromosome[i][1][-1][0]):
                j, k = 0, 0
                while (j < len(self.chromosome[i][0])) or (k < len(self.chromosome[i][1])):
                    # if tops of lists alleles match we're homozygous
                    match = (self.chromosome[i][0][j][1] == self.chromosome[i][1][k][1])
                    if (self.chromosome[i][0][j][0] < self.chromosome[i][1][k][0]):
                        end = self.chromosome[i][0][j][0]
                        j += 1
                    elif (self.chromosome[i][0][j][0] > self.chromosome[i][1][k][0]):
                        end = self.chromosome[i][1][k][0]
                        k += 1
                    else:
                        end = self.chromosome[i][0][j][0]
                        j += 1
                        k += 1
                    if match:
                        homoz += (end - start)
                    start = end
                hetro += float(ChromoLen[i] - homoz)
                total += float(ChromoLen[i])
        return hetro/total
    
    def founderFraction(self, founders):
	N = len(self.chromosome)
        founderFrac = [0 for i in xrange(len(founders))]
        total = 0.0
        for i in xrange(N):
            start = 0
            # test only chromosome pairs of the same length (not male sex XY) 
            if (self.chromosome[i][0][-1][0] == self.chromosome[i][1][-1][0]):
                j, k = 0, 0
                while (j < len(self.chromosome[i][0])) or (k < len(self.chromosome[i][1])):
                    # if tops of lists alleles match we're homozygous
                    match = (self.chromosome[i][0][j][1] == self.chromosome[i][1][k][1])
		    f1 = self.chromosome[i][0][j][1]
		    f2 = self.chromosome[i][1][k][1]		    
                    if (self.chromosome[i][0][j][0] < self.chromosome[i][1][k][0]):
                        end = self.chromosome[i][0][j][0]
                        j += 1
                    elif (self.chromosome[i][0][j][0] > self.chromosome[i][1][k][0]):
                        end = self.chromosome[i][1][k][0]
                        k += 1
                    else:
                        end = self.chromosome[i][0][j][0]
                        j += 1
                        k += 1
                    if match:
			founderFrac[founders.index(f1)]+=(end-start)*2.0
                    else:
			founderFrac[founders.index(f1)]+=(end-start)
			founderFrac[founders.index(f2)]+=(end-start)
                    start = end
                total += float(ChromoLen[i])*2.0
	for i in xrange(len(founderFrac)):
	    founderFrac[i] = founderFrac[i]/total
	return founderFrac
        
    def myHetFraction(self):
        N = len(self.chromosome)        
        totalGenomeLen = 2644100000 # total of all lengths of chromosomes from ChromoLen
        hetLen = 0
        prevPos = 0
        pos1 = 0
        pos2 = 0
        for i in xrange(N):
            j, k = 0, 0
            if (self.chromosome[i][0][-1][0] == self.chromosome[i][1][-1][0]):            # doesn't account for hets in X-chromo in males
                while (j < len(self.chromosome[i][0])) or (k < len(self.chromosome[i][1])):
                    pos1 = self.chromosome[i][0][j][0]
                    pos2 = self.chromosome[i][1][k][0]
                    f1 = self.chromosome[i][0][j][1]
                    f2 = self.chromosome[i][1][k][1]
                    #print pos1, pos2, f1, f2               
                    if (pos1 < pos2):
                        if(f1 == f2):
                            prevPos = pos1
                        else:
                            hetLen+= pos1 - prevPos
                            prevPos = pos1                        
                        j += 1
                    elif (pos1 > pos2):
                        if(f1 == f2):
                            prevPos = pos2
                        else:
                            hetLen+= pos2 - prevPos
                            prevPos = pos2
                        k += 1
                    else:
                        if(f1 == f2):
                            prevPos = pos1
                        else:
                            hetLen+= pos1 - prevPos
                            prevPos = pos1                              
                        j += 1
                        k += 1
        print hetLen, totalGenomeLen
        return (float(hetLen)/float(totalGenomeLen))                 
  
    def printPicture(self, title):
        picture = []
        picture.append("set,ppc,20")
        picture.append("title,top,%s" % title)
        picture.append("grid,2.0")
        picture.append("key")
        picture.append("set,overlap,gray")
        N = len(self.chromosome)
        for i in xrange(N):
            if (i < N-1):
                line = "chr%d,," % (i+1)
                start = 0
                for end, type in self.chromosome[i][0]:
                    line+= "%s,%d,%d," % (type, start, end)
                    start = end
                picture.append(line)
                picture.append("gap0,%d" % (i+1))
                line = "chr%d,," % (i+1)
                start = 0
                for end, type in self.chromosome[i][1]:
                    line+= "%s,%d,%d," % (type, start, end)
                    start = end
                picture.append(line)
                picture.append("gap")
            else:
                line = "chr20,X,"
                start = 0
                for end, type in self.chromosome[i][0]:
                    line+="%s,%d,%d," % (type, start, end)
                    start = end
                picture.append(line)
                picture.append("gap0")
                if (self.chromosome[i][1][-1][0] == ChromoLen[19]):
                    line = "chr20,X,"
                else:
                    line = "chr21,Y,"
                start = 0
                for end, type in self.chromosome[i][1]:
                    line+= "%s,%d,%d," % (type, start, end)
                    start = end
                picture.append(line)
        return picture

    def printPOV(self):
        N = len(self.chromosome)
        for i in xrange(N):
            for j in xrange(2):
                print "#declare C%02dH%d_List = array[%d]" % (i+1, j, len(self.chromosome[i][j])),
                segList = [type+"," for end, type in self.chromosome[i][j]]
                print "{" + "".join(segList)[0:-1] + "}"
                brkList = [float(end)/ChromoLen[i] for end, type in self.chromosome[i][j][0:-1]]
                if (len(brkList) == 0):
                    print "#declare C%02dH%d_Brks = array[%d]" % (i+1, j, 1),
                    print "{1.0}"
                else:
                    print "#declare C%02dH%d_Brks = array[%d]" % (i+1, j, len(self.chromosome[i][j])-1),
                    format = "{" + (len(brkList)*"%4.2f,")[0:-1] + "}"
                    print format % tuple(brkList)
        print "#declare Mosaic = array[%d][4] {" % (N+1)
        for i in xrange(1,N):
            print "  { C%02dH0_List, C%02dH0_Brks, C%02dH1_List, C%02dH1_Brks}," % (i,i,i,i)
        i = N
        print "  { C%02dH0_List, C%02dH0_Brks, C%02dH0_List, C%02dH0_Brks}," % (i,i,i,i)
        print "  { C%02dH1_List, C%02dH1_Brks, C%02dH1_List, C%02dH1_Brks}}" % (i,i,i,i)

#Make a certain type of mouse, known as a breeder. 
#A breeder's role is to, yep you guessed it, breed. Breeders are kept to ensure that X amount of mice are ready for experimental / timed mates. 
#They can only mate Maxmates times.
#The number of breeder pairs are optimised  

#Promote a mouse to a breeder by breeder_mouse = Breeder(mouse.chromosome, mouse.generation, mouse, Max_number_of_mates)  
class Breeder(Mouse):
    def __init__(self, arg1, arg2, Promotedmouse, Maxmates, Maxage):
        '''Had to  copy the init from Mouse'''
        if (type(arg1) is str) and (type(arg2) is str):
            # constructor for Mouse(strain, sex)
            self.generation = 0
            self.chromosome = [[[(ChromoLen[i], arg1)],[(ChromoLen[i], arg1)]] for i in xrange(len(ChromoLen)-2)]
            if (arg2 == "Female"):
                self.chromosome.append([[(ChromoLen[19], arg1)],[(ChromoLen[19], arg1)]])
            else:
                self.chromosome.append([[(ChromoLen[19], arg1)],[(ChromoLen[20], arg1)]])
        elif ((type(arg1) is list) and (type(arg2) is int) and (len(arg1) == 20)):
            # constructor for Mouse(chromoList, generation)
            self.generation = arg2
            self.chromosome = arg1
        else:
            raise TypeError("Invalid arguments for Mouse()")
        self.maxmates = Maxmates
        self.maxage =  Maxage

    def GetNumOffFromDist(self, desired_Genotype, numOffspring):
         '''Randomly generates the number of offspring of the desired breeders produced within a litter.
            A breeder is produced if the number generated random distribution - U(0, 100) is greater the desired genotype as a percentage '''
         genodist = np.random.randint(100, size=(numOffspring))
         numBreedOff = 0
         for x in xrange(numOffspring):
             if (genodist[x] < (desired_Genotype*100)):
                 numBreedOff += 1
         return numBreedOff

    def MergetoTrios(self, Pedigree):
        '''If trios are used, Odd breeder numbers are merged together'''

        OddPedigree = [breeder for breeder in Pedigree if breeder[1]%2 != 0 and breeder[1]>0] 
        
        for i in xrange(len(OddPedigree)//2):
            two_litters = hq.nsmallest(2, OddPedigree)
            
            merged_litter = (np.amax([two_litters[0][0], two_litters[1][0]]), two_litters[0][1] + two_litters[1][1] )
           
            OddPedigree.remove(two_litters[0])
            
            OddPedigree.remove(two_litters[1])

            Pedigree.remove(two_litters[0])

            Pedigree.remove(two_litters[1])
             
            Pedigree.append(merged_litter)

            
        Pedigree.sort()


        return Pedigree

    
    def MakePedigree(self, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup):
         '''Generates a pedigree where mice are optimally breed until a sufficent number of breeder mice with the rate limiting genotype
            This function works by using  the heap structure (i.e. heapq - hq) to determine optimal breeding structure. 
            Input: 
                       desired_Genotype_Percent = proportion of mice born with the rate limiting genotype (expressed as decimal e.g. 50% -> 0.5)
                       goal_num                 = number of mice with the rate limiting genotype required 
                       numOffspring             = average litter size
                       nBPrs                    = number of starting breeding pairs
                       tBirthtomate             = time to produce an adult mouse (i.e. mating + pregancy + weaning + maturity 
 
             Output: 
                       Pedigree                 = a heap of active breeders that contribute to the goal number. Each breeder is a tuple of
                                                  the time taken to produce the breeder and the number of breeder pairs. 
                       Breeder_List             = a heap of used breeders that were to generate the active breeders within the goal num 
                       numBreedOff              = number of breeders with the desired genotype that were truly produced
                       tmin                     = the time taken to produce the first breeder with the desired genotype
                       tmax                     = the time taken to produce the last breeder with the desired genotype
                       BreederArray             = number of breeders with the desired genotype produced in the G1 litter
                       TimeArray                = times taken to produce mice within the G1 litter
                     '''
         # Initialise Final Variables
         Pedigree = []
         Breeder_List = []
         

         #Peform the G1 breeding from the initial breeding pairs. 
         #Results are stored within two arrays (for simple visualisation) and a heap (for further calculations)
         TimeArray = np.zeros((self.maxmates,nBPrs))
         BreederArray = np.zeros((self.maxmates,nBPrs))
         G1Pedigree = []
         t = tBirthtoMate
         #breed for maxmates times for each BPr
         for litter in xrange(self.maxmates):
            t += tBetween
            tmp_list = []
            for pair in xrange(nBPrs):
                 numBreedOff = self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                 BreederArray[litter,pair] = numBreedOff
                 TimeArray[litter, pair] = t
                 #Add breeder and time to the heap - Pedigree
                 hq.heappush(Pedigree, (TimeArray[litter,pair], BreederArray[litter,pair]) )
                 tmp_list.append((TimeArray[litter,pair], BreederArray[litter,pair]))
            G1Pedigree.append(tmp_list)
         
         #Calculate number of breeders and the min and max for the G1 breeding
         numBreeders = sum([pair[1] for pair in Pedigree])
         tmin = np.amin([pair[0] for pair in Pedigree])
         tmax = np.amax([pair[0] for pair in Pedigree])
         #print "G1 summary: numBreedOff", numBreeders, "tmin", tmin, "tmax", tmax

         #merge cages to trios if you have trios
         if cage_setup == 2:
             Pedigree = self.MergetoTrios(Pedigree)
         #If the G1 breeding produced an insufficient breeder mice (most likely the case) - repeatedly pick all G1 breeder pairs 
         #that took the shortest time to breed and breed to produce G2 --> Gn mice. Check if the goal number is reached each time a litter is produced
         if numBreedOff < goal_num:
             not_goal = True
             

             while not_goal:

                 #merge cages
                 if cage_setup == 2:
                     Pedigree = self.MergetoTrios(Pedigree)

                 #Find all mice with tmin and breed with them
                 tminbreed = np.where([pair[0] for pair in Pedigree] == tmin)

                 #Set up all possible breeder pairs, breed the first litter and use hq.heappushpop to transfer the breeder to the breeder list (so they don't contribute
                 # to goal number) . Newborn is the count of offspring wit hthe rate limiting genotype produced. Pedigree[0][1] is number of BPrs within the Breeder
                 for x in xrange(np.shape(tminbreed)[1]):
                     newBorn = 0
                     for i in xrange(int(Pedigree[0][1])//cage_setup):
                         newBorn += self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                     if Pedigree[0][1] > (cage_setup-1):
                         Breeder_List.append(hq.heappushpop(Pedigree, (tmin+tBirthtoMate+tBetween, newBorn) ))
                     else:
                         Breeder_List.append(hq.heappop(Pedigree))


                 if cage_setup == 2:
                     Breeder_List = self.MergetoTrios(Breeder_List)
                 #Breed the remaining litter pairs with the set up breeders

                 #litter count starts at 1 for calculations to be correct
                 litter_count = 2
                 
                 for litter in xrange(self.maxmates-1):
                     
                     #All breeder pairs (BPrs) should have the same time and as such, should be prioritised first. 
                     for BPr in xrange(np.shape(tminbreed)[1]):
                         #As they are more than a single generation within the G2 litters, the number of BPrs (calculated by numMums) is not constant. As such, 
                         #the number of offspring produced needs to ccalculated inividually for each breeder pair. 
                         numBorn = 0
                         
                         #retreive the number of breeders with the desired genotype - count increments for each BPr
                         
                         #Store the newly generated  litter within the heap - if its Gn Heappushpop, if its just G2 heappush

                         if hq.nlargest((np.shape(tminbreed)[1]-BPr), Breeder_List)[-1][1] > 0: 
                             
                             numMums = hq.nlargest((np.shape(tminbreed)[1]-BPr), Breeder_List)[-1][1]

                             for mum in xrange(int(numMums)):
                                 numBorn += self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                             
                             hq.heappush(Pedigree, (tmin+tBirthtoMate+(tBetween*litter_count), numBorn) )

                         #As heap sort is not a true sort, I run heap_sort twice to ensure that the order is correct - TODO see if  Pedigree.sort() is better 
                         
                         

                         Pedigree.sort()

                         #print Pedigree
                         #calculate, the times and number of true breeders
                         #tmin = np.nanmin([pair[0] for pair in Pedigree])
                         tmax = np.amax([pair[0] for pair in Pedigree])
                         
                         #sums the number of breeders within the active  breeders, ensure that mice that will be too old to breed are not counted 
                         #TODO the if pair[0] > ... might need to be turned off in some situations - I should encode a flag for this.  
                     
                         numBreedOff = sum([pair[1] for pair in Pedigree if pair[0] > (tmax - self.maxage + tBirthtoMate)]) if allatonce else sum([pair[1] for pair in Pedigree])
                         
                         #if the goal number has been reached, optimise and return the results
                         if numBreedOff >= goal_num:

                             not_goal = False
                             #print "G2 --> Gn Make Summary: numBreedOff", numBreedOff, "tmin", tmin, "tmax", tmax, "\n"
                             return Pedigree, Breeder_List, numBreedOff, tmax, G1Pedigree
                     litter_count += 1
                 tmin = np.amin([pair[0] for pair in Pedigree])
         
                     
         return Pedigree, Breeder_List, numBreedOff, tmax, G1Pedigree

    def OptimiseBreeding(self, Breeder_List, Pedigree, G1Pedigree, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup):
        '''Replaces older generation mice with newer generation mice that take less time to produce'''
        
        #Replace G1 mice with G2
        tmax = np.nanmax([pair[0] for pair in Pedigree])
        
        numBreedOff = sum([pair[1] for pair in Pedigree if pair[0] > (tmax - 150 + tBirthtoMate)])

        tcount = tBetween


        for pup in xrange(self.maxmates):

            for breeder in Breeder_List:
                right_ages = [litter for litter in Pedigree if litter[0] == (breeder[0]+tBirthtoMate+tcount) ]
                next_ages = [litter for litter in Pedigree if litter[0] == (breeder[0]+tBirthtoMate+tcount+tBetween) ]

                if len(next_ages) < len(right_ages) and right_ages != [] and (right_ages[0][0]+tBetween) < tmax and right_ages[0][1] > 0:
                    
                    bad_litter = hq.nlargest(1, Pedigree)[0]
                    tmp_Pedigree = [litter for litter in Pedigree if litter not in hq.nlargest(1, Pedigree)]
                    ExtendedBorn = 0
                    for x in xrange(int(breeder[1])):
                        ExtendedBorn += self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                    
                    hq.heappush(tmp_Pedigree, (right_ages[0][0]+tBetween, ExtendedBorn) )
                    
                    tmp_Pedigree.sort()

                    tmp_tmax = np.amax([pair[0] for pair in Pedigree])
                    tmp_numBreedOff = sum([pair[1] for pair in Pedigree if pair[0] > (tmp_tmax - self.maxage + tBirthtoMate)]) if allatonce else sum([pair[1] for pair in Pedigree])

                        #Need to recurvisely
                    if tmp_tmax <= tmax:
                        numBreedOff, Pedigree  = tmp_numBreedOff, tmp_Pedigree

                        if cage_setup == 2:
                             Pedigree = self.MergetoTrios(Pedigree)

                        Pedigree.sort()
                        tmp_G1Pedigree = [litter for row in G1Pedigree for litter in row if litter[0] not in bad_litter]
                        G1Pedigree = [hq.nsmallest(x, tmp_G1Pedigree)[-nBPrs:] for x in range(nBPrs, len(tmp_G1Pedigree) + nBPrs, nBPrs)]
                        tmax = np.amax([pair[0] for pair in Pedigree])
                tcount += tBetween


        #Breed with G2 --> Gn mice (i.e. to reach G3 --> Gn+1) if it improves time
        
        count = 1
        for item in Pedigree:

            
            if hq.nsmallest(count,Pedigree)[-1] not in Breeder_List and hq.nsmallest(count,Pedigree)[-1][1] > 0 and (hq.nsmallest(count,Pedigree)[-1][0]+tBirthtoMate+tBetween) < tmax:

               copies = Pedigree.count(hq.nsmallest(count,Pedigree)[-1])
               
               for x in xrange(copies):
                   numMums = hq.nsmallest(count,Pedigree)[-1][1]
                   tmouse = hq.nsmallest(count,Pedigree)[-1][0]+tBirthtoMate+tBetween
                   numBorn = 0
                   for mum in xrange(int(numMums)):
                        numBorn += self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                   Breeder_List.append(hq.nsmallest(count,Pedigree)[-1])
                   Pedigree.remove(hq.nsmallest(count,Pedigree)[-1])
                   hq.heappush(Pedigree, (tmouse, numBorn) )

                   Pedigree.remove(hq.nlargest(1,Pedigree)[-1])

                   if cage_setup == 2:
                       Pedigree = self.MergetoTrios(Pedigree)
                       Breeder_List = self.MergetoTrios(Breeder_List)

                   tmax = np.amax([pair[0] for pair in Pedigree])

                   

                   #Do the remainder of breeding (i.e. self.maxmates -1) with Breederlist[-1]
                   tmouse += tBetween
                   while tmouse < tmax:
                       numMums = Breeder_List[-1][1]
                       for mum in xrange(int(numMums)):
                           numBorn += self.GetNumOffFromDist(desired_Genotype_Percent, numOffspring)
                       hq.heappush(Pedigree, (tmouse, numBorn) )
                       Pedigree.sort()
                       
                       if cage_setup == 2:
                           Pedigree = self.MergetoTrios(Pedigree)
                       
                       tmouse += tBetween


            else:
                count += 1
                 


        #print "G2 --> Gn Opt Summary: numBreedOff", numBreedOff, "tmin", tmin, "tmax", tmax, "\n"

        return Pedigree, Breeder_List, numBreedOff, tmax, G1Pedigree

    def RemoveExcess(self, Breeder_List, Pedigree, G1Pedigree, desired_Genotype_Percent, goal_num, numOffspring, nBPrs, tBirthtoMate, tBetween, allatonce, cage_setup):
        '''As G2 mice typically have higher breeder counts than G1 mice. Replacing G1 mice with G2 mice results
           in excess breeeders! This code removes excess breeders with the highest time (out of both G1 and G2)'''
        
        try:
            tmax = np.amax([pair[0] for pair in Pedigree])
        except ValueError:
           tmax = np.amax([pair[0] for pair in G1Pedigree])

        numBreedOff = sum([pair[1] for pair in Pedigree if pair[0] > (tmax - self.maxage + tBirthtoMate)]) if allatonce else sum([pair[1] for pair in Pedigree])

        excess = False if numBreedOff < (goal_num+numOffspring) else True
        while excess and Pedigree != []:

            tmaxG1 = np.amax([pair[0] for pair in G1Pedigree])
            tmaxG2 = np.amax([pair[0] for pair in Pedigree])
            

            bad_litter = hq.nlargest(1, Pedigree)[0]

            #remove bad litter and transfer bad breeder back into the pedigree if bad litter was the only  litter produced by the pedigree
            tmp_Pedigree = [litter for litter in Pedigree if litter != bad_litter ]

            if Breeder_List != []:
                if bad_litter[0] == (Breeder_List[-1][0] + tBirthtoMate + tBetween):
                
                    hq.heappush(tmp_Pedigree, Breeder_List[-1])
                    tmp_Breeder_List = Breeder_List[:-1]
                else:
                    tmp_Breeder_List = Breeder_List

            
            #calculate new tmax and numbreedoff
            try:
                tmp_tmax = np.amax([pair[0] for pair in tmp_Pedigree])
                tmp_numBreedOff = sum([pair[1] for pair in tmp_Pedigree if pair[0] > (tmp_tmax - self.maxage + tBirthtoMate)]) if allatonce else sum([pair[1] for pair in tmp_Pedigree])
            except ValueError: 
                flat_G1 = [breeder for row in G1Pedigree for breeder in row]
                tmp_tmax = np.amax([pair[0] for pair in flat_G1])
                tmp_numBreedOff = sum([pair[1] for pair in flat_G1 if pair[0] > (tmp_tmax - self.maxage + tBirthtoMate)]) if allatonce else sum([pair[1] for pair in flat_G1])
            if tmaxG1 >= tmaxG2 and tmp_numBreedOff > goal_num:

                tmp_tmp_G1Pedigree = [litter for row in G1Pedigree for litter in row if litter[0] not in bad_litter]
                tmp_G1Pedigree = [hq.nsmallest(x, tmp_tmp_G1Pedigree)[-nBPrs:] for x in range(nBPrs, len(tmp_tmp_G1Pedigree) + nBPrs, nBPrs)]
                G1Pedigree, Pedigree, tmax, numBreedOff  = tmp_G1Pedigree, tmp_Pedigree, tmp_tmax, tmp_numBreedOff
                Pedigree.sort()

            elif tmaxG1 < tmaxG2 and tmp_numBreedOff > goal_num:
                Pedigree, tmax, numBreedOff, Breeder_List  = tmp_Pedigree, tmp_tmax, tmp_numBreedOff, tmp_Breeder_List
                Pedigree.sort()

            excess = False if tmp_numBreedOff <= goal_num else True

        #print "G2 --> Gn Rem Summary: numBreedOff", numBreedOff, "tmin", tmin, "tmax", tmax, "\n"

        return Pedigree, Breeder_List, numBreedOff, tmax, G1Pedigree



class AdjacencyList(object):
    '''Graph structure used to print pedigrees and which litter came from each used breeder'''
    def __init__(self):
        self.List = {}
    
    def addEdge(self, fromVertex, toVertex):
        # check if vertex is already present
        if fromVertex in self.List.keys():
            self.List[fromVertex].append(toVertex)
        else:
            self.List[fromVertex] = [toVertex]
    
    def printList(self):
        for i  in self.List:
            #only gives litter information about successful breeders. 
            if i[1]:
                print "Breeder Pair:", i, 'Gave birth to',' -> '.join([str(j) for j in self.List[i] ])
            else:
                print "Breeder Pair: ", i, "Gave no litters"

def PrintGraph(ActiveBreeders, UsedBreeders, nBPrs, maxmates, G1Array, tBetween, tBirthtoMate):
    '''Prints G2 Breeding scenarios'''
    #initialises an empty adjacency list for eahc breeder
    BPrGraph = [AdjacencyList() for i in xrange(len(UsedBreeders))]
    #Removes unsucessful litters and G1 litters within the active breeders
    #ActiveBreeders = [breeder for breeder in ActiveBreeders if breeder[1] > 0 ]
    ActiveBreeders = [breeder for breeder in ActiveBreeders if breeder not in G1Array]
    
    
    #for each litter
    printed_Breeders = []

    for index, breeder in enumerate(UsedBreeders):
        tcount = tBetween
        #Tries to add the litter to the corresponding breeder. Create a new branch if the branch doesn't exists 
        for pup in xrange(maxmates):
            right_ages = [litter for litter in ActiveBreeders if litter[0] == (breeder[0]+tBirthtoMate+tcount) ]
            ngen_breeders = [litter for litter in UsedBreeders if litter[0] == (breeder[0]+tBirthtoMate+tcount) ]
            first_litter = right_ages[0] if right_ages !=  [] else []
            first_breeder = ngen_breeders[0] if ngen_breeders !=  [] else []
            
             

            if first_breeder != [] and first_litter == []:
                Breeder_to_tree = [mum for mum in ngen_breeders if mum not in printed_Breeders]
                
                if Breeder_to_tree != [] and breeder[1] > 0:
                    BPrGraph[index].addEdge(breeder, Breeder_to_tree[0])
                    printed_Breeders.append(Breeder_to_tree[0])

            elif first_litter != [] and breeder[1] > 0:
                try:
                   BPrGraph[index].addEdge(breeder, first_litter)
                   ActiveBreeders.remove(first_litter)
                   
                except IndexError:
                   BPrGraph[index] = BPrGraph[index]
            
            
            tcount += tBetween
    #Prints the graph
    for graph in BPrGraph:
        if graph != None:
            graph.printList()
       











