

####################--------------------IMPORT 
import numpy as np
import os, sys, argparse, random, copy, time

####################--------------------FUNCTIONS 
##########----------Read bed file containing list of TEs to simulate
##########	each element is a list containing info of TE
##########		[Contig, start, end, name, strand]
##########		e.g. ['scaffold00001', 4324, 4388, 'DNA/hAT-Ac/1', '-']
def getTE(bed):
	bedEntries = []
	with open(bed, 'r') as file:
		for line in file:
			arr=line.split()
			bedEntries.append([arr[0],int(arr[1]),int(arr[2]),arr[3],arr[5]])
	return bedEntries


#####-----Check if any other REF TE overlaps with focal (chrm, start, end)
#####	Returns 1 if Ref TE (with borders increased by L) overlaps at all with focal
#####	L is read length
#####	x[1], x[2] is start and stop of bed files
def isNested(chrm, start, stop, ID, bedEntries, L):
    for x in bedEntries:
        if ID != x[3] and chrm==x[0] and x[1]-L < start and stop < x[2]+L:
            return 1
        if ID != x[3] and chrm==x[0] and x[1]-L < start and start < x[2]+L:
            return 1
        if ID != x[3] and chrm==x[0] and x[1]-L < stop and stop < x[2]+L:
            return 1
    return 0


#####-----Check if focal (chrm, start, stop) contains any other ref TE nested within it
#####	Only returns 1 if TE completely contained in focal (with borders increased by L)
#####		if just partial overlap, still returns 0
#####	L is read length
#####	x[1], x[2] is start and stop of bed files
def containsNested(chrm, start, stop, ID, bedEntries, L):
    for x in bedEntries:
        if ID != x[3] and chrm==x[0] and start-L < x[1] and x[2] < stop+L:
            return 1
    return 0


##########----------Get number of eligible deletions
##########	Eligible if:
##########		len between min and max
##########		not overlap with other TEs
##########		no TEs nested within		
def numEligibleDeletions(mnLen, mxLen, rLen, bedEntries):
	#####-----loop through each TE
	count=0
	for x in bedEntries:
		LEN = x[2] - x[1]
		nestedTE = containsNested(x[0],x[1],x[2],x[3],bedEntries,rLen)
		overlapTE = isNested(x[0],x[1],x[2],x[3],bedEntries,rLen)
		if mxLen > LEN >= mnLen and nestedTE == 0 and overlapTE == 0:
			count+=1
	return count


##########----------Return dict of sequences in chromosome specified by TE bed file
##########	key=TE name, value=seq extracted from genome
##########	each TE must have unique name
##########		if TE share name, then only the last one will be saved
def extractSeq(chrom, bedEntries):
	extractedSeqs={}
	for x in bedEntries:
		extractedSeqs[x[3]]=""
		for i in xrange(x[1]-1,x[2]):
			extractedSeqs[x[3]]+=chrom[i]
	return extractedSeqs


##########----------Read refChrom (fasta) and extract sequences in bedEntries
##########	if refChrom contain N's, then return error
def getBedSequences(refChrom, bedEntries):
	validBase="ACGT"
	invalidIndex,ct=[],0
	rawChr,acgtN = "",[0,0,0,0,0]
	with open(refChrom, "r") as fIN:
		for line in fIN:
			#####-----get Cm name
			if line.startswith(">"):
				cmName=line.split()[0][1:]
			else:
				#####-----replace newline with blank
				seq=line.replace("\n","").upper()
				##### count num of A,C,G,T,N; record index if N
				for s in seq:
					if s in validBase:
						acgtN[validBase.index(s)]+=1
						ct+=1
					else:
						acgtN[-1]+=1
						invalidIndex.append(ct)
						ct+=1
				#####-----append onto rawChr
				rawChr += seq
	#####-----return error if there are N's
	if acgtN[-1] > 0:
		print('ERROR: Invalid refChrom, there are {} N at these indices {}'.format(acgtN[5], invalidIndex))
	else:
		#####-----extract TE seq from rawChr, using bedEntries
		teFasta=extractSeq(rawChr, bedEntries)
		return cmName, teFasta


##########----------Read a fasta file and return a list where each element is a nucleotide
def getCmAsList(fasta):
	#####-----open fasta
	raw = ''
	with open(fasta, "r") as fIN:
		for line in fIN:
			#####-----append seq (capitalized) into raw
			if not line.startswith(">"):
				seq=line.replace("\n","").upper()
				raw += seq
	#####-----convert string to list
	rawList=list(raw)
	return rawList


##########----------Generate random TSD length from possion dist with mean lambda
def getTSD(lam):
    tsd = np.random.poisson(lam, 1)[0]
    return tsd


##########----------Given a string seq, returns as fasta where each line has 70 characters
def fastaformat(seq):
    return '\n'.join(seq[i:i+70] for i in xrange(0,len(seq),70))


##########----------Given the_list, returns elements (in order) that are != val
def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]


##########---------Randomly select TEs to delete WITHOUT replacement (of course cannot delete twice)
##########	delPool is a list of [CmName, start, end, TE_name, strand]
def getRandomDeletions(nDel, mnLen, mxLen, rLen, bedEntries):
	delPool, iCheck = [], []
	#####-----choose nDel number of deletions
	while len(delPool) < nDel:
		#####random choose TE
		i = random.randint(0,len(bedEntries)-1)
		#####if i not chosen before
		if not i in iCheck:
			####if i lies on focal Cm and between min and max length
			LEN = bedEntries[i][2] - bedEntries[i][1]
			nestedTE = containsNested(bedEntries[i][0], bedEntries[i][1], bedEntries[i][2], bedEntries[i][3], bedEntries, rLen)
			overlapTE = isNested(bedEntries[i][0], bedEntries[i][1], bedEntries[i][2], bedEntries[i][3], bedEntries, rLen)
			if mnLen <= LEN <= mxLen and nestedTE == 0 and overlapTE == 0:
				#####add i to delPool
				delPool.append(bedEntries[i])
				#####append i to iCheck, so will not be chosen again
				iCheck.append(i)
	#####-----deep copy so further alterations will not affect bedEntries
	delPool = copy.deepcopy(delPool)
	#####-----Sort deletions by start position([1])
	delPool.sort(key=lambda x: x[1])
	#####-----return
	return delPool


##########---------Randomly select TEs to insert WITHOUT replacement
##########	inPool is a list of [CmName, start, end, TE_name, strand]
def getRandomInsertions(nIns, mnLen, mxLen, rLen, cmName, cmLen, bedEntries):
	####-----Get nIns number of insertions
	inPool, iCheck = [],[]
	while len(inPool) < nIns:
		#####random select TE
		i=random.randint(0,len(bedEntries)-1)
		#####if not already chosen
		if not i in iCheck:
			#####if i between min and max length and does not contain anyother TE nested within
			LEN = bedEntries[i][2] - bedEntries[i][1]
			nestedTE = containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen)
			if mnLen <= LEN <= mxLen and nestedTE == 0:
				#####add i to inPool
				inPool.append(bedEntries[i])
				#####append i to iCheck, so will not be chosen again
				iCheck.append(i)
	#####-----deep copy so further alterations will not affect bedEntries
	#####	i.e. append data below will not affect bedEntries
	inPool = copy.deepcopy(inPool)
	#####-----Get random insertion site for each TE
	insPos = []
	while len(insPos) < len(inPool):
		#####random pos from rLen to cmLen-rLen
		x = random.randint(rLen,(cmLen-rLen))
		#####if x not nested within existing TE and not chosen before
		if isNested(cmName, x, x, "NOVEL_TE", bedEntries, rLen) == 0 and not x in insPos:
			#####keep insertion position
			insPos.append(x)
	#####-----Add insertion position and length of TSD to each element in inPool
	for i in xrange(len(inPool)):
		inPool[i].append(insPos[i])
		inPool[i].append(getTSD(5))
	#####-----Sort insertions by insertion position([5])
	inPool.sort(key=lambda x: x[5])
	#####-----return
	return inPool
	

##########----------Simulate TE indels into given sequence
##########	seq: list where each element is one Nt of the chromosome
##########	teFasta: dictionary where key=TE name, value=TE sequence string
##########	delPool: [1] and [2] are start and end of TE to be deleted
##########	inPool:	[5] is insertion position of new TE
##########			[3] is the name of TE (to be referenced in teFasta)
##########			[6] is length of TSD
##########	Steps to simulate INDEL:
##########		0) seq = ['A','A','C','T','G','A','C','G','T','T,'T,'A','A']
##########		1) seq = ['A','A','C','%','%','%','C','G','T','T,'T,'A','A']		#replace deleted sites with %
##########		2) seq = ['A','A','C','%','%','%','C','G','$','T','T,'T,'A','A']	#insertion $ at insertion position 
##########		3) seq = ['A','A','C','C','G','$','T','T,'T,'A','A']				#remove all % (DELETION)
##########		4) seq = ['A','A','C','C','G','ACTGACTGACTG','T','T,'T,'A','A']		#replace $ with teSeq (INSERTION)
##########		5) joinedSeq = 'AACCGACTGACTGACTGTTTAA'								#join all elements to create string sequence
def simIndel(chrom, teFasta, inPool, delPool):
	#####-----make deep copy of chrom
	seq = copy.deepcopy(chrom)
	#####-----1) loop through delPool and replace Nt to be deleted with '%'
	for x in delPool:
		#####replace letter in seq with '%' if it is between of x[0]-1 and x[1]
		for i in xrange(x[1]-1,x[2]):
			seq[i]='%'
	#####-----2) loop through inPool and add '$' to site where new TE will be inserted
	count=0
	for i in xrange(len(inPool)):
		#####add '$' into seq at insertion site of i+count
		#####	since TEs ordered by insetion position, each additional $ need to increase the subsequent insertion position by 1
		seq.insert(inPool[i][5] + count,'$')
		count += 1
	#####-----3) DELETE TE
	#####	remove all '%' elements
	seq = remove_values_from_list(seq,'%')
	#####-----4) INSERT TE
	#####	replace '$' with TE sequence
	id=0
	for i in xrange(len(seq)):
		#####replace '$' with TE
		if seq[i] == '$':
			#####find teSeq from teFasta
			teSeq = ''
			for te in teFasta:
				if inPool[id][3] in te:
					teSeq = copy.deepcopy(teFasta[te])
					break
			#####create TSD as sequence from i-[6] to i
			tsd = ''.join(seq[i-inPool[id][6]:i])
			#####teSeq is teSeq from teFasta and tsd
			teSeq = teSeq+tsd
			#####repalce seq[i] (which is $) with teSeq
			seq[i] = teSeq
			#print(i, id, len(teSeq), len(tsd), inPool[id][3], inPool[id][6])
			#####increment id, so move to next insertion in inPool
			id += 1
	#####-----5) Join all elements in seq to create one string
	#####	all '%' were removed, all '$' were TE seqs
	joinedSeq=''.join(seq)
	#####-----return
	return joinedSeq


##########---------Randomly select TEs to insert WITH replacement
##########	insPool is a list of [CmName, start, end, TeName, strand, insPos, lenTSD, ancHomolog, maHomolog, typeMut]
##########		num of mut equal nMut + nHet
##########	nHet is the number of shared heterozygosities in ANC and MA
##########		inserts TE onto same homolog on AND and MA
##########		typeMut = 5
##########	typeMut controls how an MA mut inserts into ANC and MA
##########		1 = novel insertion; nothing in ANC, het in MA
##########		2 = novel deletion; homo in ANC, het in MA
##########		3 = LOH insertion; het in ANC, homo in MA
##########		4 = LOH deletion; het in ANC, nothing in MA
##########		5 = shared Het; het in ANC, het in MA (present on same homolog)
def getRandomMutations(nMut, typeMut, nHet, mnLen, mxLen, rLen, cmName, cmLen, bedEntries):
	#####-----Create list of mutation types and random order
	listType = [typeMut]*nMut + [5]*nHet
	random.shuffle(listType)
	#####-----Get nMut+nHet number of insertions
	insPool = []
	while len(insPool) < len(listType):
		#####random select TE
		i = random.randint(0,len(bedEntries)-1)
		#####if i between min and max length and does not contain anyother TE nested within
		LEN = bedEntries[i][2] - bedEntries[i][1]
		nestedTE = containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen)
		if mnLen <= LEN <= mxLen and nestedTE == 0:
			#####add i to insPool
			insPool.append(copy.deepcopy(bedEntries[i]))
	#####-----deep copy so further alterations will not affect bedEntries
	#####	i.e. append data below will not affect bedEntries
	insPool = copy.deepcopy(insPool)
	#####-----Get random insertion site for each TE
	insPos = []
	while len(insPos) < len(insPool):
		#####random pos from rLen to cmLen-rLen
		x = random.randint(rLen, (cmLen-rLen))
		#####if x not nested within existing TE and not chosen before
		if isNested(cmName, x, x, "NOVEL_TE", bedEntries, rLen) == 0 and not x in insPos:
			#####keep insertion position
			insPos.append(x)
	#####-----Add insertion position, length of TSD and which homolog to insert for AND and MA
	for i in xrange(len(insPool)):
		#####ins position
		insPool[i].append(insPos[i])
		#####length of TSD
		insPool[i].append(getTSD(5))
		#####type of mutation
		#####	all have a het site, so requires random insertion into homolog 1 or 2 
		randHomolog = random.randint(1,2)
		if listType[i] == 1:
			insPool[i].append(0)
			insPool[i].append(randHomolog)
		elif listType[i] == 2:
			insPool[i].append(3)
			insPool[i].append(randHomolog)
		elif listType[i] == 3:
			insPool[i].append(randHomolog)
			insPool[i].append(3)
		elif listType[i] == 4:
			insPool[i].append(randHomolog)
			insPool[i].append(0)
		elif listType[i] == 5:
			insPool[i].append(randHomolog)
			insPool[i].append(randHomolog)
		#####append typeMut
		insPool[i].append(listType[i])
	#####-----Sort insertions by insertion position([5])
	insPool.sort(key=lambda x: x[5])
	#####-----return
	return insPool


##########----------Similar to getRandomMutations, BUT all het TE will be inserted into Homolog1
def getRandomMutations_HetH1(nMut, typeMut, nHet, mnLen, mxLen, rLen, cmName, cmLen, bedEntries):
	#####-----Create list of mutation types and random order
	listType = [typeMut]*nMut + [5]*nHet
	random.shuffle(listType)
	#####-----Get nMut+nHet number of insertions
	insPool = []
	while len(insPool) < len(listType):
		#####random select TE
		i = random.randint(0,len(bedEntries)-1)
		#####if i between min and max length and does not contain anyother TE nested within
		LEN = bedEntries[i][2] - bedEntries[i][1]
		nestedTE = containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen)
		if mnLen <= LEN <= mxLen and nestedTE == 0:
			#####add i to insPool
			insPool.append(copy.deepcopy(bedEntries[i]))
	#####-----deep copy so further alterations will not affect bedEntries
	#####	i.e. append data below will not affect bedEntries
	insPool = copy.deepcopy(insPool)
	#####-----Get random insertion site for each TE
	insPos = []
	while len(insPos) < len(insPool):
		#####random pos from rLen to cmLen-rLen
		x = random.randint(rLen, (cmLen-rLen))
		#####if x not nested within existing TE and not chosen before
		if isNested(cmName, x, x, "NOVEL_TE", bedEntries, rLen) == 0 and not x in insPos:
			#####keep insertion position
			insPos.append(x)
	#####-----Add insertion position, length of TSD and which homolog to insert for AND and MA
	for i in xrange(len(insPool)):
		#####ins position
		insPool[i].append(insPos[i])
		#####length of TSD
		insPool[i].append(getTSD(5))
		#####type of mutation
		#####	all het sites will place insertionin Homolog1
		if listType[i] == 1:
			insPool[i].append(0)
			insPool[i].append(1)
		elif listType[i] == 2:
			insPool[i].append(3)
			insPool[i].append(1)
		elif listType[i] == 3:
			insPool[i].append(1)
			insPool[i].append(3)
		elif listType[i] == 4:
			insPool[i].append(1)
			insPool[i].append(0)
		elif listType[i] == 5:
			insPool[i].append(1)
			insPool[i].append(1)
		#####append typeMut
		insPool[i].append(listType[i])
	#####-----Sort insertions by insertion position([5])
	insPool.sort(key=lambda x: x[5])
	#####-----return
	return insPool
	
	
##########---------Check that expected and obserted insertion lengths match
def checkInsertions(ancH1_Exp, ancH2_Exp, maH1_Exp, maH2_Exp, ancH1_Obs, ancH2_Obs, maH1_Obs, maH2_Obs):
	if ancH1_Exp==ancH1_Obs:
		print('Inserted {} bases on ANC.H1'.format(ancH1_Obs))
	else:
		print('ERROR: expected insertion length ({}) not equal observed insertion length ({})'.format(ancH1_Exp, ancH1_Obs))
		sys.exit(0)
	if ancH2_Exp==ancH2_Obs:
		print('Inserted {} bases on ANC.H2'.format(ancH2_Obs))
	else:
		print('ERROR: expected insertion length ({}) not equal observed insertion length ({})'.format(ancH2_Exp, ancH2_Obs))
		sys.exit(0)
	if maH1_Exp==maH1_Obs:
		print('Inserted {} bases on MA.H1'.format(maH1_Obs))
	else:
		print('ERROR: expected insertion length ({}) not equal observed insertion length ({})'.format(maH1_Exp, maH1_Obs))
		sys.exit(0)
	if maH2_Exp==maH2_Obs:
		print('Inserted {} bases on MA.H2'.format(maH2_Obs))
	else:
		print('ERROR: expected insertion length ({}) not equal observed insertion length ({})'.format(maH2_Exp, maH2_Obs))
		sys.exit(0)



####################--------------------MAIN
def main():
	####################--------------------IMPORT ARGUMENTS
	parser = argparse.ArgumentParser()
	parser.add_argument('-wd', dest='wd', help='Full path to working directory (include prefix for new folders)', required=True)
	parser.add_argument('-pirs', dest='pirsPATH', help='Full path to pIRS (including /pirs at the end)', required=True)
	parser.add_argument('-c', dest='refChrom', help='Chromosome to simulate in fasta format (must contain only one sequence)', required=True)
	parser.add_argument('-b', dest='bed', help='TE annotation bed file', required=True)
	parser.add_argument('-nmut', dest='nMut', help='Number of TE mutations to simulate in MA line', type=int, required=True)
	parser.add_argument('-tmut', dest='typeMut', help='Type of mutation to simulate [1, 2, 3, 4, 5]', type=int, required=True)
	parser.add_argument('-nhet', dest='nHet', help='Number of shared heterozygous TE sites to simulate', type=int, required=True)
	parser.add_argument('-ncl', dest='nCl', help='Number of non-focal MA lines to simulate (clones of Ancestor)', type=int, required=True)
	parser.add_argument('-mnlen', dest='mnLen', help='Minimum length of TEs to insert and delete', type=int, default=400)
	parser.add_argument('-mxlen', dest='mxLen', help='Maximum length of TEs to insert and delete', type=int, default=10000)
	parser.add_argument('-r', dest='randseed', help='Seed for random number generator', type=int)
	parser.add_argument('-snp', dest='snpRate', help='Rate of heterozygosity for haploid genome', type=float, default=0.0005)
	parser.add_argument('-x', dest='cov', help='Coverage for pIRS to simulate', type=int, default=50)
	parser.add_argument('-rlen', dest='rLen', help='Read length for pIRS to simulate', type=int, default=150)
	parser.add_argument('-insz', dest='inSize', help='Insert size for pIRS to simulate', type=int,  default=350)
	args = parser.parse_args()
	
	##########----------Assign variables
	cwd = os.path.realpath(args.wd)
	pirsPATH = args.pirsPATH
	refChrom = args.refChrom
	bed = args.bed
	nMut = args.nMut
	typeMut = args.typeMut
	nHet = args.nHet
	nCl = args.nCl
	mnLen = args.mnLen
	mxLen = args.mxLen
	snpRate = args.snpRate
	cov = args.cov
	rLen = args.rLen
	inSize = args.inSize
	##########---------Set random seed
	if args.randseed == None:
		random.seed()
		np.random.seed()
	else:
		random.seed(args.randseed)
		np.random.seed(args.randSeed)

	####################--------------------PREPARE SIMULATION
	##########----------Create new directories for output
	os.system('mkdir '+cwd+'pirsAnc')
	os.system('mkdir '+cwd+'simTE')
	os.system('mkdir '+cwd+'pirsReads')
	out1 = cwd+'pirsAnc/'
	out2 = cwd+'simTE/'
	out3 = cwd+'pirsReads/'
	
	##########----------Output parameters as newline separated file
	listFlag = ['c','b','nmut','tmut','nhet','ncl','mnlen','mxlen','r','snp','x','rlen','insz']
	listValue = [refChrom, bed, nMut, typeMut, nHet, nCl, mnLen, mxLen, snpRate, cov, rLen, inSize]
	with open(cwd+'param.txt', 'w') as fOUT:
		fOUT.write('FLAG\tVALUE\n')
		for FLAG, VALUE in zip(listFlag, listValue):
			fOUT.write(FLAG+'\t'+str(VALUE)+'\n')
	
	####################--------------------PIRS: GENERATE A DIPLOID ANCESTRAL SEQUENCE CONTAINING SNPS (NO INDELS)
	####################	Generates two refChrom contianig SNPs to represent a diploid individual (with heterozygosity)
	####################		output as ANC.H1 and ANC.H2 to represents homolog 1 and 2
	print ('Simulating SNPs for ANC...\n')
	for i in range(2):
		cmd=pirsPATH +' diploid -q'+ ' -s ' +str(snpRate)+ ' -d 0.0 -v 0.0 -o ' + out1 + 'ANC.H' +str(i+1)+ ' ' + refChrom
		os.system(cmd)
	
	####################--------------------SIMULATE TE INDELS
	##########----------Read bed file of TEs
	bedEntries = getTE(bed)
	
	##########----------Extract TE sequences from refChrom given bedEntries
	print ('Extracting TE sequences from reference...\n')
	cmName, teFasta = getBedSequences(refChrom, bedEntries)
	
	##########---------Read chromosome into list where each element in a nucleotide
	#####-----ANC same as refChrom except for SNPs added by PIRS
	ancH1 = getCmAsList(out1+'ANC.H1.snp.fa')
	ancH2 = getCmAsList(out1+'ANC.H2.snp.fa')
	#####-----MA same as ANC initially
	maH1 = getCmAsList(out1+'ANC.H1.snp.fa')
	maH2 = getCmAsList(out1+'ANC.H2.snp.fa')
	#####-----cmLen is same for all chromosomes initially
	cmLen = len(ancH1)
	
	##########----------Get list of TE that will be inserted into ANC and MA
	##########	mutPool is sorted by insPos in ASCENDING order:
	##########		[CmName, start, end, TeName, strand, insPos, lenTSD, ancHomolog, maHomolog, typeMut]
	##########	ancHomolog and maHomology indicates which homolog to insert onto ANC and MA, respectively
	##########		0 = neither, 1 = homolog1, 2 = homolog2, 3 = both homologs
	##########	typeMut indicates mutation mut
	##########		1 = novel insertion; nothing in ANC, het in MA
	##########		2 = novel deletion; homo in ANC, het in MA
	##########		3 = LOH insertion; het in ANC, homo in MA
	##########		4 = LOH deletion; het in ANC, nothing in MA
	##########		5 = shared Het; het in ANC, het in MA (present on same homolog)
	print ('Selecting TEs to insert...\n')
	#mutPool = getRandomMutations(nMut, typeMut, nHet, mnLen, mxLen, rLen, cmName, cmLen, bedEntries)
	mutPool = getRandomMutations_HetH1(nMut, typeMut, nHet, mnLen, mxLen, rLen, cmName, cmLen, bedEntries)
	
	##########----------Output list of TE insertions
	##########	output = [Cm, InsPos, TeName, TeLength, TsdLen, AncHomolog, MaHomolog, TypeMut]
	with open(out2+'simTeList.txt', 'w') as fOUT:
		fOUT.write('Cm\tTeName\tInsPos\tTeLen\tTsdLen\tAncHomolog\tMaHomolog\tTypeMut\n')
		print('Cm\tTeName\tInsPos\tTeLen\tTsdLen\tAncHomolog\tMaHomolog\tTypeMut')
		for x in mutPool:
			fOUT.write(x[0]+'\t'+str(x[3])+'\t'+str(x[5])+'\t'+str(x[2]-x[1]+1)+'\t'+str(x[6])+'\t'+str(x[7])+'\t'+str(x[8])+'\t'+str(x[9])+'\n')
			print(x[0]+'\t'+str(x[3])+'\t'+str(x[5])+'\t'+str(x[2]-x[1]+1)+'\t'+str(x[6])+'\t'+str(x[7])+'\t'+str(x[8])+'\t'+str(x[9]))
	
	##########----------subset mutPool to TEs that will be added to each homolog of ANC and MA
	mutANC_H1 = [MUT for MUT in mutPool if (MUT[7]==1 or MUT[7]==3)]
	mutANC_H2 = [MUT for MUT in mutPool if (MUT[7]==2 or MUT[7]==3)]
	mutMA_H1 = [MUT for MUT in mutPool if (MUT[8]==1 or MUT[8]==3)]
	mutMA_H2 = [MUT for MUT in mutPool if (MUT[8]==2 or MUT[8]==3)]
	
	##########----------Simulate TE indels for each homolog of ANC and MA
	print ('Simulating TE indels for ANC and MA...\n')
	ancH1_indel = simIndel(ancH1, teFasta, mutANC_H1, [])
	ancH2_indel = simIndel(ancH2, teFasta, mutANC_H2, [])
	maH1_indel = simIndel(maH1, teFasta, mutMA_H1, [])
	maH2_indel = simIndel(maH2, teFasta, mutMA_H2, [])
	
	##########----------Check that expected and observed insertion lengths match
	#####-----calc expected insertion lengths
	ancH1_Exp = sum([x[2]-x[1]+1+x[6] for x in mutANC_H1])
	ancH2_Exp = sum([x[2]-x[1]+1+x[6] for x in mutANC_H2])
	maH1_Exp = sum([x[2]-x[1]+1+x[6] for x in mutMA_H1])
	maH2_Exp = sum([x[2]-x[1]+1+x[6] for x in mutMA_H2])
	#####-----Calc observed insertion lengths
	ancH1_Obs = len(ancH1_indel)-cmLen
	ancH2_Obs = len(ancH2_indel)-cmLen
	maH1_Obs = len(maH1_indel)-cmLen 
	maH2_Obs = len(maH2_indel)-cmLen
	#####-----Check expected and observed match
	checkInsertions(ancH1_Exp, ancH2_Exp, maH1_Exp, maH2_Exp, ancH1_Obs, ancH2_Obs, maH1_Obs, maH2_Obs)
	
	##########----------Output new chromosomes into out2
	with open(out2+'simInsLen.txt', 'w') as fOUT:
		fOUT.write('Line\tHomolog\texpIns\tobsIns\n')
		fOUT.write('ANC\tH1\t'+str(ancH1_Exp)+'\t'+str(ancH1_Obs)+'\n')
		fOUT.write('ANC\tH2\t'+str(ancH2_Exp)+'\t'+str(ancH2_Obs)+'\n')
		fOUT.write('MA\tH1\t'+str(maH1_Exp)+'\t'+str(maH1_Obs)+'\n')
		fOUT.write('MA\tH2\t'+str(maH2_Exp)+'\t'+str(maH2_Obs)+'\n')
	
	with open(out2+'ANC.H1.simTE.fa', 'w') as fOUT:
			fOUT.write('>ANC.H1.simTE\n'+fastaformat(ancH1_indel)+'\n')
	
	with open(out2+'ANC.H2.simTE.fa', 'w') as fOUT:
			fOUT.write('>ANC.H2.simTE\n'+fastaformat(ancH2_indel)+'\n')
	
	with open(out2+'MA.H1.simTE.fa', 'w') as fOUT:
			fOUT.write('>MA.H1.simTE\n'+fastaformat(maH1_indel)+'\n')
	
	with open(out2+'MA.H2.simTE.fa', 'w') as fOUT:
			fOUT.write('>MA.H2.simTE\n'+fastaformat(maH2_indel)+'\n')
	
	####################--------------------PIRS: SIMULATE READS FROM DIPLOID ANC AND MA
	print ('Simulating reads for ANC and {} CLONES...'.format(nCl))
	HOMOLOG1 = out2 + 'ANC.H1.simTE.fa'
	HOMOLOG2 = out2 + 'ANC.H2.simTE.fa'
	listName = ['ANC']+['CL'+str(x) for x in range(1,nCl+1)]
	for NAME in listName:
		print('\tSimulating reads for {}\n'.format(NAME))
		cmd=(pirsPATH + ' simulate --diploid --no-substitution-errors --no-indel-errors --no-gc-content-bias -t 1 -q -z' +
			' -l '+str(rLen)+' -x '+str(cov)+' -m '+str(inSize)+
			' -s ' +NAME+ ' -o ' +out3+ ' '+HOMOLOG1+ ' '+HOMOLOG2)
		os.system(cmd)
	
	print ('Simulating reads for one focal MA line...\n')
	HOMOLOG1 = out2 + 'MA.H1.simTE.fa'
	HOMOLOG2 = out2 + 'MA.H2.simTE.fa'
	NAME = 'MA'
	cmd=(pirsPATH + ' simulate --diploid --no-substitution-errors --no-indel-errors --no-gc-content-bias -t 1 -q -z' +
		' -l '+str(rLen)+' -x '+str(cov)+' -m '+str(inSize)+
		' -s ' +NAME+ ' -o ' +out3+ ' '+HOMOLOG1+ ' '+HOMOLOG2)
	os.system(cmd)
	
	print('DONE!')


##########----------Run main()
if __name__ == "__main__":
    main()
	
