import numpy as np
def parseMSAFile(fastapath):
	"""
	Parses a FASTA multiple sequence alignment file into a dictionary.
	---------------------------------------------------------------------------
	Parameters:
	- fastapath: full path of the FASTA MSA file.
	Returns:
	- A dictionary, where:
		- The keys of the dicitionary depend on the description of the sequences.
		Assuming the format of the description follows the
		<SpeciesTreeID>_<LocusTreeID>_<GeneTreeID> convention. The key will then
		correspond to the gene family of the sequence <SpeciesTreeID>_<LocusTreeID>.
		- Each element under a gene family, will be another dictionary which key is
		the <GeneTreeID> and the elements are the full description and the sequence
		itself.
	"""
	fastafile=open(fastapath, 'r')
	lines=fastafile.readlines()
	fastafile.close()
	seqdata=[]
	seqDict=dict()
	description=""; seq=""; tmp="";count=1
	for line in lines:
		if line.strip().startswith(">"):
			seqdata+=[line.strip()]
		else:
			if(seqdata[-1].startswith(">")):
				seqdata+=[line.strip()]
			else:
				seqdata[-1]+=line.strip()

	for line in seqdata:
		if not (line.strip()==''):
			if (count%2==0):
				seq=line
				try:
					test=seqDict[tag]
				except:
					seqDict[tag]={}

				try:
					seqDict[tag].update({tmp[2]:{\
						'description':description,\
						'sequence':seq\
					}})
				except:
					seqDict[tag][tmp[2]]={}
					seqDict[tag].update({tmp[2]:{\
						'description':description,\
						'sequence':seq\
					}})
				seq=None
				description=None
				tmp=None
			else:
				description=line[1:len(line)]
				tmp=description.split("_")
				tag="{0}_{1}".format(tmp[0], tmp[1])
		count+=1
	return seqDict


def parseMSAFileWithDescriptions(fastapath):
	"""
	Parses a FASTA multiple sequence alignment file into a dictionary.
	---------------------------------------------------------------------------
	Parameters:
	- fastapath: full path of the FASTA MSA file.
	Returns:
	- A dictionary, where:
		- The keys of the dicitionary are the description of the sequences.
		- Each element will contain the correspondin sequences.
	"""
	fastafile=open(fastapath, 'r')
	lines=fastafile.readlines()
	fastafile.close()
	seqDict=dict()
	description=""; seq=""; tmp="";count=1
	seqdata=[]
	for line in lines:
		if line.strip().startswith(">"):
			seqdata+=[line.strip()]
		else:
			if(seqdata[-1].startswith(">")):
				seqdata+=[line.strip()]
			else:
				seqdata[-1]+=line.strip()
	for line in seqdata:
		if not (line.strip()==''):
			if (count%2==0):
				seq=line[0:-1].strip()
				seqDict[tmp]=seq
				seq=None
				description=None
				tmp=None
			else:
				description=line.strip()
				tmp=description[1:len(description)]
		count+=1
	return seqDict

def isFasta(filepath):
	"""
	Checks whether the format of the filepath introduced corresponds to FASTA
	---------------------------------------------------------------------------
	Parameters:
	- fastapath: full path of the FASTA MSA file.
	Returns:
	- A bollean. TRUE if it is a proper fasta, FALSE otherwise.
	"""
	fastaOk=True
	f=open(filepath)
	lines=f.readlines()
	f.close()
	seqdata=[]
	seqDict=dict()
	description=""; seq=""; tmp="";count=1
	for line in lines:
		if line.strip().startswith(">"):
			seqdata+=[line.strip()]
		else:
			if(seqdata[-1].startswith(">")):
				seqdata+=[line.strip()]
			else:
				seqdata[-1]+=line.strip()

	numLinesFile=len(seqdata)
	if (numLinesFile % 2 == 0): # It looks like a fasta, has even num of lines
		numSeqs=len(seqdata)/2
		index=0
		while (index < numSeqs and seqdata[index*2][0]==">") and not (seqdata[index*2][0]=="+"):
			index=index+1

		if not index >= numSeqs: # there must have been an error
			return False
		if seqdata[index*2][0]=="+":
			return False
	else:
		return False
	return True
