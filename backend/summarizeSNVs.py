import sys 
from subprocess import Popen,PIPE
import os

#With a Folder containing al BED files, a file name for the result, a metadata table and a BED file with genomic coordinates of SNVs, and the genome reads were mapped to in FASTA format this script calculates the BAF in every cell for each locus. 

parentFolder=sys.argv[1]
res_name=sys.argv[2]
metadata=sys.argv[3]
bedFile=sys.argv[4]
genome=sys.argv[5]

outf=open(res_name+"_baf.txt","w")
outf1=open(res_name+"_af.txt","w")
outf2=open(res_name+"_bf.txt","w")

nucs={}
valid_cells=[]

for line in open(metadata,"r"):
	line=line.split("\t")
	if not line[0]=="ID":
		if line[-3]=="PASS":
			valid_cells.append(line[0])
			
for line in open(bedFile,"r"):
	line=line.split("\t")
	id=line[0]+":"+line[1]
	ref=line[3].split("_")[1]
	alt=line[3].split("_")[2]
	nucs[id]=[ref,alt,line[3]]

i=0
dirs = [d for d in os.listdir(parentFolder) if os.path.isdir(os.path.join(parentFolder, d))]
for folder in dirs:  
	cellName=folder.split("-")[-1]
	if cellName in valid_cells:
		bamName="accepted_hits.bam"
		bamPath=os.path.join(parentFolder,folder, bamName)
		molemonkey={}
		molemonkey1={}
		molemonkey2={}
		if os.path.exists(bamPath): 
			#if cellName in valid_cells:
			p1=Popen(["bam-readcount","-q","20","-b","20","-w","0","-l",bedFile,"-f",genome,bamPath],stdout=PIPE)
			res=p1.communicate()[0]
			for line in res.split("\n")[:-1]:
				line=line.split("\t")
				idf=line[0]+":"+line[1]
				hasi={}
				hasi["A"]=float(line[5].split(":")[1])
				hasi["C"]=float(line[6].split(":")[1])
				hasi["G"]=float(line[7].split(":")[1])
				hasi["T"]=float(line[8].split(":")[1])
				refallele=nucs[idf][0]
				varallele=nucs[idf][1]
				total=hasi[refallele]+hasi[varallele]
				ba=hasi[varallele]
				aa=hasi[refallele]
				if total>0:
					baf=ba/total
					molemonkey[idf]=str(baf)
					molemonkey1[idf]=str(aa)
					molemonkey2[idf]=str(ba)
				else:
					molemonkey[idf]="NA"
					molemonkey1[idf]="NA"
					molemonkey2[idf]="NA"
			if i==0:
				for key in nucs.keys():outf.write("\t"+key+"_"+nucs[key][2]);outf2.write("\t"+key+"_"+nucs[key][2]);outf1.write("\t"+key+"_"+nucs[key][2])
				outf.write("\n")
				outf1.write("\n")
				outf2.write("\n")
				i+=1
			outf.write(cellName)
			outf1.write(cellName)
			outf2.write(cellName)
			for key in nucs.keys():
				if molemonkey.has_key(key):
					outf.write("\t"+molemonkey[key])
					outf1.write("\t"+molemonkey1[key])
					outf2.write("\t"+molemonkey2[key])
				else:
					outf.write("\tNA")
					outf1.write("\tNA")
					outf2.write("\tNA")
			outf.write("\n")
			outf1.write("\n")
			outf2.write("\n")
					
outf1.close()
outf2.close()
outf.close()

	
	
