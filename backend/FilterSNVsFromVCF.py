import sys

#This script takes 2 VCF files, eg generated with GATK, containing all SNVs detected in exome-seq of a tumor and a control sample. Additionally, it takes a BED file which contains regions with somatic CNVs as an input. 

normalVcf=open(sys.argv[1],"r")
tumorVcf=open(sys.argv[2],"r")
region_bed=open(sys.argv[3],"r")
base=sys.argv[4]
res_file=open(base+"_germline_snvs.bed","w")

regions={}
for line in region_bed:
	line=line.split("\t")
	if not regions.has_key(line[0]):
		regions[line[0]]=[[line[1],line[2],line[0]+":"+line[1]+"-"+line[2]]]
	else:
		regions[line[0]].append([line[1],line[2],line[0]+":"+line[1]+"-"+line[2]])

cands={}
for line in normalVcf:
	if not line[0]=="#":
		f=line.rstrip().split("\t")
		dp=int(f[-1].split(":")[2])
		ref=float(f[-1].split(":")[1].split(",")[0])
		alt=float(f[-1].split(":")[1].split(",")[1])
		gt=f[-1].split(":")[0]
		qscore=f[6]
		vaf=alt/(ref+alt)
		rs=f[2]
		chr=f[0]
		pos=f[1]
		if gt=="0/1" and dp>7 and qscore=="PASS" and len(f[3])==1 and len(f[4])==1:
			#Focus on heterozygous germline SNVs
			if vaf >0.2 and vaf <0.8:
				if regions.has_key(chr):
					writeit=False
					idf=""
					for reg in regions[chr]:
						if  int(pos) > int(reg[0]) and int(pos) < int(reg [1]): 
							writeit=True
							idf=reg[2]
					if writeit==True:
						cands[f[0]+"\t"+f[1]]=[vaf,idf]
				
#print(cands)
#print(median(cands.values()))			

for line in tumorVcf:
	if not line[0]=="#":
		f=line.rstrip().split("\t")
		dp=int(f[-1].split(":")[2])
		ref=float(f[-1].split(":")[1].split(",")[0])
		alt=float(f[-1].split(":")[1].split(",")[1])
		gt=f[-1].split(":")[0]
		qscore=f[6]
		vaf=alt/(ref+alt)
		if cands.has_key(f[0]+"\t"+f[1]):
			if dp>7 and qscore=="PASS" and len(f[3])==1 and len(f[4])==1:
				#Thrshold for reliable identification of minor (B) allele
				if vaf >0.6 or vaf< 0.4:
					#At least 10 percent difference in VAF
					if abs(cands[f[0]+"\t"+f[1]][0]-vaf)>0.1:
						a=f[3]
						b=f[4]
						if vaf>0.6:
							a=f[4]
							b=f[3]
							vaf=1-vaf
							cands[f[0]+"\t"+f[1]][0]=1-cands[f[0]+"\t"+f[1]][0]
						res_file.write(f[0]+"\t"+f[1]+"\t"+f[1]+"\t"+cands[f[0]+"\t"+f[1]][1]+"_"+a+"_"+b+"_"+str(cands[f[0]+"\t"+f[1]][0])+"_"+str(vaf)+"\t"+str(int(ref+alt))+"\t+\n")

