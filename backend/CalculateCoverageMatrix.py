
import sys
import os
import re
from subprocess import Popen,PIPE


parentFolder=sys.argv[1]
cnvBedFile=sys.argv[2]
mappingCutoff=sys.argv[3]
readCutoff=sys.argv[4]
samtool=sys.argv[5]
bedtool=sys.argv[6]
base=sys.argv[7]
cpm_output=base+"_cpmMatrix.csv"
raw_output=base+"_rawMatrix.csv"

mat={}
mat2={}
files_all = [d for d in os.listdir(parentFolder) if os.path.isfile(os.path.join(parentFolder, d))]
files = [d for d in files_all if bool(re.search('bam$', d))]
for file in files:  
    cell=file.replace(".bam", "")
    print(file+"\t")
    mapfile=os.path.join(parentFolder,file)
    p1=Popen([samtool,"view","-q",mappingCutoff,mapfile],stdout=PIPE)
    p2=Popen(["grep","-v","ERCC"],stdin=p1.stdout,stdout=PIPE)
    p1.stdout.close()
    p3=Popen(["wc","-l"],stdin=p2.stdout,stdout=PIPE)
    p2.stdout.close()
    cpm=int(p3.communicate()[0].rstrip())
    p2.stdout.close()
    print("\t"+"count: "+str(cpm)),
    if cpm>int(readCutoff):
        print("PASS")
        mat[cell]={}
        mat2[cell]={}
        p0=Popen([samtool,"view","-b","-q",mappingCutoff,mapfile],stdout=PIPE)
        p1=Popen([bedtool,"coverage","-a",cnvBedFile, "-b", "stdin","-counts"],stdin=p0.stdout,stdout=PIPE)
        p0.stdout.close()
        cnvResult=p1.communicate()[0]
        p1.stdout.close()
        for line in cnvResult.split("\n")[:-1]:
            line=line.rstrip()
            cpmExpr=str(float(line.split("\t")[-1])/cpm*1000000)
            idf=line.split("\t")[3]
            mat[cell][idf]=cpmExpr
            mat2[cell][idf]=line.split("\t")[-1]
    else:
        print("FAIL")
if len(mat.keys())>0:
    cpmMatrix=open(cpm_output,"w")
    rawMatrix=open(raw_output,"w")
    for idf in mat[mat.keys()[0]].keys():
        cpmMatrix.write("\t"+idf)
        rawMatrix.write("\t"+idf)
    cpmMatrix.write("\n")
    rawMatrix.write("\n")
    for cell in mat.keys():
        cpmMatrix.write(cell)
        rawMatrix.write(cell)
        for idf in mat[mat.keys()[0]].keys():        
            cpmMatrix.write("\t"+mat[cell][idf])
            rawMatrix.write("\t"+mat2[cell][idf])
        cpmMatrix.write("\n")
        rawMatrix.write("\n")
    cpmMatrix.close()
    rawMatrix.close()
    

    
            
