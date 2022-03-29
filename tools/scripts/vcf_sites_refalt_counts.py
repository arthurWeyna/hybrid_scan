##This scripts preprocesses vcf files to output ref and alt allelic depths in a simpler format. This is so snp_refalt_plot.R works faster.

import re, sys

vcf=sys.argv[1]


for l in open(vcf, "r"):
    if not re.search("^#", l):
        l2=l.strip().split("\t")
        gt=l2[9].split(":")
        if not re.search("\.", gt[0]): 
            dps=sorted([int(n) for n in gt[2].split(",")], reverse=True)
            print(l2[0]+" "+l2[1]+" "+str(dps[0])+" "+str(dps[1]))



