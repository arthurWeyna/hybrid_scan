##This script outputs the number of snp per contig, using vcf as input

import scipy.stats
import re, sys

vcf=sys.argv[1]


conts={}
for l in open(vcf, "r"):
    if re.search("^##contig", l):
        conts[re.sub(",length=.*$", "", re.sub("##contig=<ID=", "", l)).strip()]=[re.sub("^.*,length=", "", re.sub(">$", "", l)).strip(), 0]
    elif not re.search("^#", l):
        l2=l.strip().split("\t")
        conts[l2[0]][1]+=1



for k in conts.keys():
    print(k+" "+conts[k][0]+" "+str(conts[k][1]))
