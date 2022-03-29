##This script applies a filter for contamination on a vcf. The principle is to remove sites that are more likely to be true homozygous sites than true heterozygous sites. This is done by using the binomial distribution and assuming that the expected frequency of the alt allele is either 0.5 (true heterozygous sites) or user-chosen frequency prob_alt_null < 0.5 (true homozygous sites). prob_alt_null is non-null because of sequencing error and contamination. Higher values will result in harder filtering.


import scipy.stats
import re, sys

vcf=sys.argv[1]


prob_alt_null=float(sys.argv[2])


conts={}
for l in open(vcf, "r"):
    if re.search("^#", l):
        print(l.rstrip())
    else:
        l2=l.strip().split("\t")
        gt=l2[9].split(":")
        if not re.search("\.", gt[0]): 
            dps=sorted([int(n) for n in gt[2].split(",")], reverse=True)
            probhet=scipy.stats.binom.pmf(dps[0], dps[0]+dps[1], 0.5)
            probhomo=scipy.stats.binom.pmf(dps[1], dps[0]+dps[1], prob_alt_null)
            if probhet > probhomo:
                print(l.rstrip())



