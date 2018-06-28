#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import vcf

if len(sys.argv) == 2:
	vcf_file = sys.argv[1]
else:
	sys.exit("Enter a vcf file")


vcf_reader = vcf.Reader(open(vcf_file, "r"))
# vcf_writer = vcf.Writer(open(out, "w"), vcf_reader)


out_tab = []


header = ["CHROM","POS","QD","REF","ALT"]
header.extend(vcf_reader.samples)


for record in vcf_reader:
	if len(record.FILTER) > 0:
		continue
	try:
		qd = str(record.INFO["QD"])
	except:
		qd = "-" 
	out_tab.append([record.CHROM,str(record.POS),qd,record.REF,str(record.ALT)])
	alts = list(record.ALT)

	for sample in record.samples:
		try:
			### get GQ and DP to print
			gq=str(sample["GQ"])
			dp=str(sample["DP"])

			if sample.is_variant:
				afs=[]
				### get AF for each allele (REF and all the ALT)
				for af in  sample["AD"]:
					if sample["DP"] > 0:
						afs.append(str(round(float(af/float(sample["DP"])),2)))
					else:
						afs.append(str(0))
				max_af = np.argmax(afs)

				if max_af == 0:
					alleles = str(record.REF)
				else:
					alleles = str(alts[max_af-1])

				out_tab[-1].extend([gq+"/"+dp+"/"+alleles+":"+afs[max_af]])
			else:
				out_tab[-1].extend([gq+"/"+dp+"/."])

		except: ### there's no information available
			out_tab[-1].extend(["././."])

print "\t".join(header)
for line in out_tab:
	print  "\t".join(line)
