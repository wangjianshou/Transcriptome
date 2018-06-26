#usage python3 qcSum.py qcdir sampleName1 sampleName2...
import os
from os import path
import sys
import re
import glob
import collections

#sample RawReads CleanReads CleanReadsRate NReadsRate AdatperReadsRate RawQ30BaseRate CleanQ30BaseRate rRNARate Exon Intron Intergenic
header = ["sample", "RawReads", "CleanReads", "CleanReadsRate", "NReadsRate",
          "AdatperReadsRate", "RawQ30BaseRate", "CleanQ30BaseRate", "rRNARate",
          "MappingRate", "Exon", "Intron", "Intergenic"
         ]
qcdir = os.path.abspath(sys.argv[1])
qc = open(path.join(qcdir, "qcSum.xls"), "w")
qc.write('\t'.join(header)+'\n')
stat = collections.OrderedDict()
for i in header:
  stat[i] = ''
for sample in sys.argv[2:]:
  sampleDir = os.path.join(qcdir, sample)
  stat_clean = glob.glob(sampleDir+"/CleanData/*clean.rRNA.fq.gz.stat")
  assert len(stat_clean)==1, "statfile should only one"
  with open(stat_clean[0], "r") as f:
    stat_1 = dict([line.strip().split('\t') for line in f])
  stat['sample'] = sample
  stat['RawReads'] = "%.3f" % (int(stat_1["Original reads number"])/10**6)
  stat["CleanReads"] = "%.3f" % (int(stat_1["Clean reads number"])/10**6)
  stat["CleanReadsRate"] = "%.3f" % float(stat_1["Clean reads rate(%)"])
  stat["NReadsRate"] = "%.3f" % float(stat_1["Ns reads rate(%)"])
  stat["AdatperReadsRate"] = "%.3f" % float(stat_1["Adapter polluted reads rate(%)"])
  stat["RawQ30BaseRate"] = "%.3f" % float(stat_1["Original Q30 bases rate(%)"])
  stat["CleanQ30BaseRate"] = "%.3f" % float(stat_1["Clean Q30 bases rate(%)"])
  rRNA_file = glob.glob(sampleDir+"/CleanData/*.rRNA.stat")
  assert len(rRNA_file)==1, "statfile should only one"
  with open(rRNA_file[0], "r") as f:
    stat["rRNARate"] = re.sub(r"% overall alignment rate", '', f.readlines()[-1].strip())
  bamfile = glob.glob(sampleDir+"/BAM/*.bam.stat")
  assert len(rRNA_file)==1, "statfile should only one"
  with open(bamfile[0], "r") as f:
    stat["MappingRate"] = re.sub(r"% overall alignment rate", '', f.readlines()[-1].strip())
  mappingDistribute = glob.glob(sampleDir+"/ReadCount/*.reads.distributes")
  assert len(mappingDistribute)==1, "statfile should only one"
  with open(mappingDistribute[0], "r") as f:
    stat["Exon"], stat["Intron"], stat["Intergenic"] = ['%.3f' % (float(i)*100) for i in f.readlines()[-1].strip().split()]
  qc.write('\t'.join([stat[i] for i in stat])+'\n')
  





