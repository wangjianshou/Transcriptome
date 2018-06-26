#! /usr/bin/python3
'''
__version__=v0.0
__author__=wangjianshou
__email__=jianshouwang@scisoon.cn
'''

import argparse
import configparser
import os
from os import path
import glob
from functools import wraps
import itertools
import collections
from sjm import *

def Decorator(cls):
  class Wrapper:
    def __init__(self, outdir='.', *, config, iterInfo, **kwargs):  # iterInfo should be iterable object
      self.wrapped = cls(**kwargs) if kwargs else ''
      self.iterInfo = iterInfo
      self.__wrapped__ = cls
      self.idx = 0
      self.outdir = outdir
      self.config = config
    def __getattr__(self, attr):
      return getattr(self.wrapped, attr)
    def __iter__(self):
      self.iterobj = itertools.dropwhile(lambda x: x.startswith('#'), self.iterInfo)
      return self
    def __next__(self):
      infoLine = next(self.iterobj)
      while infoLine.strip()=='':
        infoLine = next(self.iterobj)
      sample = infoLine.split()
      qc = self.__wrapped__(fqdir=sample[2], outdir=self.outdir, name=sample[0])
      qc.getCleanData(self.config)
      qc.Mapping(self.config)
      qc.ReadCount(self.config)
      return (sample[1], sample[0])
  return Wrapper

#@Decorator
class QC:
  '''
  QC class
  '''
  def __init__(self, *, fqdir, outdir='.', name=None, R1=None, R2=None, R1adapter=None, R2adapter=None):
    self.R1 = R1 is None and glob.glob(path.join(fqdir, "*R1.fq.gz"))[0] or R1
    self.R2 = R2 is None and glob.glob(path.join(fqdir, "*R2.fq.gz"))[0] or R2
    self.R1adapter = R1adapter is None and glob.glob(path.join(fqdir, "*R1.adapter.txt.gz"))[0] or R1adapter
    self.R1adapter = R2adapter is None and glob.glob(path.join(fqdir, "*R2.adapter.txt.gz"))[0] or R2adapter
    self.name = name is None and path.basename(self.R1).rstrip("_R1.fq.gz") or name
    self.outdir = path.abspath(path.join(outdir, self.name))
    path.isdir(self.outdir) or os.mkdir(self.outdir)
  def getCleanData(self, config):
    self.CleanDir = path.join(self.outdir, "CleanData")
    path.isdir(self.CleanDir) or os.mkdir(self.CleanDir)
    R1_before_rRNA = path.join(self.CleanDir, path.basename(self.R1).rstrip('.fq.gz')+'.clean.rRNA.fq.gz')
    R2_before_rRNA = path.join(self.CleanDir, path.basename(self.R2).rstrip('.fq.gz')+'.clean.rRNA.fq.gz')
    before_rRNA_sh = " \\\n  ".join([config.get("software", "FQTOOLS"),
                                     config.get("parameter", "FQTOOLS_PAR"),
                                     self.R1, self.R2, self.R1adapter, self.R1adapter,
                                     R1_before_rRNA, R2_before_rRNA,
                                     ])
    rm_rRNA_sh = " \\\n  ".join([config.get("software", "BOWTIE2") + " -p 4",
                                 "-x " + config.get("database", "rRNA"),
                                 "-1 " + R1_before_rRNA,
                                 "-2 " + R2_before_rRNA,
                                 "--un-conc-gz " + path.join(self.CleanDir, self.name+".clean.fq.gz"),
                                 "> /dev/null 2> %s.rRNA.stat" % path.join(self.CleanDir, self.name)
                                 ])
    sh = path.join(self.CleanDir, self.name+"_getCleanData.sh")
    with open(sh, "w") as f:
      f.write(before_rRNA_sh + "\n" + rm_rRNA_sh)
    Job(name=self.name+"_getCleanData", shell=sh, vf="4G", proc="1")
    return self.name
  def Mapping(self, config):
    self.BamDir = path.join(self.outdir, "BAM")
    path.isdir(self.BamDir) or os.mkdir(self.BamDir)
    mapping_sh = " \\\n  ".join([config.get("software", "HISAT2") + " -p 4",
                                 "-x " + config.get("database", "HISAT2_INDEX"),
                                 "-1 " + path.join(self.CleanDir, self.name+".clean.fq.1.gz"),
                                 "-2 " + path.join(self.CleanDir, self.name+".clean.fq.2.gz"),
                                 "2> " + path.join(self.BamDir, self.name + ".bam.stat |"),
                                 config.get("software", "SAMTOOLS") + " view -bS - |",
                                 config.get("software", "SAMTOOLS") + " sort -n - " + path.join(self.BamDir, self.name + ".sort")
                                 ])
    sh = path.join(self.BamDir, self.name+"_Mapping.sh")
    with open(sh, "w") as f:
      f.write(mapping_sh)
    Job(name=self.name+"_Mapping", shell=sh, after=self.name+"_getCleanData", vf="4G", proc="1")
  def ReadCount(self, config):
    self.ReadCountDir = path.join(self.outdir, "ReadCount")
    path.isdir(self.ReadCountDir) or os.mkdir(self.ReadCountDir)
    exonCount = "%s.exon.count.txt" % path.join(self.ReadCountDir, self.name)
    exon_htseq_sh = " \\\n  ".join([config.get("software", "PYTHON3"),
                                    "-m HTSeq.scripts.count -f bam -r name -s no -t exon -i gene_id -m union",
                                    path.join(self.BamDir, self.name + ".sort.bam"),
                                    config.get("database", "GTF"),
                                    ">" + exonCount
                                    ])
    transcript_htseq_sh = " \\\n  ".join([config.get("software", "PYTHON3"),
                                          "-m HTSeq.scripts.count -f bam -r name -s no -t transcript -i gene_id -m union",
                                          path.join(self.BamDir, self.name + ".sort.bam"),
                                          config.get("database", "GTF"),
                                          ">%s.transcript.count.txt" % path.join(self.ReadCountDir, self.name)
                                          ])
    read_distribute_sh = " \\\n  ".join([config.get("software", "Rscript"),
                                         path.join(config.get("software", "BIN"), "count_intron_exon.r"),
                                         "%s.transcript.count.txt" % path.join(self.ReadCountDir, self.name),
                                         exonCount,
                                         "%s.reads.distributes" % path.join(self.ReadCountDir, self.name),
                                         ])
    sh = path.join(self.ReadCountDir, self.name+"_ReadCount.sh")
    with open(sh, "w") as f:
      f.write('\n'.join([exon_htseq_sh, transcript_htseq_sh, read_distribute_sh]))
    Job(name=self.name+"_ReadCount", shell=sh, after=self.name+"_Mapping", vf="4G", proc="1")
    return self.name + ":" + exonCount
  @staticmethod
  def qcSum(qcdir, sample):
    qcSum_sh = ' '.join([config.get("software", "PYTHON3"),
                         path.join(config.get("software", "BIN"), "qcSum.py"),
                         path.abspath(qcdir)] + sample)
    sh = path.join(qcdir, "qcSum.sh")
    after = [i+"_ReadCount" for i in sample]
    with open(sh, "w") as f:
      f.write(qcSum_sh)
    Job(name="qcSum", shell=sh, after=after, vf="1G", proc="1")
def GeneratorQC(sample_info, config, qcdir):
  sampleName = set()
  for infoLine in itertools.dropwhile(lambda x: x.startswith('#'), sample_info):
    if infoLine.strip()=='':
      continue
    sample = infoLine.split()
    if sample[0] in sampleName:
      continue
    else:
      sampleName.add(sample[0])
    qc = QC(fqdir=sample[2], outdir=qcdir, name=sample[0])
    qc.getCleanData(config)
    qc.Mapping(config)
    exonCount = qc.ReadCount(config)
    yield (sample[0], sample[1], exonCount)

class DiffExpr:
  def __init__(self, readCountMatrix, group, compare, venn=None, outdir='.'):
    self.readCountMatrix = readCountMatrix
    self.group = group
    self.compare = compare
    self.venn = venn
    self.outdir = path.abspath(outdir)
    path.isdir(self.outdir) or os.mkdir(self.outdir)
  @staticmethod
  def getSumReadCount(outdir, readcount, config):
    outdir = path.abspath(outdir)
    outfile = path.join(outdir, "ReadCountMatrix.xls")
    ReadCountMatrix_sh = " \\\n  ".join([path.join(config.get("software", "Rscript")),
                                         path.join(config.get("software", "BIN"), "getReadCountMatrix.r"),
                                         outdir] + readcount
                                       )
    sh = path.join(outdir, "getReadCountMatrix.sh")
    with open(sh, "w") as f:
      f.write(ReadCountMatrix_sh)
    after = map(lambda aa: aa.split(':')[0]+"_ReadCount", readcount)
    Job(name="getReadCountMatrix", shell=sh, after=after, vf="2G", proc="1")
    return outfile
  def rundiff_AvsB(self, A, B, config):
    AvsB_outdir = path.join(self.outdir, ''.join([A, 'vs', B]))
    path.isdir(AvsB_outdir) or os.mkdir(AvsB_outdir)
    AvsB = A + ':' + ','.join(self.group[A]) + ' ' + B + ':' + ','.join(self.group[B])
    diff_AvsB_sh = " \\\n  ".join([path.join(config.get("software", "Rscript")),
                                   path.join(config.get("software", "BIN"), "diff_A_vs_B.r"),
                                   self.readCountMatrix, AvsB,
                                   path.join(AvsB_outdir, A+'vs'+B+"_result.xls")
                                  ])
    sh = path.join(AvsB_outdir, A+'vs'+B+'.sh')
    with open(sh, "w") as f:
      f.write(diff_AvsB_sh)
    Job(name="diff_"+A+'vs'+B, shell=sh, after="getReadCountMatrix", vf="2G", proc="1")
  def rundiff(self, config):
    for i in self.compare:
      A, B = i.split(':')
      self.rundiff_AvsB(A, B, config)

if __name__=='__main__':
  parser = argparse.ArgumentParser(description="transcriptome pipline")
  parser.add_argument("--sample_info", "-s", required=True, help="sample information file")
  parser.add_argument("--compare", "-c", required=False, nargs="*", default=None, help="comparable group, : as separator")
  parser.add_argument("--venn", "-vn", required=False, default=None, help="venn information, for example, 2:3:4, int represent index of --compare")
  parser.add_argument("--config", "-con", required=False, default="/USCIMD/usr/wangjianshou/tr/bin/config", help="config file path")
  parser.add_argument("--out", "-o", required=False, default=".", help="output directory, default is current directory")

  args = parser.parse_args()
  config = configparser.ConfigParser()
  config.read(args.config)
  qcdir = path.abspath(path.join(args.out, "QC"))
  path.isdir(qcdir) or os.mkdir(qcdir)
  group = collections.defaultdict(list)
  readcount = []
  with open(args.sample_info, "rt") as f:
    for i in GeneratorQC(f, config, qcdir):
      group[i[1]].append(i[0])
      readcount.append(i[2])
  sample = [i.split(':')[0] for i in readcount]
  QC.qcSum(qcdir, sample)
  outdir = path.abspath(args.out)
  diffDir = path.join(outdir, "DiffExpr")
  path.isdir(diffDir) or os.mkdir(diffDir)
  ReadCountMatrix = DiffExpr.getSumReadCount(diffDir, readcount, config)
  if args.compare is not None:
    diff = DiffExpr(ReadCountMatrix, group, args.compare, outdir=diffDir)
    diff.rundiff(config)
  JobDir = path.join(outdir, "Job")
  path.isdir(JobDir) or os.mkdir(JobDir)
  Job = iter(Job)
  with open(path.join(JobDir, "Trans.JOB"), "w") as f:
    f.write('\n'.join(Job) + '\n\n')
    f.write('\n'.join(Job.order).strip())


