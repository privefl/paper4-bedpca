wd <- setwd("data/ukbb_bed")
for (chr in 1:22) {
  # bed
  system(paste0("./ukbgene cal -c", chr))
  # fam
  system(paste0("./ukbgene cal -c", chr, " -m"))
}
# bim
download.file("https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_bim.tar", "ukb_snp_bim.tar")
untar("ukb_snp_bim.tar")

write(sapply(1:22, function(chr) {
  c(paste0("ukb_cal_chr", chr, "_v2.bed"),
    paste0("ukb_snp_chr", chr, "_v2.bim"),
    paste0("ukb25589_cal_chr", chr, "_v2_s488282.fam"))
}), tmp <- tempfile(), ncolumns = 3)

library(bigsnpr)
snp_plinkQC(
  plink.path = download_plink("."),
  prefix.in = tmp,
  file.type = "--merge-list",
  prefix.out = "ukbb_488282",
  geno = 0.01,
  autosome.only = TRUE
)
# PLINK v1.90b6.10 64-bit (17 Jun 2019)          www.cog-genomics.org/plink/1.9/
# (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to ukbb_488282.log.
# Options in effect:
#   --autosome
#   --geno 0.01
#   --hwe 1e-50
#   --maf 0.01
#   --make-bed
#   --merge-list /tmp/RtmpHoVtFZ/file32bc56c9faf9b
#   --mind 0.1
#   --out ukbb_488282
#
# 513440 MB RAM detected; reserving 256720 MB for main workspace.
# Performing single-pass merge (488377 people, 784256 variants).
# Merged fileset written to ukbb_488282-merge.bed + ukbb_488282-merge.bim +
#   ukbb_488282-merge.fam .
# 784256 variants loaded from .bim file.
# 488377 people (223474 males, 264808 females, 95 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to ukbb_488282.nosex .
# 6 people removed due to missing genotype data (--mind).
# IDs written to ukbb_488282.irem .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 488371 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate in remaining samples is 0.969389.
# 144436 variants removed due to missing genotype data (--geno).
# --hwe: 50742 variants removed due to Hardy-Weinberg exact test.
# 84939 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 504139 variants and 488371 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to ukbb_488282.bed + ukbb_488282.bim + ukbb_488282.fam ... done.

setwd(wd)
