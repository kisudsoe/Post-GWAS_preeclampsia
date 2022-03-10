#!/bin/bash

# 11/11/2021 VERSION
# This file is maintained by Seungsoo Kim (kisudsoe@gmail.com).

WORK_DIR="/data"
BASE_BED=$WORK_DIR"/candidate_1197.bed"
SH_FILE=$WORK_DIR"/dist_data.sh"
ANN_PATH="/data/db_gwas"

ROAD_FILE=$ANN_PATH"/roadmap_meta.tsv"
REG_DIR=$ANN_PATH"/regulome"
ENCODE_FILE=$ANN_PATH"/wgEncodeRegTfbsClusteredWithCellsV3.bed"


### RUN FUNCTIONS ###
# DO NOT CHANGE BELLOW THIS CODE.

#ANN_GTEX=$WORK_DIR"/gtex_signif_3938.tsv"
export WORK_DIR BASE_BED SH_FILE GWAS_IN ANN_PATH
export RAOD_FILE REG_DIR ENCODE_FILE

printf "\n1. Overlap the UCSC, ENCODE, and Roadmap annotations\n"
printf "  Writing "$SH_FILE"\n"
Rscript postgwas-exe.r \
    --utils bash \
    --base $BASE_BED \
    --out  $WORK_DIR \
    --ann  $ANN_PATH

printf "  Run bedtools commands for genome_dist, Roadmap_dist_15/18/25, and custom_dist. This step is taking a while. "
time cat $SH_FILE | parallel

printf "\n2. Annotations for UCSC, Roadmap, ENCODE, RegulomeDB, lncRNASNP, and GTEx\n\n"
time cat /tsv-analysis.sh | parallel

printf "\n3. Convert BED for UCSC, ENCODE, Regulome, and lncRNASNP annotations\n\n"
time cat /bed-analysis.sh | parallel

printf "\n4. GWAS summary "
printf "4-1. Roadmap 18 Enh_Prom & Encode Tfbs to BED"
RD25ENH=`find $WORK_DIR"/summary_bed" -name "roadmap_18-enh_prom_*.bed"`
ENCODE=`find $WORK_DIR"/summary_bed" -name "snp_encode_tfbs_*.bed"`

Rscript postgwas-exe.r \
    --dbvenn uniset \
    --base $RD25ENH,$ENCODE \
    --out $WORK_DIR"/genome_tsv/enh_prom_tfbs" \
    --out_bed $WORK_DIR"/summary_bed/enh_prom_tfbs" \
    --intersect TRUE \
    > $WORK_DIR"/log-enhancer_tfbs.txt"
printf " - done\n"

printf "4-2. Overlap-1 to BED"
ENH_TFBS=`find $WORK_DIR"/summary_bed" -name "enh_prom_tfbs_*.bed"`
REGULOME=`find $WORK_DIR"/summary_bed" -name "snp_regulome2b_*.bed"`

Rscript postgwas-exe.r \
    --dbvenn uniset \
    --base $ENH_TFBS,$REGULOME \
    --out $WORK_DIR"/genome_tsv/overlap-1" \
    --out_bed $WORK_DIR"/summary_bed/overlap-1" \
    --intersect TRUE \
    > $WORK_DIR"/log-overlap-1.txt"
printf " - done\n"

printf "4-3. Overlap-2 to BED"
ENH_TFBS=`find $WORK_DIR"/summary_bed" -name "enh_prom_tfbs_*.bed"`
GTEX=`find $WORK_DIR"/summary_bed" -name "gtex_eqtl_*.bed"`

Rscript postgwas-exe.r \
    --dbvenn uniset \
    --base $ENH_TFBS,$GTEX \
    --out $WORK_DIR"/genome_tsv/overlap-2" \
    --out_bed $WORK_DIR"/summary_bed/overlap-2" \
    --intersect TRUE \
    > $WORK_DIR"/log-overlap-2.txt"
printf " - done\n"

printf "4-4. Summarize GWAS candidates, Nearest genes, and GTEx eQTLs"
ANN_GTEX=`find $WORK_DIR"/genome_tsv" -name 'gtex_signif_*.tsv'`

Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR"/summary_bed" $WORK_DIR"/gtex_eqtl" $WORK_DIR"/encode_bed" \
  --out    $WORK_DIR"/genome_tsv" \
  --sub_dir FALSE \
  --uni_save FALSE \
  --ann_gwas $GWAS_IN \
  --ann_near $WORK_DIR"/genome_tsv/nearest_gene.tsv" \
  --ann_gtex $ANN_GTEX \
  > $WORK_DIR"/log-gwas_summary.txt"
printf "done\n"

### END FUNCTIONS ###
