#  chmod +x 0_cut_adapters.sh

# INPUT: 
# folder with RAW reads
BF_RAW_FOLDER=~/Projects_R/twins_microbiome_pipeline/data_set_bodyfl/raw/raw/
# adapters


#cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq
# -q 20,20 - quality trimming
# -e - error filtering
# -o - output / -p - paired output
# -a - 3’ adapter, forwardPrimer1InverseSequence
# -g - 5’ adapter  (front), forwardPrimer


echo "--- Start cutting adapters with cutadapt version ------"
cutadapt --version

for SAMPLE_F in ${BF_RAW_FOLDER}*R1.fastq;
do
  # form reverse reads
  SAMPLE_R=$(echo ${SAMPLE_F} | sed "s/R1/R2/") 
  echo ${SAMPLE_F}
  echo ${SAMPLE_R}
  
  # form output path
  OUT_F=$(echo ${SAMPLE_F} | sed "s/raw\/raw/raw\/no_primers/") 
  OUT_R=$(echo ${SAMPLE_R} | sed "s/raw\/raw/raw\/no_primers/") 
  echo ${OUT_F}
  echo ${OUT_R}
  
  cutadapt \
  -g "CCAGCAGCYGCGGTAAN" \
  -G "CCGTCAATTCNTTTRAGT" -G "CCGTCAATTTCTTTGAGT"  -G "CCGTCTATTCCTTTGANT" \
  -o ${OUT_F} -p ${OUT_R} \
  ${SAMPLE_F} ${SAMPLE_R}

done
