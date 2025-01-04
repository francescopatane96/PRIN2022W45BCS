#!/bin/bash

# Dichiarazione dell'array
ids=("17_EI" "13_FS_SCS" "10_VG" "8_SC" "7_EHZ")
samples=("17_EI" "13_FS_SCS" "10_VG" "8_SC" "7_EHZ")

# Loop per ogni campione
for i in "${!samples[@]}"; do
  sample="${samples[$i]}"
  # Rimuovi eventuali caratteri non ammessi (se presenti) e valida l'id
  id=$(echo "${ids[$i]}" | tr -cd '[:alnum:]_-')
  
  echo "Executing cellranger count for sample: $sample"
  
  cellranger count \
    --id="$id" \
    --fastqs=/home/francesco/Desktop/mnt2/bcl_fastq/RUN374NS/RUN374NS/outs/fastq_path/HKHVMDRX5/ \
    --sample="$sample" \
    --transcriptome=/home/francesco/Desktop/refdata-gex-GRCh38-2020-A \
    --localcores=23 \
    --localmem=60 
done




