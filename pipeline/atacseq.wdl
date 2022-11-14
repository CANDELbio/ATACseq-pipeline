## Copyright Parker Institute for Cancer Immunotherapy, 2022
## 
## 
## Runs ATAC-Seq peak calling on preprocessed bam files or fastq files.
##
## This version specifically calls peaks on individual samples,
## INSTEAD of doing a pooled peak calling. 
##
## Requirements:
## input_type: "bam" or "fastq_list"
##
## Input:
## - input_fastq_list_r1: comma separated list of fastq files broken up by lane
## - input_fastq_list_r2: 
## OR
## - input_bam: a single bam file, mapped or unmapped
##
## References:
## - Blacklist
## - Bowtie2 References (BTDIX1-4, BT_CHROM/REV1/REV2)
##
## Entity Type: Sample
##
## Licensing :
## This script is released under the Apache 2.0 License
## Note however that the programs it calls may be subject to different licenses.
## Users are responsible for checking that they are authorized to run all programs
## before running this script.
##
## Maintainer:
## Robin Kageyama (rkageyama@parkerici.org)
##
## 

workflow atacseq {

  # Bowtie 2 Index
  File BTIDX_1
  File BTIDX_2
  File BTIDX_3
  File BTIDX_4
  File BT_CHROM
  File BT_REV1
  File BT_REV2

  File BLACKLIST
  File TSS #Transcriptional start site map.

  File pepatac_zip # Required for QC. Zip file containing pepatac directory.

  String bowtie2_docker = "biocontainers/bowtie2"
  String genomics_docker = "pici/genomics" # public
  String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.8.1" # public
  String macs2_docker = "regevlab/macs2"
  String qc_docker = "databio/pepatac:latest"


  # Must be uploaded to the data model, comma separated bucket locations, no quotations. 
  String? input_fastq_list_r1
  String? input_fastq_list_r2

  File? input_bam

  String SampleId

  String input_type

  Float disk_multiplier = 2.5
  Float additional_disk = 20


  if (input_type == "bam"){
    Float bam_size = size(input_bam, "GB")
    call bam2fastq {
      input:
        BAM=input_bam,
        # Needs to be genomics_docker as the GATK_docker samtools does not have bam->fastq
        docker = genomics_docker,
        disk_size = bam_size + (disk_multiplier * bam_size) + additional_disk,
        preemptible_tries = 3,
        basename = SampleId
    }
    Array[File] bam2fastq_fqr1 = [bam2fastq.fqr1]
    Array[File] bam2fastq_fqr2 = [bam2fastq.fqr2]
  }

  if (input_type == "fastq_list"){

    call list_to_filearray as list_to_filearray_r1{
      input:
        str = input_fastq_list_r1
    }

    call list_to_filearray as list_to_filearray_r2{
      input:
        str = input_fastq_list_r2
    }
  }

  Array[File] fqr1 = select_first([list_to_filearray_r1.out, bam2fastq_fqr1])
  Array[File] fqr2 = select_first([list_to_filearray_r2.out, bam2fastq_fqr2])



  # Input File
    call bowtie {
      input:
        fastq_r1 = fqr1,
        fastq_r2 = fqr2,
        BTIDX_1=BTIDX_1,
        BTIDX_2=BTIDX_2,
        BTIDX_3=BTIDX_3,
        BTIDX_4=BTIDX_4,
        BT_CHROM=BT_CHROM,
        BT_REV1=BT_REV1,
        BT_REV2=BT_REV2,
        SampleId = SampleId,
        docker=bowtie2_docker
    }
    
    call sam_to_bam{
      input: 
        SAM=bowtie.sam,
        docker=gatk_docker
    }

    call remove_unmapped_and_mito {
      input:
        BAM=sam_to_bam.bam,
        BAI=sam_to_bam.bai,
        docker=gatk_docker

    }
    call remove_duplicates {
      input:
        BAM=remove_unmapped_and_mito.bam,
        BAI=remove_unmapped_and_mito.bai,
        docker=gatk_docker
     }

    call deblacklist{
      input:
        BAM=remove_duplicates.bam,
        BAI=remove_duplicates.bai,
        docker=gatk_docker,
        BLACKLIST=BLACKLIST
     }

    call qc_fragment_length_distribution{
      input:
        BAM=remove_duplicates.bam,
        BAI=remove_duplicates.bai,
        docker=qc_docker,
        pepatac_zip=pepatac_zip
    }

    call qc_tss{
      input:
        BAM=deblacklist.bam,
        BAI=deblacklist.bai,
        docker=qc_docker,
        TSS=TSS,
        pepatac_zip=pepatac_zip
     }

  call peaks{
    input:
      BAM = deblacklist.bam,
      BAI = deblacklist.bai,
      SampleId = SampleId,
      docker=macs2_docker
  }

  call peak_count{
    input:
      BAM=deblacklist.bam,
      BAI=deblacklist.bai,
      PEAKS=peaks.narrow,
      docker=gatk_docker
  }

  call qc_frip{
    input:
      BAM=deblacklist.bam,
      BAI=deblacklist.bai,
      PEAKS=peaks.narrow,
      docker=gatk_docker
  }

output {
  bowtie.mapstat
  sam_to_bam.bam
  sam_to_bam.bai
  deblacklist.bam
  deblacklist.bai
# make_normalized_bigwig.bigwig

  peaks.narrow
  peaks.xls
  peaks.summits

  qc_tss.tss_pdf
  qc_tss.tss_png
  qc_fragment_length_distribution.fld_pdf
  qc_fragment_length_distribution.fld_stats
  qc_frip.frip
}
}

task list_to_filearray{
  # Convert a list of fastqs to an Array
  String str
  String arr = "{ARR[@]}"
  command <<<
    IFS=',' read -ra ARR <<< "${str}"
    for i in $${arr}; do
      echo "$i" >> out
    done
  >>>

  runtime {
    docker: "ubuntu:latest"
  }
  output{
    Array[String] out = read_lines("out")
  }
}


task bam2fastq {
  # Convert a (paired end sequenced) bam to fastq_r1 and fastq_r2 files
  File BAM
  String basename

  String docker
  Int preemptible_tries
  Float disk_size


  command <<<
    samtools fastq -1 ${basename}_r1.fastq -2 ${basename}_r2.fastq ${BAM}
  >>>

  output {
    File fqr1 = "${basename}_r1.fastq"
    File fqr2 = "${basename}_r2.fastq"
  }

  runtime {
    preemptible : preemptible_tries
    memory : "14 GB"
    cpu : "2"
    disks : "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker : docker
  }
}


task bowtie {
  # Align Array[Fastq] to reference genome with bowtie2. 
  # Returns SAM file and mapstats with alignmnet metrics
  Array[File] fastq_r1
  Array[File] fastq_r2
  File BTIDX_1
  File BTIDX_2
  File BTIDX_3
  File BTIDX_4
  File BT_CHROM
  File BT_REV1
  File BT_REV2
  String SampleId

  String docker

  #String sample_base1=basename(fastq_r1)
  #String sample_base = sub(sample_base1,"\\_R1.*fastq.*","")
  String genome_base=basename(BTIDX_1,".1.bt2")


  command <<<
    # Merge Fastq files
    v=$(dirname ${BTIDX_1})/${genome_base}
    bowtie2 -k1 -N1 -p32 -x "$v" -1 ${sep="," fastq_r1} -2 ${sep="," fastq_r2} -S ${SampleId}.sam 2> ${SampleId}_mapstats.txt;
  >>>
  output {
    File sam = "${SampleId}.sam"
    File mapstat = "${SampleId}_mapstats.txt"
  }


  runtime {
    preemptible: "3"
      memory: "32 GB"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "32"
    }

}



task sam_to_bam{
  # Convert sam to bam and index.
  # Returns bam and bai.
  File SAM
  String sample_base=basename(SAM,".sam")

  String docker

  command <<<
    samtools view -bS ${SAM} | samtools sort - -m 32G -@32 ${sample_base}
    samtools index ${sample_base}.bam
  >>>

  output {
    File bam = "${sample_base}.bam"
    File bai = "${sample_base}.bam.bai"
  }

  
  runtime {
    preemptible: "3"
      memory: "32G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "32"
    }
} 

task remove_unmapped_and_mito{
  # Remove mitochondrial reads from BAM
  # Remove unmapped (flag 0x2)
  # Return Bam_deMito
  File BAM
  File BAI
  String sample_base=basename(BAM,".bam")

  String docker

  command<<<
    samtools idxstats ${BAM} | cut -f 1 | grep -v chrM | xargs samtools view -f 2 -b -@ 32 ${BAM} > ${sample_base}_deMito.bam
    samtools index ${sample_base}_deMito.bam
  >>>

  output{
    File bam = "${sample_base}_deMito.bam"
    File bai = "${sample_base}_deMito.bam.bai"
  }
  runtime {
    preemptible: "3"
      memory: "32G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "32"
    }

}

task remove_duplicates{
  # Remove duplicate reads with GATK.
  # Returns Bam_deMito_dedup
  File BAM
  File BAI
  String sample_base=basename(BAM,".bam")

  String docker

  command<<<
    /gatk/gatk  --java-options "-Xmx3000m" MarkDuplicates --INPUT ${BAM}\
     --OUTPUT ${sample_base}_dedup.bam \
     --METRICS_FILE ${sample_base}_dedupMetrics.txt \
     --VALIDATION_STRINGENCY LENIENT \
     --ASSUME_SORTED true \
     --REMOVE_DUPLICATES true
    samtools index ${sample_base}_dedup.bam
  >>>

  output{
    File bam = "${sample_base}_dedup.bam"
    File bai = "${sample_base}_dedup.bam.bai"
  }
  runtime {
    preemptible: "3"
      memory: "16G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "4"
    }

}

task deblacklist{
  # Remove blacklisted regions.
  # - Convert to bed
  # - Adjust read posistion based on strand (due to TN5 offset)
  # - Subtract blacklisted regions
  # - Convert back to bed.
  #
  # This step removes A LOT of read information because the bam is converted (with loss) 
  # to a bed before converted back into a bam. 
  # Example Read in output bam:
  # J00118:442:HJ2CVBBXX:5:1228:2158:27637/1        16      chr1    9998    255     76M     *       0       0       *       *
  # PE info is /1 at the end of the readname. 
  File BAM
  File BAI
  String sample_base=basename(BAM,".bam")
  File BLACKLIST

  String docker

  command<<<
    bedtools bamtobed -i  ${BAM} | awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 5, $3 + 5, $4, $5, $6; else print $1, $2 - 4, $3 - 4, $4, $5, $6}' > ${sample_base}.bed
    bedtools subtract -A -a ${sample_base}.bed -b ${BLACKLIST} > ${sample_base}_debl.bed
    bedtools bedtobam -i ${sample_base}_debl.bed -g /usr/share/bedtools/genomes/human.hg19.genome | samtools sort - -m $16G  -@ 8 ${sample_base}_debl
    samtools index ${sample_base}_debl.bam
  >>>

  output{
    File bam = "${sample_base}_debl.bam"
    File bai = "${sample_base}_debl.bam.bai"
  }
  runtime {
    preemptible: "3"
      memory: "16G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "8"
    }
}


task peaks{
  # Call Peaks with macs2.
  # Returns:
  # - XLS excel file.
  # - narrowpeak:
  #   - [chr start stop peakname int(-10*log10qvalue) strand(none) fold-change -log10pvalue -log10qvalue summit_offset]
  # - summits.bed
  #   - [chr start start+1 peakname -log10qvalue]
  File BAM
  File BAI
  String SampleId

  String docker

  command <<<
  macs2 callpeak -t ${BAM} -f BAM -g hs -n ${SampleId} -q 0.01
  >>>


  output{
    File xls = "${SampleId}_peaks.xls"
    File narrow = "${SampleId}_peaks.narrowPeak"
    File summits = "${SampleId}_summits.bed"
  }

  runtime {
      memory: "16G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "2"
    }
}

task peak_count{
  # Get the number of reads at each peak. 
  # Returns _deMito_dedup_debl.bed
  # - [chr start stop peakname int(-10*log10qvalue) strand(none) fold-change -log10pvalue -log10qvalue summit_offset
  #       #ofReads-in-peak #of-bases-with-reads #length-of-peak %peak-covered-with-reads]
  File BAM
  File BAI
  File PEAKS
  String sample_base=basename(BAM,".bam")

  String docker

  command <<<
    # Re-sort the bed file to match bam chromosome order.
    samtools idxstats ${BAM} | cut -f 1 > chromosome_order.txt
    bedtools sort -faidx chromosome_order.txt -i ${PEAKS} > resorted_peaks.narrowPeak
    bedtools coverage -sorted -a resorted_peaks.narrowPeak -b ${BAM} > ${sample_base}.bed
  >>>

  output {
    File bed = "${sample_base}.bed"
    File chromosome_order = "chromosome_order.txt"
  }
  runtime {
    preemptible: "3"
      memory: "16G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "4"
    }
}



########### QC #############

task qc_tss{
  File BAM
  File BAI
  String docker
  File TSS
  File pepatac_zip

  String sample_base=basename(BAM,".bam")

  command <<<
    unzip ${pepatac_zip}
    python pepatac-0.8.6/tools/pyTssEnrichment.py -a ${BAM} -b ${TSS} -o ${sample_base}_tss -p ends -e 2000 -u -v -s 4
    Rscript pepatac-0.8.6/tools/PEPATAC_TSSenrichmentPlot.R ${sample_base}_tss
  >>>

  output {
    File tss_pdf = "${sample_base}_tss.pdf"
    File tss_png = "${sample_base}_tss.png"
  }
  runtime {
    preemptible: "3"
    memory: "16G"
    docker: docker
    disks: "local-disk 100 SSD"
    cpu: "4"
  }
}

task qc_fragment_length_distribution {
  File BAM
  File BAI
  String docker
  File pepatac_zip
  
  command <<<
    unzip ${pepatac_zip}
    perl pepatac-0.8.6/tools/fragment_length_dist.pl ${BAM} fld_length  
    sort -n fld_length | uniq -c  > fld_count
    Rscript pepatac-0.8.6/tools/fragment_length_dist.R fld_length fld_count fld.pdf fld_stats.txt
  >>>

  output {
    File fld_stats = "fld_stats.txt"
    File fld_pdf = "fld.pdf"
  }

  runtime {
    preemptible: "3"
    memory: "16G"
    docker: docker
    disks: "local-disk 100 SSD"
    cpu: "2"
  }
}

task qc_frip{
  # Calculate fraction of reads in peaks. 
  File BAM
  File BAI
  File PEAKS

  String docker

  command <<<
    # Re-sort the bed file to match bam chromosome order.
    samtools idxstats ${BAM} | cut -f 1 > chromosome_order.txt
    bedtools sort -faidx chromosome_order.txt -i ${PEAKS} > resorted_peaks.narrowPeak
    
    # For each read, figure out how many peaks it overlaps
    bedtools coverage -sorted -a ${BAM} -b resorted_peaks.narrowPeak > frip.bed
    # Get percent of reads that overlap at least one peak. This prints to stdout.
    awk '{ if ($13 > 0) total += 1 } END { print total/NR }' frip.bed
  >>>

  output{
    String frip = read_string(stdout())
  }

  runtime {
    preemptible: "3"
      memory: "8G"
      docker: docker
      disks: "local-disk 100 SSD"
      cpu: "2"
    }
}










