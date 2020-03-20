#!/usr/bin/env nextflow

params.outdir = "."
params.datadir = "."
params.gbrs_data = "/projects/churchill-lab/projects/CubeProject/GBRS_pipeline/AllSuppFiles/R84-REL1505"
params.model = 4

def getLibraryId( prefix ){
  // Return the ID number
  prefix.split("_")[0]
}

// Gather the pairs of R1/R2 accordsing to mouse ID
Channel
     .fromFilePairs(params.datadir + '/**R{1,2}*.fastq.gz', flat: true)
     .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
     .groupTuple() 
     .filter { it[0] =~ /^[0-9]/}
     .into { fastq_ch; fastq_ch_sal; fastq_ch_star }

// Read the metadata file with Sex and generation
Channel
     .fromPath(params.metadata)
     .splitCsv(header:true)
     .map {row -> tuple (row."Mouse.ID", row.Sex, row.Generation)}
     .set { metadata }

process bowtie{
  publishDir path:params.outdir, mode:'copy', pattern:"*.log"
  label 'bowtie1'
  label 'high_mem'
  input:
    tuple id, file(reads1), file(reads2) from fastq_ch

  output:
    tuple id, file ("${id}.1.bam"), file ("${id}.2.bam") into bams
    file "${id}.1.log" into log1
    file "${id}.2.log" into log2
  script:
  """
  zcat $reads1 | bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.1.log |samtools view -bS - > ${id}.1.bam 
  zcat $reads2 | bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.2.log |samtools view -bS - > ${id}.2.bam  
  """
}

process bamtoemase{
  label 'gbrs'
  label 'high_mem'
  input:
    tuple id, file (align1), file (align2) from bams
  
  output:
    tuple id, file ("${id}.merged_compressed.h5") into aln_compressed

  script:
  """
  gbrs bam2emase -i $align1 \
                 -m ${params.gbrs_data}/ref.transcripts.info \
                 -s A,B,C,D,E,F,G,H \
                 -o ${id}.emase1.h5

  gbrs compress -i ${id}.emase1.h5 \
                -o ${id}.compressed.emase1.h5

  gbrs bam2emase -i $align2 \
                 -m ${params.gbrs_data}/ref.transcripts.info \
                 -s A,B,C,D,E,F,G,H \
                 -o ${id}.emase2.h5

  gbrs compress -i ${id}.emase2.h5 \
                -o ${id}.compressed.emase2.h5

  gbrs compress -i ${id}.compressed.emase1.h5,${id}.compressed.emase2.h5\
                -o ${id}.merged_compressed.h5
  """
}

// joins the metadata with the merged into one channel
metadata.join(aln_compressed).into {compm1; compm2}
process quantify{
  publishDir path:params.outdir, mode:'copy'
  label 'gbrs'
  label 'high_mem'
  input:
    tuple id, sex, generation, file (comp) from compm1
    val model from params.model
    env GBRS_DATA from params.gbrs_data
  output:
    file "*" into publish
    file "${id}.gbrs.interpolated.genoprobs.npz" into genoprobs
    tuple id, file("${id}.multiway.isoforms.tpm") into genes_tpm
  script:
  """
  gbrs quantify -i ${comp} \
                -g ${params.gbrs_data}/ref.gene2transcripts.tsv \
                -L ${params.gbrs_data}/gbrs.hybridized.targets.info \
                -M ${model}  --report-alignment-counts -o ${id}

  gbrs reconstruct -e ${id}.multiway.genes.tpm \
                  -t ${params.gbrs_data}/tranprob.DO.${generation}.${sex}.npz \
                  -x ${params.gbrs_data}/avecs.npz \
                  -g ${params.gbrs_data}/ref.gene_pos.ordered.npz -o ${id}

  gbrs quantify -i ${comp} \
                -G ${id}.genotypes.tsv \
                -g ${params.gbrs_data}/ref.gene2transcripts.tsv \
                -L ${params.gbrs_data}/gbrs.hybridized.targets.info \
                -M ${model}  --report-alignment-counts -o ${id}

  gbrs interpolate -i ${id}.genoprobs.npz \
               -g ${params.gbrs_data}/ref.genome_grid.69k.noYnoMT_KBEdit.txt \
               -p ${params.gbrs_data}/ref.gene_pos.ordered_0.1.6.npz \
               -o ${id}.gbrs.interpolated.genoprobs.npz

  gbrs plot -i ${id}.gbrs.interpolated.genoprobs.npz \
               -o ${id}.gbrs.plotted.genome.pdf \
               -n ${id}
  """

}

process exportgeno{
  publishDir path:params.outdir, mode:'copy'
  label 'export'
  input:
    file(genop) from genoprobs
  output:
    file "*" into export_out
  script:
  """
  export-genoprob-file -i ${genop} \
				-s A,B,C,D,E,F,G,H \
				-g ${params.gbrs_data}/ref.genome_grid.69k.noYnoMT_KBEdit.txt
  """
}


Channel
  .fromPath(params.gbrs_data + "/suppFiles/C57BL6J/C57BL6J.fa")
  .into{genome_fa; genome_fa2; star_genome; rsem_genome}
Channel
  .fromPath(params.gbrs_data + "/suppFiles/C57BL6J/C57BL6J.gtf")
  .into{genome_gtf; trans_gtf; star_gtf; rsem_gtf}


// Generate STAR index
process makeSTARindex{
  label 'high_mem'
  label 'STAR'
  tag "$fasta"
  publishDir path: params.outdir, mode: 'copy'

  input:
    file fasta from star_genome
    file gtf from star_gtf

  output:
    file "star/*" into star_index

  script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
      --runMode genomeGenerate \\
      --runThreadN ${task.cpus} \\
      --sjdbGTFfile $gtf \\
      --sjdbOverhang 149 \\
      --genomeDir star \\
      --genomeFastaFiles $fasta \\
      $avail_mem
    """
 }


process runSTAR{
  label 'STAR'
  label 'high_mem'
  tag "$id"
  
  input:
    tuple id, file (reads1), file (reads2) from fastq_ch_star   
    file staridx from star_index.collect()

  output:
    tuple id, file ("${id}Aligned.toTranscriptome.out.bam") into starout

  script:
  """
  STAR --runThreadN ${task.cpus} \\
  --genomeDir . \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --quantTranscriptomeBAMcompression -1 \
  --quantTranscriptomeBan IndelSoftclipSingleend \
  --outFileNamePrefix $id \
  --outSAMattributes NH HI AS nM \
  --readFilesIn ${reads1.join(",")} ${reads2.join(",")}
  """
}

process RSEMprep{
  label 'RSEM'
  label 'high_mem'
  tag "$fasta"
  
  input:
    file fasta from rsem_genome
    file gtf from rsem_gtf

  output:
    file '*' into rsem_index
    file fasta into rsem_index_name
  
  script:
    """
    rsem-prepare-reference --gtf $gtf $fasta $fasta
    """
}

process RSEMit{
  publishDir path:params.outdir, mode:'copy'
  label 'RSEM'
  label 'high_mem'
  tag "$id"
  input:
    tuple id, file (bamf) from starout
    file fasta from rsem_index_name.collect()
    file idxf from rsem_index.collect()
  output:
    tuple id, file ("${id}_RSEM_count.isoforms.results") into rsem_iso_count
    tuple id, file ("${id}_RSEM_count.genes.results") into rsem_genes_count
  script:
  """
  rsem-calculate-expression --bam --paired-end -p ${task.cpus} \\
  --estimate-rspd --append-names \\
  $bamf \
  $fasta \
  ${id}_RSEM_count
  """
}
/*
 * Salmon get transcripts and index
 */
process transcriptsToFasta {
  tag "$fasta"
  label 'gffread'

  input:
    file fasta from genome_fa
    file gtf from genome_gtf

  output:
    file "transcripts.fa" into ch_fasta_for_salmon_index

  script:
    """
    grep ">" $fasta | cut -d " " -f 1 | sed 's/>//' > chrs.txt
    awk 'BEGIN{ while (getline line < "chrs.txt") {chrs[i++]=line}} \$1 in chrs' $gtf > ${gtf}.2
    gffread -F -w transcripts.fa -g $fasta ${gtf}.2
    """
}


process makeSalmonIndex {
  label 'salmon'
  label 'high_mem'
  tag "$fasta"

  input:
    file fasta from ch_fasta_for_salmon_index
    file genome from genome_fa2

  output:
    file 'salmon_index' into salmon_index

  script:
    def gencode = params.gencode  ? "--gencode" : ""
    """
    cat $fasta $genome > gentrome.fa
    grep "^>" $genome | cut -d " " -f 1 > decoys.txt
    salmon index --threads $task.cpus -t gentrome.fa $gencode -d decoys.txt -i salmon_index -k 21
    """
}


process salmon{
  publishDir path:params.outdir, mode:'copy'
  label 'high_mem'
  label 'salmon'
  tag "$id"

  input:
    tuple id, file(reads1), file(reads2) from fastq_ch_sal
    file(index) from salmon_index.collect()
    
  output:
    tuple id, file("salmon_quant_${id}/quant.sf") into salmon_out

  script:
  """
  salmon quant -l A -p ${task.cpus} --seqBias --gcBias --validateMappings --useEM -i $index -1 $reads1 -2 $reads2 -o salmon_quant_${id}
  """
}


process transtogene{
  label 'gffutils'
  publishDir path:params.outdir, mode:'copy'
  input:
    file trans from trans_gtf
  output:
    file "isoforms_to_genes.txt" into isogenes
  script:
  """
  awk '\$3=="transcript"' $trans |cut -f 9 | cut -d " " -f 2,6 | tr -d "\\";" >isoforms_to_genes.txt
  """
}

genes_tpm.join(salmon_out).join(rsem_genes_count).set{both_out}
process compare_tpm{
  publishDir path:params.outdir, mode:'copy'
  label 'R'
  tag "$id"
  input:
    tuple id, file(emase), file(salmon), file(rsem) from both_out
    file isogenes from isogenes.collect()
  output:
    file "${id}_genes_corr.txt" into isocorr
    file "${id}_salmon_genes.pdf" into corrpdf
    file "${id}_rsem_genes.pdf" into rsempdf
    file "${id}_rsem_salmon_genes.pdf" into rsemsalpdf
    file "${id}_genes_count.csv" into genelist
  script:
  """
  #!/usr/bin/env Rscript
  ems <- read.delim("$emase")
  sal <- read.delim("$salmon")
  rsem <- read.delim("$rsem")
  rsem\$Name <- gsub("_.*", "", rsem\$gene_id)
  rsem\$RSEM <- rsem\$TPM
  iso <- read.delim("$isogenes", sep=" ", header=FALSE)
  colnames(iso) <- c("gene", "Name")
  ems <- data.frame(Name = ems\$locus, emase = rowSums(ems[,2:9]))
  both <- merge(ems, sal[,c("Name", "TPM")], by="Name")
  both <- merge(both, iso, by="Name")
  gboth <- aggregate(both[,2:3], by=list(Name=both\$gene), sum)
  gboth <- merge(gboth, rsem[, c("Name", "RSEM")])
  cval <- cor(log10(gboth\$TPM+1), log10(gboth\$emase+1))
  pdf(paste0("${id}_salmon_genes.pdf"))
  plot(log10(gboth\$TPM+1), log10(gboth\$emase+1), cex=0.5, pch=20, main=paste0("$id", " ", cval), xlab="Salmon", ylab="gbrs") 
  dev.off()
  cval2 <- cor(log10(gboth\$RSEM+1), log10(gboth\$emase+1))
  pdf(paste0("${id}_rsem_genes.pdf"))
  plot(log10(gboth\$RSEM+1), log10(gboth\$emase+1), cex=0.5, pch=20, main=paste0("$id", " ", cval2), xlab="STAR+RSEM", ylab="gbrs") 
  dev.off()
  cval3 <- cor(log10(gboth\$RSEM+1), log10(gboth\$TPM+1))
  pdf(paste0("${id}_rsem_salmon_genes.pdf"))
  plot(log10(gboth\$RSEM+1), log10(gboth\$TPM+1), cex=0.5, pch=20, main=paste0("$id", " ", cval3), xlab="STAR+RSEM", ylab="Salmon") 
  dev.off()
  write.csv(gboth, file="${id}_genes_count.csv")

  x <- data.frame(salmon_corrval=cval, RSEM_corrval=cval2, RSEM_Salmon=cval3)
  row.names(x) <- "$id"
  write.table(x, file="${id}_genes_corr.txt", col.names = F, sep = "\t")
  """
}

process col_corr{
  publishDir path:params.outdir, mode:'copy'
  input:
    file corrs from isocorr.collect()

  output:
    file "all_genes_corr.txt" into allcorr

  script:
  """
  cat $corrs > all_genes_corr.txt
  """
}
