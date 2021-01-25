#!/usr/bin/env nextflow

params.outdir = "."
params.datadir = "."
params.gbrs_data = "/projects/churchill-lab/projects/CubeProject/GBRS_pipeline/AllSuppFiles/R84-REL1505"
params.genoprobs = "/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/apr_DO2816.RData"
params.covar = "/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/covar.csv"
params.model = 4
params.rdata = "rdata"
params.minreads = 10
params.minsamples = 0.8 // At least 0.8 of the samples have at least 5 reads

def getLibraryId( prefix ){
  // Return the ID number
  prefix.split("_")[0]
}

// Gather the pairs of R1/R2 accordsing to mouse ID
Channel
     .fromFilePairs(params.datadir + '/**R{1,2}*.fastq.gz', flat: true)
     .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
     .groupTuple() 
//     .filter { it[0] =~ /^[0-9]/}
     .into { fastq_ch; fastq_ch_sal; fastq_ch_star }

// Read the metadata file with Sex and generation
Channel
     .fromPath(params.covar)
     .splitCsv(header:true)
     .map {row -> tuple (row.id.split("_")[5], row.sex, "G${row.ngen}")}
     .set { metadata }
// Also keep the metadata as a channel for thr R eQTL
Channel
     .fromPath(params.covar)
     .set { QTL_metadata }
// Generate a channel for rdata
Channel
     .fromPath(params.rdata)
     .set { rdataCh }
//genotypes
Channel
      .fromPath(params.genoprobs)
      .set { genopCh }
/*
 * Run the classical GBRS pipeline with bowtie1
 */


process bowtie{
  publishDir path:params.outdir, mode:'copy'
  label 'bowtie1'
  label 'high_mem'
  input:
    tuple id, file(reads1), file(reads2) from fastq_ch

  output:
    tuple id, file ("${id}.1.bam"), file ("${id}.2.bam") into bams
    file "${id}.combined.R1.fastq" into fqall1
    file "${id}.combined.R2.fastq" into fqall2
    file "${id}.1.log" into log1
    file "${id}.2.log" into log2
  script:
  """
  zcat $reads1 |tee ${id}.combined.R1.fastq | bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.1.log |samtools view -bS - > ${id}.1.bam 
  zcat $reads2 |tee ${id}.combined.R2.fastq |bowtie -p ${task.cpus} -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.2.log |samtools view -bS - > ${id}.2.bam  
  """
}

process bamtoemase{
  publishDir path:params.outdir, mode:'copy'
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
