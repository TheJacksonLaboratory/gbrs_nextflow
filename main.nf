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
     .set { fastq_ch }

// Read the metadata file with Sex and generation
Channel
     .fromPath(params.metadata)
     .splitCsv(header:true)
     .map {row -> tuple (row."Mouse.ID", row.Sex, row.Generation)}
     .set { metadata }

process bowtie{
  publishDir path:params.outdir, mode:'copy', pattern:"*.log"
  label 'bowtie1'
  input:
    tuple id, file(reads1), file(reads2) from fastq_ch

  output:
    tuple id, file ("${id}.1.sam"), file ("${id}.2.sam") into sams
    file "${id}.1.log" into log1
    file "${id}.2.log" into log2
  script:
  """
  zcat $reads1 | bowtie -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.1.log > ${id}.1.sam 
  zcat $reads2 | bowtie -q -a --best --strata --sam -v 3 ${params.gbrs_data}/transcripts - 2>${id}.2.log > ${id}.2.sam  
  """
}

process samtobam{
  publishDir path:params.outdir, mode:'copy'
  label 'samtools'
  input:
    tuple id, file (align1), file (align2) from sams

  output:
    tuple id, file ("${align1.baseName}.bam"), file ("${align2.baseName}.bam") into bams

  script:
  """
  samtools view -bS $align1 > ${id}.1.bam
  samtools view -bS $align2 > ${id}.2.bam
  """
}

process bamtoemase{
  label 'gbrs'
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
  input:
    tuple id, sex, generation, file (comp) from compm1
    val model from params.model
    env GBRS_DATA from params.gbrs_data
  output:
    file "*" into publish
    file "${id}.gbrs.interpolated.genoprobs.npz" into genoprobs
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

