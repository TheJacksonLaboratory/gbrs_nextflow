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
     .filter { it[0] =~ /^[0-9]/}
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

/*
process bowtie{
  publishDir path:params.outdir, mode:'copy', pattern:"*.{log,fastq}"
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

*/
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

/* 
 * run STAR + RSEM to be consistent with Mike Saul
 */
process runSTAR{
  publishDir path: params.outdir, mode: 'copy'
  label 'STAR'
  label 'high_mem'
  tag "$id"
  
  input:
    tuple id, file (reads1), file (reads2) from fastq_ch_star   
    file staridx from star_index.collect()

  output:
    tuple id, file ("${id}.Aligned.toTranscriptome.out.bam") into starout

  script:
  """
  STAR --runThreadN ${task.cpus} \\
  --genomeDir . \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --quantTranscriptomeBAMcompression -1 \
  --quantTranscriptomeBan IndelSoftclipSingleend \
  --outFileNamePrefix ${id}. \
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
    file "${id}_RSEM_count.genes.results" into rsem_genes_count2, RSEM_genes
    val id into RSEM_ids
  script:
  """
  rsem-calculate-expression --bam --paired-end -p ${task.cpus} \\
  --estimate-rspd --append-names \\
  $bamf \
  $fasta \
  ${id}_RSEM_count
  """
}

process RSEM_table{
  publishDir path:params.outdir, mode:'copy'
  label 'RSEM'
  input:
    file genes from rsem_genes_count2.collect()

  output:
    file "RSEM_genes_table.txt" into rsem_table

  script:
  """
  rsem-generate-data-matrix $genes >  RSEM_genes_table.txt
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
/*

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

*/
/*
 * Run the eQTL:
 *   - Read the RSEM results with txImport
 *   - Use DESeq2 to generate rlog transformation of the counts
 *   - Run eQTL with the rlog transformations. Use rdata from the qtlviewer (Striatum) for genotypes
 *   
 * TODO: Add QC GBRS match genotypes
 */
process prepeQTL{
  publishDir path:params.outdir, mode: 'copy'
  label 'R'
  label 'const_mem'
  input:
    file rsems from RSEM_genes.collect()
    file metadata from QTL_metadata
    file geno from genopCh

  output:
    file "Cecum_data_prep.rdata" into readydat, readydat2
    file "tsv_count_table.csv" into count_tbl
    file "raw_count_table.csv" into raw_tbl
  script:
  """
  #!/usr/bin/env Rscript

  library(tximport)
  library(DESeq2)
  library(qtl2)
  library(qtl2convert)
  load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
  snps = GM_snps
  # Read the apr 
  load("$geno")
  # Change rownames to match mouse ID
  for (i in 1:20){ rownames(apr[[i]]) = sapply(sapply(rownames(apr[[i]]), strsplit, "_"), "[[", 6)}
  chrorder=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
  snps <- snps[snps\$chr %in% chrorder, ]
  # Leave markers that are in both snps and apr
  aprmar <- c()
  for (i in 1:20) aprmar <- c(aprmar, dimnames(apr[[i]])[[3]])
  valmar <- intersect(snps\$marker, aprmar)
  snps <- snps[snps\$marker %in% valmar,]
  for (i in 1:20){apr[[i]] = apr[[i]][,,intersect(valmar, dimnames(apr[[i]])[[3]])]}
  map = map_df_to_list(map = snps, pos_column="pos")
  pmap = map_df_to_list(map = snps, pos_column="cM")
  fnames <- strsplit("${rsems}", " ")[[1]]
  metadata <- read.csv("$metadata", stringsAsFactors = FALSE, row.names = NULL)
  metadata\$id <- sapply(sapply(metadata\$id, strsplit, "_"), "[[", 6)
  metadata <- metadata[!duplicated(metadata\$id),]
  rownames(metadata) <- metadata\$id
  names(fnames) <- sapply(sapply(fnames, strsplit, "_"), "[[", 1)
  miceid <- intersect(rownames(metadata), intersect(names(fnames), rownames(apr[[1]])))
  fnames <- fnames[miceid]
  txi.rsem <- tximport(fnames, type="rsem", txIn=FALSE, txOut=FALSE)
  # remove genes with zero counts
  ing <- (rowMin(txi.rsem\$length)>0) & (rowSums(txi.rsem\$counts >= ${params.minreads}) >= ${params.minsamples} * dim(txi.rsem\$counts)[2])
  txi.rsem\$length <- txi.rsem\$length[ing,,drop=F]
  txi.rsem\$counts <- txi.rsem\$counts[ing,,drop=F]
  txi.rsem\$abundance <- txi.rsem\$abundance[ing,,drop=F]
  # import into DESeq2
  sampleTable <- metadata[miceid, "sex", drop=F]
  sampleTable\$sex <- factor(sampleTable\$sex)
  dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~sex)
  #dds <- DESeq(dds)
  ctable <- as.data.frame(assay(vst(dds)))

  # Leave only mouse IDs in the dataset
  for (i in 1:20){apr[[i]] = apr[[i]][miceid,,]}
  # calculate kinship matrix with loco
  K = calc_kinship(probs = apr, type = "loco", use_allele_probs = TRUE)
  K.overall = calc_kinship(probs = apr, type = "overall", use_allele_probs = TRUE)
  cbout <- ctable
  cbout\$gene <- rownames(cbout)
  write.csv(cbout, file="tsv_count_table.csv", quote=FALSE, row.names = FALSE)
  write.csv(txi.rsem\$counts, file="raw_count_table.csv", quote=FALSE)
  ctable <- base::t(ctable)
  # Generate a new dataset for rdata
  save(txi.rsem, ctable, apr, metadata, K, K.overall, map, pmap, snps, file="Cecum_data_prep.rdata")
  #dataset.DO_cecum_416 <- list(annot.mrna = dataset.DO_Striatum_416\$annot.mrna,
                               
  print(fnames)
  print(txi.rsem)
  """
}

/*
 * Run each individual gene eQTL
 */
count_tbl.splitCsv(header:true, by:100)
     .map {row -> row.gene}
     .set { genes }

process rungene{
  publishDir path:params.outdir, mode:'copy'
  label 'R'
  label 'single_cpu'
  input:
    each file (rdata) from readydat
    val gene from genes

  output:
    file "*.pdf" optional true into pdfs
    file "${gene[0]}_eQTL_results.csv" into id_eQTL
    file "${gene[0]}_h2g.csv" into h2g_gene
    file "${gene[0]}_scan1perm.csv" into gene_perm
  script:
  """
  #!/usr/bin/env Rscript
  library("qtl2")
  library("parallel")
  # Includes metadata, apr, K, and ctable - each gene is column, rows are individuals
  load("$rdata")
  genes <- strsplit(gsub(",|\\\\[|\\\\]", "","${gene}"), " ")[[1]]
  covar <- metadata[,c("sex", "ngen")]
  covar\$sex <- (covar\$sex=="M")+0
  for (gene in genes){
    qtl_cis_i = scan1(genoprobs = apr,
                      pheno = ctable[,gene, drop = FALSE],
                      kinship = K,
                      addcovar = covar,
                      cores = ${task.cpus})
    qtl_peaks_7 = find_peaks(qtl_cis_i, map=map,
                             threshold=7)
    # Take only the maximal peak for each gene. This is done for statistical reasons
    # Since the empirical p-value use only the top score
    qtl_peaks_7 <- qtl_peaks_7[which.max(qtl_peaks_7\$lod),]
    h2g = est_herit(ctable[,gene, drop=FALSE], K.overall, covar)
    hout <- c(gene, h2g[1], attr(h2g, "sample_size"), attr(h2g, "log10lik"))
    if (gene == genes[1]){
      all_qtl7 <- qtl_peaks_7
      allh2g <- hout
      qtl_perm_i = scan1perm(genoprobs = apr,
                             pheno = ctable[,gene, drop = FALSE],
                             kinship = K,         
                             addcovar = covar,
                             cores = ${task.cpus},
                             n_perm = 100)
    }else{    
      all_qtl7 <- rbind(all_qtl7, qtl_peaks_7)
      allh2g <- rbind(allh2g, hout)
    }

  }
  write.table(all_qtl7, file="${gene[0]}_eQTL_results.csv", sep=",", quote=F, col.names=F, row.names=F)
  write.table(allh2g, paste0(genes[1], "_h2g.csv"), col.names=F, row.names=F, sep=",", quote=F)
  write.table(qtl_perm_i, paste0(genes[1], "_scan1perm.csv"), col.names=F, row.names=F, sep=",", quote=F)  
  save(qtl_perm_i, file = "${gene[0]}_qtl_out.rdata")
  print(allh2g)
  print(all_qtl7)
  """
}
h2g_gene.collectFile().set { h2g_all }
gene_perm.collectFile().set { perm_all }
process collect_h2g{
  publishDir path:params.outdir, mode:'copy'
  input:
    file h2g from h2g_all.collect()
    file perm from perm_all.collect()

  output:
    file "all_genes_Cecum_h2g_estimates.csv" into h2g_af
    file "all_genes_scan1perm_results.csv" into pall

  script:
  """
  echo "gene,h2g,sample_size,loglike" > all_genes_Cecum_h2g_estimates.csv
  cat ${h2g} >> all_genes_Cecum_h2g_estimates.csv
  cat ${perm} > all_genes_scan1perm_results.csv
  """
}

process eQTL_rdata{
  publishDir path:params.outdir, mode:'copy'
  label 'R'
  input:
    file rdata from readydat2
    file id_res from id_eQTL.collect()
    file permall from pall
  output:
    file "qtlviewer_DO_Cecum_416_*.RData" into viewerdata
    file "all_genes_eQTL_peaks_results.csv" into alleqtl
  script:
  """
  #!/usr/bin/env Rscript
  # Use biomaRt to get gene info
  library(biomaRt)
  library(parallel)
  library(doParallel)
  library(tibble)
  load("${rdata}")
  
  # Load the empirical shuffled LOD scores
  shuflod <- read.csv("${permall}", header=F)\$V1
  getpv <- function(LOD) {
    (sum(shuflod >= LOD)+1)/(length(shuflod)+1)
  }
  maRt = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               host = "dec2015.archive.ensembl.org",
               dataset = "mmusculus_gene_ensembl")
  maRt_filter = "ensembl_gene_id"
  maRt_attributes = c("mgi_symbol","chromosome_name","start_position",
                      "end_position","strand","ensembl_gene_id")
  colnames(ctable) <- sapply(sapply(colnames(ctable), strsplit, "_"), "[[", 1)
  eQTL_maRt = getBM(maRt_attributes, maRt_filter, colnames(ctable), maRt)
  eQTL_maRt\$marker_start = ifelse(sign(eQTL_maRt\$strand) == 1,
                                  eQTL_maRt\$start_position,
                                  eQTL_maRt\$end_position) 
  annotation <- data.frame(gene.id = eQTL_maRt\$ensembl_gene_id,
                           symbol  = eQTL_maRt\$mgi_symbol,
                           chr     = eQTL_maRt\$chromosome_name,
                           start   = round(eQTL_maRt\$start_position/10^6, digits=1),
                           end     = round(eQTL_maRt\$end_position/10^6, digits=1),
                           strand  = eQTL_maRt\$strand, stringsAsFactors = F)
  annotation\$middle <- round((annotation\$start + annotation\$end)/2, digits=1)
  annotation <- annotation[annotation\$chr %in% c(1:19,"X"),]
  annotation <- annotation[!duplicated(annotation\$gene.id),]
  rownames(annotation) <- annotation\$gene.id
  #nearest marker
  cl <- makeCluster(${task.cpus})
  registerDoParallel(cl)
  snps <- snps[snps\$chr %in% c(1:19,"X"),]
  idx <- foreach(i=1:nrow(annotation), .combine='c') %dopar% {
    dist.to <- abs(snps\$pos - annotation\$middle[i])
    min.dist <- min(dist.to[snps\$chr == annotation\$chr[i]])
    which(snps\$chr == annotation\$chr[i] & dist.to==min.dist)[1]
  }
  annotation\$nearest.marker.id <- snps[idx,"marker"]
  stopCluster(cl)

  gnames <- intersect(colnames(ctable), annotation\$gene.id)
  ctable <- ctable[, gnames]
  annotation <- annotation[gnames,]
  
  covar <- metadata[rownames(ctable),c("sex", "ngen")]
  covar\$sex <- (covar\$sex=="M")+0
  cv <- covar
  cv\$sex <- factor(cv\$sex)
  cv\$ngen <- factor(cv\$ngen)
  covar.matrix <- model.matrix(~sex + ngen, cv)[, -1]
  covar\$mouse.id <- rownames(covar)

  #covar.info
  covar.info <- tibble(sample.column   = colnames(covar[,c("sex", "ngen")]),
                           display.name    = c("Sex", "Generation"),
                           interactive     = rep(FALSE, 2),
                           primary         = c(TRUE, TRUE),
                           lod.peaks       = c(NA,NA), stringsAsFactors = FALSE)

  #data
  data <- ctable # you need the expression profile matrix
  #datatype
  datatype <- "mRNA"

  #display.name
  display.name <- "DrugsNaiveCecum_DO"


  # Read eQTLs results
  qfiles = strsplit("${id_res}", " ")[[1]] 
  alleqtls <- NULL
  for (qf in qfiles){
    qr <- read.table(qf, col.names = c("X", "gene", "chr", "pos", "LOD"), stringsAsFactors = FALSE, sep = ",")[,-1]
    alleqtls <- rbind(alleqtls, qr)
  }
  alleqtls\$pvalue <- sapply(alleqtls\$LOD, getpv)
  alleqtls\$qvalue <- p.adjust(alleqtls\$pvalue, method="BH")
  alleqtls <- alleqtls[alleqtls\$qvalue < 0.1, ]
  alleqtls\$gene <- sapply(sapply(alleqtls\$gene, strsplit, "_"), "[[", 1)
  alleqtls <- merge(alleqtls, snps[!duplicated(paste(snps\$chr, snps\$pos, sep="-")),], all.x = TRUE, all.y = FALSE, by = c("chr", "pos"))
  alleqtls <- merge(alleqtls, annotation, by.x = "gene", by.y = "gene.id", suffixes=c("", ".gene"))
  alleqtls\$cis <- alleqtls\$chr == alleqtls\$chr.gene & alleqtls\$start >= alleqtls\$pos-10 & alleqtls\$end <= alleqtls\$pos+10
  write.csv(alleqtls, "all_genes_eQTL_peaks_results.csv")
  lod.peaks <- tibble(gene.id = alleqtls\$gene, marker.id = alleqtls\$marker, lod = alleqtls\$LOD)

  #dataset list
  assign(paste0("dataset.","DO_Cecum_416"),list(annot.mrna    = as_tibble(annotation),
                                            annot.samples      = covar,
                                            covar.matrix       = covar.matrix,
                                            covar.info         = covar.info, 
                                            data               = as.matrix(ctable),      
                                            datatype           = datatype, 
                                            display.name       = display.name, 
                                            lod.peaks          = list(additive=lod.peaks)))
  genoprobs = apr
  markers = tibble(marker.id = snps\$marker, chr = snps\$chr, pos = snps\$pos)
  ensembl.version = 83
  save(genoprobs, 
       K, 
       map, 
       markers, 
       dataset.DO_Cecum_416, 
       ensembl.version,
       file = paste0("qtlviewer_DO_Cecum_416_06302020",".RData")
  )



  
  """
}
