Run GBRS (https://github.com/churchill-lab/gbrs) on paired-end RNAseq data using nextflow pipelines. This is adapted to Sumner cluster (slurm, singularity) and uses the jaxreg.jax.org containers repository, follow the User Guide there to add it to your library. To test it install nextflow and run:
```bash
nextflow run TheJacksonLaboratory/gbrs_nextflow --datadir <input_dir> --metadata <metadata_file> -profile singularity,slurm -resume --outdir <outdir_dir>
```
Make sure the input files are in the format: ID_*R{1,2}*.fastq.gz
The metadata file is a csv file that contains the columns: Mouse.ID, Sex and Generation

