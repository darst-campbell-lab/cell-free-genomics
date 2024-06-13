# cell-free-genomics
The Jupyter Notebooks, Python and R scripts in this repo were used to analyze raw fastq data, identify promoters and terminators in the M. tuberculosis genome (i.e. genomic coordinates where transcript ends significantly accumulated), further narrow in on putative TF targets (promoters/terminators that were differentially expressed in the presence of a TF), and generate all subsequent analysis figures. \
<br>
The core cell-free genomics pipeline relies on the following Jupyter notebooks:\
<br>
## CellFreeGenomics_readPreparation 
**General purpose:** Prepare and align sequencing reads to reference genomes.   \
*Note:* separate de-multiplexing and quality-control pipelines were developed for RNA and genomic DNA samples, since these libraries were prepared differently.\
Pipeline steps include:
* Check fastq quality with fastqc  
* De-multiplex fastq files based on inline barcode sequence (code written by Peter Culviner, PhD)
* UMI removal from read sequence and addition to read ID (necessary for subsequent read de-duplication)
* Paired-end read alignment to concatenated reference genome (E. coli/Eco + M. tuberculosis/Mtb)
* Remove PCR and optical duplicates from alignments using umi_tools
* Generate alignments containing transcript end reads only (read 2)
* 2 separate alignments: one for spike genome (E. coli) & one for experimental genome (M. tuberculosis)
<!------>
## CellFreeGenomics_identifyEnrichedEnds 
**General purpose:** Identify TSSs and TTSs in each replicate.\
Pipeline steps include:
* Downsample the Eco and Mtb alignments containing only transcript end reads to equivalent sequencing depths
* Generate single-bp resolution .txt files, containing both transcript end read counts only and total coverage counts at every genomic position (important for later steps)
* Generate bigWig files (separate + and – strand files) as inputs for the nonparametric resampling script
* Call TSSs and TTSs using a nonparametric resampling approach at every position in the genome (developed by Mike Wolfe, PhD). This requires that the NETseq_pause_calling.py script, with dependencies arraytools.py & bwtools.py, is in the working directory.
* Generate consensus TSS/TTS calls for the Eco spike alignments (i.e. TSSs or TTSs found in all three replicates of a given condition)
<!------>
## CellFreeGenomics_thresholdSelection
**General purpose:** Identify putative transcription factor targets (promoters or terminators) and de novo motifs.\
There are two possible pipelines, depending on the experimental design.
<br>
* Pairwise (+/– TF): the main function is `cpmThreshold_motifRecovery_pairwise`. This function iterates through a range of user-provided CPM thresholds to compare motif discovery and number    of differentially-expressed TSSs/TTSs identified at each threshold. Specifically, at each CPM threshold, the function:
    * Identifies high-confidence TSSs and TTSs (those present in all three replicates of a given condition that exceed the CPM threshold and achieve a nonparametric resampling q-value ≤ 0.05)
    * Filters out TSSs and TTSs found in blacklist samples
    * Performs quality-control motif discovery from the high-confidence TSSs and TTSs (for your reference, this doesn't affect the TSSs/TTSs included in subsequent steps)
      * Control sequence input: randomly-selected sequences from the genome
      * Quality-control motif for TSSs: –10 element
      * Quality-control motif for TTSs: U-tract
    * Takes the union of all TSSs/TTSs identified in both conditions (+TF and –TF)
    * Filters out any TSSs/TTSs with significant changes in background transcription levels (total coverage – transcript ends)
    * Performs DESeq2 on remaining TSS/TTSs to identify putative TF targets
    * Performs de novo motif analysis on putative TF targets
      * XSTREME: MEME/STREME for discovery, SEA for enrichment, FIMO for scanning
<br>
* Multifactor (two or more TFs): the main function is `cpmThreshold_multifactor`. A similar procedure is performed as described above to calculate numbers of differentially expressed           TSSs/TTSs at certain thresholds, except that no motif discovery is performed because this is best analyzed using a pairwise design.

Results from each CPM threshold tested are compiled in a .csv file using the `compile_results` function. \
<br>
If relevant (i.e. a de novo motif is discovered), DESeq2 results can be integrated with FIMO motif scanning results using the `integrate_DESeq2_motif` function.


