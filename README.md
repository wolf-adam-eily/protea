# protea

<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
<li><a href="#First_Point_Header">1 Overview and directory layout</>
<li><a href="#Second_Point_Header">2 Previous steps taken and paths</a></li>
<li><a href="#Third_Point_Header">3 Quality control using sickle</a></li>
<li><a href="#Fourth_Point_Header">4 Aligning reads to a genome using hisat2</a></li>
<li><a href="#Fifth_Point_Header">5 Predicting gene models with BRAKER</a></li>
<li><a href="#Sixth_Point_Header">6 Pairwise differential expression with counts in R with DESeq2</a></li>
	<ol><li><a href="#types_of_plots">1 Common plots for differential expression analysis</a></li>
		<li><a href="#using_deseq2">2 Using DESeq2</a></li></ol>
<li><a href="#EnTAP">7 EnTAP: Functional Annotation for Genomes</a></li>
 <li><a href="#Integration">8 Integrating the DE Results with the Annotation Results</a></li>
<li><a href="#Citation">Citations</a></li>
</ul>
</div>

<h2 id="Second_Point_Header">Previous steps taken</h2>

The Protea transcriptome was assembled and masked (masked assembly located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/masked_genome`). The  `repeatmodeler` output is located at `/UCHC/LABS/Wegrzyn/proteaBraker/protea/repeatmodeler`.

<h2 id="Third_Point_Header">Quality control using sickle</h2>

The raw reads were paired-end. The forward reads are located at `/UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allforward.fastq`. The reverse reads are located at `/UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allreverse.fastq`. 

The reads were trimmed using sickle with the following code:

<pre style="color: silver; background: black;">module load sickle
sickle pe -f /UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allforward.fastq \
-r /UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allreverse.fastq -t sanger \
-o trimmed_output_RIV21_UCBdata_allforward.fastq -p trimmed_output_RIV21_UCBdata_allreverse.fastq \
-s trimmed_singles_RIV21_UCBdata_all.fastq \
-q 30 \
-l 40</pre>

The output is broken down as follows:

Trimmed forward reads: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/trim_data/trimmed_output_RIV21_UCBdata_allforward.fastq`

Trimmed reverse reads: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/trim_data/trimmed_output_RIV21_UCBdata_allreverse.fastq`

Trimmed singles: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/trimmed_singles_RIV21_UCBdata_all.fastq`

These are the stats from sickle:
<pre style="color: silver; background: black;">
PE forward file: /UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allforward.fastq
PE reverse file: /UCHC/LABS/CBC/vj_projects/protea_repens/raw_reads/RIV21_UCBdata_allreverse.fastq

Total input FastQ records: 403479910 (201739955 pairs)

FastQ paired records kept: 238570808 (119285404 pairs)
FastQ single records kept: 52593028 (from PE1: 48353068, from PE2: 4239960)
FastQ paired records discarded: 59723046 (29861523 pairs)
FastQ single records discarded: 52593028 (from PE1: 4239960, from PE2: 48353068)</pre>

<h2 id="Fourth_Point_Header">Aligning the reads with hisat2</h2>

First, using the masked transcriptome at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/masked_genome`, an index was created with the following code:

<pre style="color: silver; background: black;">
module load hisat2
#hisat2-build genome.fa.masked.filtered  protea_index</pre>

The index is located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/indexing_and_alignment/index`.

Next, the trimmed reads were aligned using the following code:

<pre style="color: silver; background: black;">
module load hisat2
hisat2 -x protea_index -1 trimmed_output_RIV21_UCBdata_allforward.fastq -2 trimmed_output_RIV21_UCBdata_allreverse.fastq -S RIV21_UCBdata_all.sam -p 25</pre>

The singles were omitted from the rest of the analysis.

The aligned reads are located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/indexing_and_alignment/alignment/RIV21_UCBdata_all.sam`.

Next, the aligned reads were sorted with the following code:

<pre style="color: silver; background: black;">module load samtools

samtools sort -@ 8 -n -o sorted_reads.bam RIV21_UCBdata_all.sam</pre>

The sorted reads are located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/alignment_sorting/sorted_reads.bam`.

<h2 id="Fifth_Point_Header">Predicting gene models with BRAKER</h2>

After sorting the alignments, gene models were predicted with BRAKER using those alignments. This was done with the following code:

<pre style="color: silver; background: black;">
module load BRAKER/2.0.5
module load bamtools/2.4.1
export AUGUSTUS_CONFIG_PATH=$HOME/3.2.3/config
export TMPDIR=/home/CAM/$USER/tmp/
export BAMTOOLS_PATH=/isg/shared/apps/bamtools/2.4.1/bin/
export GENEMARK_PATH=/UCHC/LABS/Wegrzyn/local_software/gm_et_linux_64/gmes_petap/
module load perl/5.24.0
export PERL5LIB=/UCHC/LABS/Wegrzyn/perl5/lib/perl5/
export PERLINC=/UCHC/LABS/Wegrzyn/perl5/lib/perl5/

braker.pl --cores 16 --species=protea --softmasking 1 --GENEMARK_PATH=/UCHC/LABS/Wegrzyn/local_software/gm_et_linux_64/gmes_petap/ \
--genome=/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/masked_genome/genome.masked.filtered.fa \
--bam /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/alignment_sorting/sorted_reads.bam --gff3 </pre>

The output is located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gene_modeling_with_BRAKER/braker/protea`.

Of note are a few particular files from the output:

The predicted amino acid sequences of the models: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa`

The corresponding gff3 file to the predicted AA sequences: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa`
