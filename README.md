# protea

<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
<li><a href="#First_Point_Header">1 Overview and directory layout</>
<li><a href="#Second_Point_Header">2 Previous steps taken and paths</a></li>
<li><a href="#Third_Point_Header">3 Quality control using sickle</a></li>
<li><a href="#Fourth_Point_Header">4 Aligning reads to a genome using hisat2</a></li>
<li><a href="#Fifth_Point_Header">5 Predicting gene models with BRAKER</a></li>
<li><a href="#Sixth_Point_Header">6 Quality control of gene models using gfacs</a></li>
	<li><a href="#Seventh_Point_Header">7 Functional annotation using EnTAP</a></li>
	<ol><li><a href="#uniprot">Statistics from uniprot database search</a></li>
		<li><a href="#plantfaa">Statistics from ref-seq plant faa 87 database search</a></li>
		<li><a href="#combined">Integrated statistics for i & ii</a></li>
		<li><a href="#taxonomics">Taxonomic breakdown of EnTAP run</a></li></ol>
<li><a href="#Eighth_Point_Header">8 Further statistical breakdown of EnTAP output</a></li>
	<li><a href="#Ninth_Point_Header">9 Final GTF check</a></li>
</ul>
</div>

<h2 id="Second_Point_Header">Previous steps taken</h2>

The Protea transcriptome was assembled and masked (masked assembly located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome`). The  `repeatmodeler` output is located at `/UCHC/LABS/Wegrzyn/proteaBraker/protea/repeatmodeler`.

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

Trimmed forward reads: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/trim_data/trimmed_output_RIV21_UCBdata_allforward.fastq`

Trimmed reverse reads: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/trim_data/trimmed_output_RIV21_UCBdata_allreverse.fastq`

Trimmed singles: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/trimmed_singles_RIV21_UCBdata_all.fastq`

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

First, using the masked transcriptome at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome`, an index was created with the following code:

<pre style="color: silver; background: black;">
module load hisat2
#hisat2-build genome.fa.masked.filtered  protea_index</pre>

The index is located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/indexing_and_alignment/index`.

Next, the trimmed reads were aligned using the following code:

<pre style="color: silver; background: black;">
module load hisat2
hisat2 -x protea_index -1 trimmed_output_RIV21_UCBdata_allforward.fastq -2 trimmed_output_RIV21_UCBdata_allreverse.fastq -S RIV21_UCBdata_all.sam -p 25</pre>

The singles were omitted from the rest of the analysis.

The statistics of the alignment are:
<pre style="color: silver; background: black;">
119285404 reads; of these:
  119285404 (100.00%) were paired; of these:
    48394548 (40.57%) aligned concordantly 0 times
    69075616 (57.91%) aligned concordantly exactly 1 time
    1815240 (1.52%) aligned concordantly >1 times
    ----
    48394548 pairs aligned concordantly 0 times; of these:
      1673017 (3.46%) aligned discordantly 1 time
    ----
    46721531 pairs aligned 0 times concordantly or discordantly; of these:
      93443062 mates make up the pairs; of these:
        70385893 (75.32%) aligned 0 times
        22302298 (23.87%) aligned exactly 1 time
        754871 (0.81%) aligned >1 times
70.50% overall alignment rate
</pre>

The aligned reads are located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/indexing_and_alignment/alignment/RIV21_UCBdata_all.sam`.

Next, the aligned reads were sorted with the following code:

<pre style="color: silver; background: black;">module load samtools

samtools sort -@ 8 -n -o sorted_reads.bam RIV21_UCBdata_all.sam</pre>

The sorted reads are located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/alignment_sorting/sorted_reads.bam`.

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
--genome=/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome/genome.masked.filtered.fa \
--bam /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/alignment_sorting/sorted_reads.bam --gff3 </pre>

The output is located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea`.

Of note are a few particular files from the output:

The predicted amino acid sequences of the models: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa`

The corresponding gff3 file to the predicted AA sequences: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa`

There were 35896 predicted amino acid sequences determined via: `grep -c ">" augustus.hints.aa`.

However, there were only 33270 genes predicted, determined via: `grep -c "gene" augustus.hints.gff3`.

This implies that a few genes have multiple isoforms.

<h2 id="Sixth_Point_Header">Quality control of gene models using gfacs</h2>

The gene models were quality checked using gfacs. There were no flags used to eliminate any partials, including double partials, as the quality control thresholds for the partials are the same to the completes save for the actual partials flags (`--remove-without-start-codon`, `--remove-without-stop-codon`). For better statistics, the gene models were split into multi-exonic and mono-exonic groups using gfacs (with the flags `--remove-multiexonics`, `--remove-monoexonics` on two separate runs). Therefore, no gene appears in both sets and a combination of the two sets compiles all passing gene models. After each gfacs run, the amino acid fasta was run through checking software (included in this github) to separate the models into four categories: double partials, 5p partials, 3p partials, and completes. The double partials were then removed.

MONOEXONICS:
<pre style="color: silver; background: black;">
module load perl/5.24.0
cd /UCHC/LABS/Wegrzyn/gFACs/
perl gFACs.pl \
-f braker_2.05_gff3 \
--statistics \
--statistics-at-every-step \
-p mono \
--splice-rescue \
--unique-genes-only \
--rem-start-introns \
--rem-end-introns \
--splice-table \
--min-exon-size 20 \
--min-intron-size 20 \
--rem-multiexonics \
--get-fasta-without-introns \
--get-fasta-with-introns \
--create-gtf \
--get-protein-fasta \
--fasta /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome/genome.masked.filtered.fa \
-O /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/mono_exonics/ \
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.gff3</pre>

The results are at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/mono_exonics/`.

Here are some statistics (from `*statistics.txt`:
<pre style="color: silver; background: black;">
Number of genes:        6332
Number of positive strand genes:        3210
Number of negative strand genes:        3122

Average size of mono-exonic genes:	657.741
Median size of mono-exonic genes:	432
Largest monoexonic gene:        4955
Smallest monoexonic gene:	74
</pre>

The monoexonics amino acid fasta was run through the checking software with the following results:
<pre style="color: silver; background: black;">
Double Partials		5p Partials		3p Partials	 		Complete Genes	 	Total
83			37			113				6099			6332</pre>

Only the complete genes were taken from this sample. We see agreement.

MULTIEXONICS:
<pre style="color: silver; background: black;">
module load perl/5.24.0
cd /UCHC/LABS/Wegrzyn/gFACs/
perl gFACs.pl \
-f braker_2.05_gff3 \
--statistics \
--statistics-at-every-step \
-p multi \
--splice-rescue \
--unique-genes-only \
--rem-start-introns \
--rem-end-introns \
--splice-table \
--min-exon-size 20 \
--min-intron-size 20 \
--rem-monoexonics \
--get-fasta-without-introns \
--get-fasta-with-introns \
--create-gtf \
--get-protein-fasta \
--fasta /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome/genome.masked.filtered.fa \
-O /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/multi_exonics/ \
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.gff3</pre>

The results are at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/multi_exonics/ `

Here are some statistics (from `*statistics.txt`):
<pre style="color: silver; background: black;">
Number of genes:        17306
Number of positive strand genes:        8579
Number of negative strand genes:        8727

Average overall gene size:	6212.739
Median overall gene size:	3397
Average overall CDS size:	1067.629
Median overall CDS size:        834
Average overall exon size:	203.940
Median overall exon size:	124

Average number of exons per gene:   5.235
Median number of exons per gene:    4
Largest exon:	8028
Smallest exon:	20
Most exons in one gene: 47

Average number of introns per gene: 4.235
Median number of introns per gene:  3
Largest intron: 60462
Smallest intron:        42</pre>

The multiexonics amino acid fasta was run through the checking software with the following results:
<pre style="color: silver; background: black;">
Double Partials		5p Partials		3p Partials	 	Complete Genes	 	Total
147			63			179			16917			17306</pre>

The double partials were removed. We see agreement.

The multiexonic and monoexonic checked faa's and gtfs were combined into `all_genes.faa` and `all_genes.gtf` at the location:
`/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/all_genes`.

Our final table is as follows:
<pre style="color: silver; background: black;">
WITH_MONO_EXONICS_PARTIALS

5p Partials		3p Partials		Complete Genes		Total
100			292			23016			23408

WITHOUT_MONO_EXONICS_PARTIALS
all_genes.faa

5p Partials		3p Partials		Complete Genes		Total
63			179			23016			23258</pre>

For reference, the complete BRAKER output (`/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa`) was run through the checking script. Because the checking script requires `*` as a stop codon to determine partials and `augustus.hints.aa` does not contain stop codons, the following code was executed to re-write `augustus.hints.aa` with stop codons:
<pre style="color: silver; background: black;">
cd /UCHC/LABS/Wegrzyn/gFACs/
perl gFACs.pl -f braker_2.05_gff3 \
--statistics \
--splice-rescue \ 
--get-protein-fasta \
--fasta /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/masked_genome/genome.masked.filtered.fa \
-O /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/braker_out/ \
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gene_modeling_with_BRAKER/braker/protea/augustus.hints.gff3</pre>

The final `faa` is located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/braker_out/genes_without_introns.fasta.faa`.

This fasta was run through the checking software to yield the following statistics:
<pre style="color: silver; background: black;">
WITH_DOUBLES

Double Partials		5p Partials	3p Partials	Complete Genes	 	Total
2882			1797		3019		28196			35894

WITHOUT_DOUBLES

5p Partials	3p Partials	Complete Genes		Total
1797		3019		28196			33012

NUMBER_OF_REMOVED_GENES_WITH_GFACS

5p Partials	3p Partials	Complete Genes		Total
1734		2840		5180			9754</pre>

<h2 id="Seventh_Point_Header">Functional annotation using EnTAP</h2>

The surviving and unique gene models were then annotated using EnTAP with the following code:

<pre style="color: silver; background: black;">
module load eggnog-mapper/0.99.1
module load anaconda2/4.4.0
module load perl/5.24.0
module load diamond/0.9.19
module load python/2.7.9

/UCHC/LABS/Wegrzyn/EnTAP/EnTAP --runP -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
-d /isg/shared/databases/Diamond/RefSeq/plant.protein.faa.87.dmnd \
-i /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/gfacs_stats_and_cleaning/all_genes/all_genes.faa \
-c fungi -c bacteria --taxon Protea \
--out-dir /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/functional_annotation_with_EnTAP/entap_out \
--tcoverage 70 --qcoverage 70 \
-t 16
</pre>

The output is located at `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/annotation/functional_annotation_with_EnTAP/entap_out`.

<h2 id="uniprot">Statistics from uniprot database search</h2>
Here are the statistics for the uniprot database hits:
<pre style="color: silver; background: black;">
 Flagged contaminants (all % based on total contaminants):
                        bacteria: 68(59.65%)
                        fungi: 46(40.35%)
                Top 10 contaminants by species:
                        1)Schizosaccharomyces pombe (strain 972 / ATCC 24843): 35(30.70%)
                        2)Bacillus subtilis (strain 168): 7(6.14%)
                        3)Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv): 4(3.51%)
                        4)Escherichia coli (strain K12): 4(3.51%)
                        5)Saccharomyces cerevisiae (strain ATCC 204508 / S288c): 4(3.51%)
                        6)Synechocystis sp. (strain PCC 6803 / Kazusa): 4(3.51%)
                        7)Bacillus licheniformis: 3(2.63%)
                        8)Methylophilus methylotrophus: 2(1.75%)
                        9)Bacillus halodurans (strain ATCC BAA-125 / DSM 18197 / FERM 7344 / JCM 9153 / C-125): 2(1.75%)
                        10)Klebsiella pneumoniae: 2(1.75%)
        Top 10 alignments by species:
                        1)Arabidopsis thaliana: 5572(75.55%)
                        2)Oryza sativa subsp. japonica: 267(3.62%)
                        3)Nicotiana tabacum: 110(1.49%)
                        4)Solanum lycopersicum: 73(0.99%)
                        5)Homo sapiens: 71(0.96%)
                        6)Glycine max: 57(0.77%)
                        7)Mus musculus: 52(0.71%)
                        8)Pisum sativum: 47(0.64%)
                        9)Oryza sativa subsp. indica: 46(0.62%)
                        10)Solanum tuberosum: 39(0.53%)</pre>

<h2 id="plantfaa">Statistics from ref-seq plant faa 87 database search</h2>
And here are the statistics for the ref-seq plant protein 87 database in the taxon Protea:
<pre style="color: silver; background: black;">
        Top 10 alignments by species:
                        1)Nelumbo nucifera: 6774(56.52%)
                        2)Quercus suber: 595(4.96%)
                        3)Hevea brasiliensis: 377(3.15%)
                        4)Durio zibethinus: 342(2.85%)
                        5)Vitis vinifera: 333(2.78%)
                        6)Manihot esculenta: 301(2.51%)
                        7)Herrania umbratica: 279(2.33%)
                        8)Jatropha curcas: 218(1.82%)
                        9)Prunus avium: 204(1.70%)
                        10)Citrus clementina: 193(1.61%)
</pre>

<h2 id="combined">Integrated statistics from uniprot and ref-seq searches</h2>
Here are the statistics integrating both of the previous two searches:
<pre style="color: silver; background: black;">
    Flagged contaminants (all % based on total contaminants):
                        bacteria: 4(80.00%)
                        fungi: 1(20.00%)
                Top 10 contaminants by species:
                        1)Bacillus subtilis (strain 168): 1(20.00%)
                        2)Chlorobium phaeobacteroides (strain DSM 266): 1(20.00%)
                        3)Klebsiella pneumoniae: 1(20.00%)
                        4)Magnetospirillum magneticum (strain AMB-1 / ATCC 700264): 1(20.00%)
                        5)Schizosaccharomyces pombe (strain 972 / ATCC 24843): 1(20.00%)
        Top 10 alignments by species:
                        1)Arabidopsis thaliana: 4239(35.26%)
                        2)Nelumbo nucifera: 4161(34.61%)
                        3)Quercus suber: 261(2.17%)
                        4)Vitis vinifera: 253(2.10%)
                        5)Hevea brasiliensis: 166(1.38%)
                        6)Nicotiana tabacum: 156(1.30%)
                        7)Durio zibethinus: 147(1.22%)
                        8)Juglans regia: 131(1.09%)
                        9)Manihot esculenta: 113(0.94%)
                        10)Herrania umbratica: 107(0.89%)
</pre>
<h2 id="taxonomics">Taxonomic breakdown of EnTAP run</h2>
Lastly, here is the taxonomic breakdown of the complete run:
<pre style="color: silver; background: black;">
Total unique sequences with family assignment: 18140
Total unique sequences without family assignment: 5118
Top 10 Taxonomic Scopes Assigned:
        1)Viridiplantae: 17570(96.86%)
        2)Eukaryotes: 487(2.68%)
        3)Ancestor: 76(0.42%)
        4)Animals: 3(0.02%)
        5)Bacteria: 3(0.02%)
        6)Fungi: 1(0.01%)
</pre>

<h2 id="Eighth_Point_Header">Further statistical breakdown of EnTAP output</h2>

As there are only five contaminants, they are not included in the rest of the analysis. 

The non-contaminant alignments were extracted by combining the non-contaminant csv's first. This was done in the directory `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/functional_annotation_with_EnTAP/entap_out/`:

<pre style="color: silver; background: black;">
cat *no_contam*tsv >> no_contaminants.tsv</pre>

After this, `gFACs` was run with the `--entap-annotations` and `--annotated-all-genes-only` flags in the following script:

<pre style="color: silver; background: black;">
module load perl/5.24.0
cd /UCHC/LABS/Wegrzyn/gFACs/
perl gFACs.pl -f braker_2.05_gff3 \
--statistics \
--splice-rescue \
--entap-annotation /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/functional_annotation_with_EnTAP/entap_out/no_contaminants.tsv \
--annotated-all-genes-only \
--get-protein-fasta \
--create-gtf \
--fasta /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/masked_genome/genome.masked.filtered.fa \
-O /UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants/ \
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gene_modeling_with_BRAKER/braker/protea/augustus.hints.gff3
</pre>

The output is in `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants`.

Next, the checking software was run on `genes_without_introns.fasta.faa` in `entap_no_contaminants`. The following statistics were determined:

<pre style="color: silver; background: black;">
5p Partials	3p Partials	Complete Genes	Total
37		144		17076		17257</pre>

Here is a flow of the statistics through the annotation process:

<strong>BRAKER_OUTPUT</strong>
<pre style="color: silver; background: black;">
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gene_modeling_with_BRAKER/braker/protea/augustus.hints.aa

5p Partials	3p Partials	Complete Genes		Total
1797		3019		28196			33012
</pre>

<strong>TRIMMED_BRAKER_OUTPUT_FROM_GFACS</strong>
<pre style="color: silver; background: black;">
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/all_genes/all_genes.faa

5p Partials		3p Partials		Complete Genes		Total
63			179			23016			23258
</pre>


<strong>TRIMMED_BRAKER_OUTPUT_FROM_GFACS --> EnTAP --> ANNOTATED_GENES_ONLY_NO_CONTAMINANTS --> gFACs</strong>
<pre style="color: silver; background: black;">
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants/genes_without_introns.fasta.faa

5p Partials	3p Partials	Complete Genes	Total
37		144		17076		17257
</pre>

<strong>TOTAL_GENE_MODELS_REMOVED_THROUGH_ANNOTATION</strong>
<pre style="color: silver; background: black;">
/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants/genes_without_introns.fasta.faa
5p Partials	3p Partials	Complete Genes	Total
1760		2875		11120		15755</pre>

<strong>PERCENT_DECREASE_THROUGH_ANNOTATION</strong>
<pre style="color: silver; background: black;">
5p Partials	3p Partials	Complete Genes	Total
0.979		0.952		0.394		0.477</pre>

<h2 id="Ninth_Point_Header">Final GTF check</h2>
Before creating the `HISAT2` index, the final gtf ( `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants/out.gtf` ) was compared to the <strong>BRAKER_OUTPUT --> gFACs</strong> ( `UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/all_genes/all_genes.gtf` ) using the following code:

`bedtools intersect -s -v -a out.gtf -b all_genes.gtf >> check.gtf`

The output is located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/gtf_comparisons/`

The command above will print all features in `out.gtf` which have different coordinates in `out.gtf` and `all_genes.gtf`, with strandedness a requirement. We see from:

<pre style="color: silver; background: black;">
head check.gtf

</pre>

That there are no features in `out.gtf` which have different coordinates in `out.gtf` and `all_genes.gtf`.

Let's now see if there are any exons and introns which are overlapping. First, we need to isolate the exons and introns into separate files:

<pre style="color: silver; background: black;">sort -k1,1 -k4,4n -k5,5n -k7,7 out.gtf >> ordered.out.gtf
grep "CDS" ordered.out.gtf >> exons
grep "intron" ordered.out.gtf >> introns
sort -k1,1 -k4,4n -k5,5n -k7,7 exons >> exons.ordered
sort -k1,1 -k4,4n -k5,5n -k7,7 introns >> introns.ordered</pre>

After this, the two files were compared using the `bedtools overlap` function with a window size of `0`:

<pre style="color: silver; background: black;">windowBed -a exons.ordered -b introns.ordered -w 0 | bedtools overlap -i stdin -cols 4,5,13,14
<strong>scaffold129144	GFACS	CDS	84431	84492	0.63	+	.	g16205.t1	scaffold129144	GFACS	intron	83426	90950	1	-	.	.	61
scaffold129144	GFACS	CDS	85080	85174	0.86	+	.	g16205.t1	scaffold129144	GFACS	intron	83426	90950	1	-	.	.	94
scaffold129144	GFACS	CDS	85252	85649	0.89	+	.	g16205.t1	scaffold129144	GFACS	intron	83426	90950	1	-	.	.	397
scaffold147195	GFACS	CDS	25997	26063	0.87	+	.	g23384.t1	scaffold147195	GFACS	intron	24781	36276	1	-	.	.	66
scaffold147195	GFACS	CDS	26179	26450	0.89	+	.	g23384.t1	scaffold147195	GFACS	intron	24781	36276	1	-	.	.	271
scaffold148095	GFACS	CDS	28928	29239	0.97	-	.	g16614.t1	scaffold148095	GFACS	intron	7114	30682	1	+	.	.	311
scaffold148870	GFACS	CDS	5954	6010	0.88	+	.	g24792.t1	scaffold148870	GFACS	intron	5726	6534	0.53	+	.	.	56
scaffold148870	GFACS	CDS	6154	6285	1	+	.	g24792.t1	scaffold148870	GFACS	intron	5726	6534	0.53	+	.	.	131
scaffold169703	GFACS	CDS	63854	63992	1	-	.	g31111.t1	scaffold169703	GFACS	intron	62155	64182	1	+	.	.	138
scaffold169703	GFACS	CDS	64183	64251	0.92	+	.	g31110.t1	scaffold169703	GFACS	intron	63993	65499	1	-	.	.	68
scaffold177559	GFACS	CDS	39689	40318	0.99	-	.	g9338.t1	scaffold177559	GFACS	intron	38101	44840	1	+	.	.	629
scaffold188205	GFACS	CDS	39535	39867	0.46	-	.	g11059.t1	scaffold188205	GFACS	intron	36556	42966	1	+	.	.	332
scaffold193932	GFACS	CDS	67656	67997	0.18	-	.	g31213.t1	scaffold193932	GFACS	intron	67638	68481	0.45	-	.	.	341
scaffold34067	GFACS	CDS	69355	69623	0.77	-	.	g25544.t1	scaffold34067	GFACS	intron	61354	70850	1	+	.	.	268
scaffold34067	GFACS	CDS	69888	70035	0.79	-	.	g25544.t1	scaffold34067	GFACS	intron	61354	70850	1	+	.	.	147
scaffold36178	GFACS	CDS	9982	10167	0.51	-	.	g29143.t1	scaffold36178	GFACS	intron	10153	16942	1	-	.	.	14
scaffold5113	GFACS	CDS	11436	12026	0.97	+	.	g23104.t1	scaffold5113	GFACS	intron	7440	13074	1	+	.	.	590
</strong></pre>

Because strandedness was not considered, we grab columns which are on the same strand:
<pre style="color: silver; background: black;"><strong>scaffold148870	GFACS	CDS	5954	6010	0.88	+	.	g24792.t1	scaffold148870	GFACS	intron	5726	6534	0.53	+	.	.	56
scaffold148870	GFACS	CDS	6154	6285	1	+	.	g24792.t1	scaffold148870	GFACS	intron	5726	6534	0.53	+	.	.	131
scaffold193932	GFACS	CDS	67656	67997	0.18	-	.	g31213.t1	scaffold193932	GFACS	intron	67638	68481	0.45	-	.	.	341
scaffold36178	GFACS	CDS	9982	10167	0.51	-	.	g29143.t1	scaffold36178	GFACS	intron	10153	16942	1	-	.	.	14
scaffold5113	GFACS	CDS	11436	12026	0.97	+	.	g23104.t1	scaffold5113	GFACS	intron	7440	13074	1	+	.	.	590
</strong></pre>

Let's see if any of these overlaps occur in nested models. To determine if there are any nested models we first use `bedtools merge` to collapse our gtf into the gene models:

<pre style="color: silver; background: black;">bedtools merge -s -i out.ordered.gtf -c 2,7,9 -o distinct,distinct,distinct >> merged.out.gtf
head merged.out.gtf
<strong>C67174424	61	2723	GFACS	-	.,g32988,g32988.t1
C67174986	967	1439	GFACS	-	.,g29054,g29054.t1
C67175022	263	2282	GFACS	-	.,g27556,g27556.t1
C67175416	710	1988	GFACS	+	.,g10622,g10622.t1
C67175522	359	1961	GFACS	+	.,g11052,g11052.t1
C67177288	820	3024	GFACS	+	.,g33243,g33243.t1
C67177542	1088	1580	GFACS	-	.,g19881,g19881.t1
C67178404	542	2757	GFACS	-	.,g23770,g23770.t1
C67178796	492	2142	GFACS	-	.,g1325,g1325.t1
C67179430	56	2864	GFACS	+	.,g24443,g24443.t1
</strong>

wc -l merged.out.gtf
<strong>17251</strong>

grep -c "gene" merged.out.gtf
<strong>17257</strong></pre>

We already that six genes have been lost in the merge. To identify these lost genes the last column was split into two columns with:

<pre style="color: silver; background: black;">head merged.out.gtf
<strong>C67174424	61	2723	GFACS	-	.,g32988,g32988.t1
C67174986	967	1439	GFACS	-	.,g29054,g29054.t1
C67175022	263	2282	GFACS	-	.,g27556,g27556.t1
C67175416	710	1988	GFACS	+	.,g10622,g10622.t1
C67175522	359	1961	GFACS	+	.,g11052,g11052.t1
C67177288	820	3024	GFACS	+	.,g33243,g33243.t1
C67177542	1088	1580	GFACS	-	.,g19881,g19881.t1
C67178404	542	2757	GFACS	-	.,g23770,g23770.t1
C67178796	492	2142	GFACS	-	.,g1325,g1325.t1
C67179430	56	2864	GFACS	+	.,g24443,g24443.t1</strong>

sed -i 's/\.\,//g' merged.out.gtf
sed -i 's/\.t1//g' merged.out.gtf
sed -i 's/\,/        /g' merged.out.gtf
column -t merged.out.gtf >> formatted.merged.out.gtf
head formatted.merged.out.gtf 
<strong>C67174424                      61      2723    GFACS  -  g32988  g32988
C67174986                      967     1439    GFACS  -  g29054  g29054
C67175022                      263     2282    GFACS  -  g27556  g27556
C67175416                      710     1988    GFACS  +  g10622  g10622
C67175522                      359     1961    GFACS  +  g11052  g11052
C67177288                      820     3024    GFACS  +  g33243  g33243
C67177542                      1088    1580    GFACS  -  g19881  g19881
C67178404                      542     2757    GFACS  -  g23770  g23770
C67178796                      492     2142    GFACS  -  g1325   g1325
C67179430                      56      2864    GFACS  +  g24443  g24443</strong>

awk -F"\t" '$6 != $7 { print $0 }' formatted.merged.out.gtf
</pre>

No genes were collapsed together based on this. But to make sure, the following code was executed:
<pre style="color: silver; background: black;">
bedtools merge -s -i out.ordered.gtf -c 3,7,9 -o collapse,distinct,distinct >> merged.out.gtf.copy
sed -i 's/CDS//g' *copy
sed -i 's/intron//g' *copy
sed -i 's/start_codon//g' *copy
sed -i 's/stop_codon//g' *copy
head *copy
<strong>C67174424	61	2723	,,gene,,,,,,,,,,,	-	.,g32988,g32988.t1
C67174986	967	1439	,,gene,,,	-	.,g29054,g29054.t1
C67175022	263	2282	,,gene,,,,,	-	.,g27556,g27556.t1
C67175416	710	1988	,,gene,	+	.,g10622,g10622.t1
C67175522	359	1961	,,gene,	+	.,g11052,g11052.t1
C67177288	820	3024	,,gene,,	+	.,g33243,g33243.t1
C67177542	1088	1580	,,gene,,,,,	-	.,g19881,g19881.t1
C67178404	542	2757	,,gene,,,,,	-	.,g23770,g23770.t1
C67178796	492	2142	,,gene,,,,,	-	.,g1325,g1325.t1
C67179430	56	2864	,,gene,,,	+	.,g24443,g24443.t1
</strong>
sed -i 's/\,//g' *copy
grep "genegene" *copy
<strong>scaffold121444	10355	11838	genegene	-	.g14240g14240.t1g14241g14241.t1
scaffold148870	714	6690	genegene	+	.g24790g24790.t1g24792g24792.t1
scaffold193932	67294	70638	genegene	-	.g31212g31212.t1g31213g31213.t1
scaffold36178	4431	23074	genegene	-	.g29143g29143.t1g29144g29144.t1
scaffold44960	9707	18008	genegene	+	.g31619g31619.t1g31620g31620.t1
scaffold5113	6708	31552	genegene	+	.g23103g23103.t1g23104g23104.t1</strong></pre>

This is why we always second check! We see that there are some genes which were combined, and the result is _four_ gene names combined together. Because of this, when we split the table earlier the result for the first row would have been 

<pre style="color: silver; background: black;">
scaffold121444	10355	11838	genegene	-	g14240	g14240	g14241	g14241</pre>

This row would have 9 columns, and when we check if columns 6 and 7 are equal these particular rows would pass. However, we see that genes 14241 and 14240 are actually one gene nested inside of another. Therefore, the genes which are nested together are:

<pre style="color: silver; background: black;">g14240
g14241
g24790
g24792
g29143
g29144
g31619
g31620
g23103
g23104</pre>

We want to remove these genes from both the fasta and the gtf. Let's remove them from the gtf first:

<pre style="color: silver; background: black;">grep -vwE "(g14240|g14241|g24790|g24792|g29143|g29144|g31619|g31620|g23103|g23104)" out.ordered.gtf >> nested_genes_removed.gtf</pre>

And re-do our steps:
<pre style="color: silver; background: black;">bedtools merge -s -i nested_genes_removed.gtf -c 3,7,9 -o distinct,distinct,distinct >> merged_nested_genes_removed.gtf

wc -l merged_nested_genes_removed.gtf 
<strong>17308 merged_nested_genes_removed.gtf</strong>
grep -c "gene" nested_genes_removed.gtf 
<strong>17247</strong>
grep -c "gene" merged_nested_genes_removed.gtf 
<strong>17246</strong></pre>

Any merged features which did not represent a gene were removed:

<pre style="color: silver; background: black;">grep "gene" merged_nested_genes_removed.gtf >> merged_nested_genes_removed.genes
wc -l merged*genes
17246 merged_nested_genes_removed.genes</strong></pre>

The numbers check out. Our final gtf is `nested_genes_removed.gtf`.

The last step in our `GTF` check is to make sure that there are no merged modules which are multiple introns, exons, or both concatenated (that is, to make sure we aren't actually missing what could be a gene):
<pre style="color: silver; background: black;">
grep -v "gene" merged_nested_genes_removed.gtf >> non_gene_modules
head non_gene_modules
<strong>scaffold121444	10355	10358	stop_codon	-	.
scaffold121444	10672	10675	stop_codon	-	.
scaffold121444	10775	10778	start_codon	-	.
scaffold121444	11167	11679	intron	-	.
scaffold121444	11835	11838	start_codon	-	.
scaffold148870	714	717	start_codon	+	.
scaffold148870	777	5154	intron	+	.
scaffold148870	5327	5696	intron	+	.
scaffold148870	5725	6534	intron,start_codon	+	.
scaffold148870	6623	6626	stop_codon	+	.</strong><

grep "intron" non_gene_modules >> non_gene_introns
wc -l non_gene_introns
46 non_gene_introns
nano non_gene_introns

<strong>OUTPUT OMITTED, NO POSSIBLE GENE MODULES IN OUTPUT</strong>

grep "CDS" non_gene_modules >> non_gene_exons
wc -l non_gene_exons
<strong>0 non_gene_exons</strong></pre>

Lastly, we need to remove the genes from the fasta. We place our nested gene ids into the file `nested_list`. Next, we use `awk` to retrieve the sequences, and `grep` to check that the actual protein sequences do not appear in any other genes:

<pre style="color: silver; background: black;">id='(cat nested_list)'
for gene in $id; 
do awk '/'$gene'/{flag=1;print $0;next}/^>/{flag=0}flag' genes_without_introns.fasta.faa >> no_nested_gene_models.fasta ;
done;
cat no_nested_gene_models.fasta | xargs -Ivar grep var genes_without_introns.fasta.faa >> no_nested_gene_models.check;
diff no_nested_gene_models.fasta no_nested_gene_models.check</pre>

After determining the protein sequences were unique, they were removed from the fasta and checked to ensure that the remaining sequences are undisturbed:

<pre style="color: silver; background: black;">
grep -vFf no_nested_gene_models.fasta genes_without_introns.fasta.faa >> genes_without_introns_or_nests.fasta.faa
awk 'NR==FNR{array[$1];next}!($1 in array){print $1}' no_nested_gene_models.fasta genes_without_introns_or_nests.fasta.faa >> genes_without_introns_or_nests.check

sort genes_without_introns_or_nests.check >> sorted_genes_without_introns_or_nests.check
sort genes_without_introns_or_nests.fasta.faa >> sorted_genes_without_introns_or_nests.fasta.faa
diff sorted_genes_without_introns_or_nests.check sorted_genes_without_introns_or_nests.fasta.faa</pre>

The fasta passes the check. The final fasta file is located at: `/UCHC/LABS/Wegrzyn/proteaBraker/braker/protea/wolfo_analysis/gfacs_stats_and_cleaning/entap_no_contaminants/genes_without_introns_or_nests.fasta.faa`.

Lastly:
<pre style="color: silver; background: black;">
grep -c ">" genes_without_introns_or_nests.fasta.faa
<strong>17247</strong>
grep -c "gene" nested_genes_removed.gtf
<strong>17247</strong></pre>



