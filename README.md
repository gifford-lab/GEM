<h2>GEM</h2><br>
 High resolution peak calling and motif discovery for ChIP-seq and ChIP-exo data</h1> 
      <font size='+2'><b>G</b>enome wide <font size='+2'><b>E</b>vent finding and <font size='+2'><b>M</b>otif discovery
      <p><font size="-1"><b>Citation</b>:<br>
      <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002638">High Resolution Genome Wide Binding Event Finding and Motif Discovery Reveals Transcription Factor Spatial Binding Constraints</a>. Yuchun Guo, Shaun Mahony &amp; David K Gifford, (2012) PLoS Computational Biology 8(8): e1002638.

<p><font size="-1"><b>Abstract</b>:<br>

An essential component of genome function is the syntax of genomic regulatory elements that determine how diverse transcription factors interact to orchestrate a program of regulatory control. A precise characterization of in vivo spacing constraints between key transcription factors would reveal key aspects of this genomic regulatory language. To discover novel transcription factor spatial binding constraints in vivo, we developed a new integrative computational method, genome wide event finding and motif discovery (GEM). GEM resolves ChIP data into explanatory motifs and binding events at high spatial resolution by linking binding event discovery and motif discovery with positional priors in the context of a generative probabilistic model of ChIP data and genome sequence. GEM analysis of 63 transcription factors in 214 ENCODE human ChIP-Seq experiments recovers more known factor motifs than other contemporary methods, and discovers six new motifs for factors with unknown binding specificity. GEM's adaptive learning of binding-event read distributions allows it to further improve upon previous methods for processing ChIP-Seq and ChIP-exo data to yield unsurpassed spatial resolution and discovery of closely spaced binding events of the same factor. In a systematic analysis of in vivo sequence-specific transcription factor binding using GEM, we have found hundreds of spatial binding constraints between factors. GEM found 37 examples of factor binding constraints in mouse ES cells, including strong distance-specific constraints between Klf4 and other key regulatory factors. In human ENCODE data, GEM found 390 examples of spatially constrained pair-wise binding, including such novel pairs as c-Fos:c-Jun/USF1, CTCF/Egr1, and HNF4A/FOXA1. The discovery of new factor-factor spatial constraints in ChIP data is significant because it proposes testable models for regulatory factor interactions that will help elucidate genome function and the implementation of combinatorial control.      
<p><font size="-1"><b>News</b>:<br>
<li><a href="https://scholar.google.com/scholar?oi=bibs&hl=en&cites=8022792170454236947&as_sdt=5">Papers citing GEM</a>
<li>GEM has been selected to be part of the <a href="https://www.encodeproject.org/software/gem/">ENCODE TF ChIP-seq analysis pipeline.</a>
<li> MIT NEWS <a href="http://web.mit.edu/newsoffice/2012/deciphering-the-language-of-transcription-factors-0910.html"><b>Deciphering the language of transcription factors</b></a> (MIT News article on the GEM paper).
</ul>
      <p><font style="font-family: arial,sans-serif; font-size: 10pt;"><b>GEM is a scientific software for studying protein-DNA interaction at high resolution using ChIP-seq/ChIP-exo data. It can also be applied to CLIP-seq and Branch-seq data.</b> <br>GEM links binding event discovery and motif discovery with positional priors in the context of a generative probabilistic model of ChIP data and genome sequence, resolves ChIP data into explanatory motifs and binding events at unsurpassed spatial resolution. GEM reciprocally improves motif discovery using binding event locations, and binding event predictions using discovered motifs.  </p>
      <p><font style="font-family: arial,sans-serif; font-size: 10pt;">GEM has following features:
      </p>
      <ol> 
	<li><font style="font-family: arial,sans-serif; font-size: 10pt;">Exceptionally high spatial resolution on binding event prediction (aka peak calling)</li>
	<li>Highly accurate <i>de novo</i> motif discovery</li>
	<li>Resolve closely spaced (less than 500bp) homotypic events that appear as a single cluster of reads</li>
	<li>Enable analysis of spatial binding constraints of multiple transcription factors, for predicting TF dimer binding sites, enhanceosomes, etc. </li>
	<li>Analyze ChIP-seq, ChIP-exo, CLIP-seq and Branch-seq data, single-end or paired-end</li>
        <li>Run in single-condition mode  or multi-condition mode</li>
</ol>
GEM is implemented in Java, which comes with all the major operating systems.<p>
<font color="red"><a href='http://cgs.csail.mit.edu/gem/ES_Sox2_example/ES_Sox2_2_result.htm'>See an example of GEM output for ES cell Sox2 binding</a></font>, including binding events, K-mer Set motifs (KSMs), PWM motifs, and motif spatial distribution plots.


<p></p>


<h3 style="color: RoyalBlue;"><a class="anchor" id="download"></a>Download</h3>
<p>Download, unzip, and run ... see command line <a href="#examples">examples</a>.<p>
<a href="http://cgs.csail.mit.edu/gem/versions.html"> Download GEM software (version 2.6) and test data </a>

<h3 style="color: RoyalBlue;"><a class="anchor" id="gemvsgps"></a>GEM vs. GPS</h3>
<p>GEM is a superset of <a href="http://cgs.csail.mit.edu/gps/">GPS</a>. GEM uses both ChIP-seq read data and genome sequence to perform de novo motif discovery and binding event calling, while GPS uses only  ChIP-seq read data for binding event calling. <p>
GEM can be activated by giving a genome sequence (<code>--genome</code>) and using any one of the following command line options:
<ul>
<li><code>--k</code>: the length of the k-mers
<li><code>--k_min</code> and <code>--k_max</code>: the range for the length of k-mers
<li><code>--seed</code>: the seed k-mer to jump start k-mer set motif discovery. The length of the seed k-mer will be used to set k.
</ul>
If these three options are not used, GEM will just run GPS and stop.
<p>

<h3 style="color: RoyalBlue;"><a class="anchor" id="requirements"></a>System requirements</h3>

<p>GEM is a java software. Not installation is required. It works across most of computer systems. 
<p>Java 7 is required to execute the JAR. For analysis with mamalian genome, GEM requires about 5-15G memory. It can be specified at the command line with the option <code>-Xmx</code> (i.e. <code>java
-jar -Xmx10G gem.jar</code> allocates 10GB of memory). </p>


<h3 style="color: RoyalBlue;"><a class="anchor" id="readdistrib"></a>Read distributions</h3>


<p>As GPS, 
read distribution file is required for GEM.  The
user can use the default read distribution file provided with the
software as starting point (<a href="http://cgs.csail.mit.edu/gem/download/Read_Distribution_default.txt">ChIP-seq default read distribution file</a> <a href="download/Read_Distribution_ChIP-exo.txt">ChIP-exo default read distribution file</a>). After one round of prediction, GEM will 
re-estimate the read distribution using the predicted events.<br>
<br>
The
read distribution file specifies
the empirical spatial distribution of reads for a given binding event.
The file contains tab-delimited position/density pairs. The first field
is the position relative to the binding position (i.e. binding event is
at position 0). The second field contains the corresponding read
density at that position. For example,</p>
<p>-344&nbsp;&nbsp;&nbsp;
1.42285E-4<br>
-343&nbsp;&nbsp;&nbsp; 1.42275E-4<br>
-342&nbsp;&nbsp;&nbsp; 1.42265E-4<br>
... ... </p>
<p>Alternatively,
it can be estimated directly from the ChIP-seq data. Given a set of
events, we count all the reads at each position (the 5' end of the
reads) relative to the corresponding event positions. The initial set
of events for estimating the empirical spatial distribution can be
defined by using known motifs or by finding the center of the forward
and reverse read profiles. GEM has a tool to calculate
the read distribution from a user provided file (coords.txt) containing
the coordinates (format chr:coord, e.g. 1:3498832, or chr:coord:strand, e.g. 1:3498832:+, each coordinate on a separate line).</p>
<p><code>java -Xmx1G -cp gem.jar edu.mit.csail.cgs.deepseq.analysis.GPS_ReadDistribution --g "mm8.chrom.sizes" --coords "coords.txt" --chipseq "SRX000540_mES_CTCF.bed" --f "BED" --name "CTCF" --range 250 --smooth 5 --mrc 4</code>
</p>
<font color='red' style="font-family: arial,sans-serif; font-size: 10pt;">For ChIP-exo data, we suggest  <code>--smooth 3 --mrc 20</code>.</font>
<p>

After GPS/GEM makes the prediction, it will re-estimate the read
distribution using the predicted events. 
<p> If
the data are too noisy or too few events are used for re-estimation,
the new read distribution may not be accurate. The users are encouraged
to examine the read distribution using the plot of read distributions
(X_All_Read_Distributions.png) output by GEM. </p>


<h3 style="color: RoyalBlue;"><a class="anchor" id="io"></a>Input and output</h3>


<p>
GEM
takes an alignment file of ChIP-seq reads and a genome sequence as input and reports a list
of predicted binding events and the explanatory binding motifs. </p>
<p><b>ChIP-seq alignment file</b><br>

ChIP-seq alignment file formats that are supported: </p>
<ul>
<li><a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">BED</a></li>
<li><a href="http://samtools.sourceforge.net/SAM1.pdf">SAM</a></li>
<!--li><a href="http://bowtie-bio.sourceforge.net/index.shtml">Bowtie</a></li>
<li><a href="http://www.illumina.com/documents/products/datasheets/datasheet_genomic_sequence.pdf">ELAND</a></li>
<li><a href="http://www.novocraft.com/products.html">NovoAlign</a></li-->
</ul>
Support for some file formats may not be the most updated. If you have an error related to FileReadLoader, please do the following:
<ol>
<li>Send me a few example lines (and the version of the aligner software), I will update it in the next release.
<li>A quick fix is to convert the read alignment file into BED format, and then try again.
</ol>

<b>Genome sequence</b><br>

To run de novo motif discovery, a genome sequence (<a href='http://hgdownload.cse.ucsc.edu/downloads.html'>UCSC download</a>) is needed. The path to directory containing the genome sequence files (by chromosome, *.fa or *.fasta files, with the prefix "chr") can be specified using option <code>--genome</code> (for example, <code>--genome your_path/mm8/</code>). Note that the chromosome name should match those in the "<code>--g</code>" genome.chrom.sizes file, as well as those in your read alignment file.<p>
For example, the fasta file for chromosome 2 is chr2.fa:
<code><br>
>chr2<br>
TAATTGTAATAGTATATACTTGTATGTACTTAAAATAttttatcatagtt<br>
ATCTGGATTTTTGATGGCTATCATGACCTCTGAATGACTAGGGAATCTTG<br>
... ...
</code>

<p>
<b>Output files</b><br>
GEM outputs both the binding event files and the motif files. 
Because of the read distribution re-estimation, GEM outputs event
prediction and read distribution files for multiple rounds. <a href="#QA">(See more details)</a></p>
<ol>
<li>GEM event text file (GEM_events.txt, <a href="#gem">see more details</a>) 
<li>K-mer set motifs (KSM.txt, <a href="#ksm">see more details</a>)
<li>PFM file of PWM motifs (PFM.txt)
<li>HTML file summarizing the GEM event and motif results <a href='http://cgs.csail.mit.edu/gem/ES_Sox2_example/ES_Sox2_2_result.htm'>(see an example and explanations)</a>
<li>GEM output folder containing more detailed result files (Round 1 and 2 are GPS and GEM results, respectively <a href="#QA">(See more details)</a>)
<ul>
    <li>GEM event text files (significant, insignificant and filtered)
    <li>K-mer set motifs
    <li>PFM file of PWM motifs
    <li>Read distribution file
    <li>Spatial distribution between primary motif and all the secondary motifs in the 61bp around the GEM events (this is based on PWM motif match, not based on the GEM binding sites). You will have the files such as Name_Spatial_dist_0_1.txt/png showing the spatial distribution of secondary motifs (#1) with respect to the primary motif (#0) in text format or as a PNG image. If you click on the PNG image on the HTML output page, you also get the txt file with the values.
    </ul>
</ol>


<p>
GEM also outputs the list of insignificant events (those do not pass
the statistical test), and the filtered events (those would pass the
statistical test using the read count, but have a low IP/Ctrl ratio, or
the distribution of reads are quite different from the empirical
distribution).


<p>Optionally, GEM can be set to output BED files (using option <code>--outBED</code>)
for loading the GEM results to Genome Browser as custom tracks for
visualization. Note that in BED file, the coordinates are offset to
left/right 100bp to give a region for visualization. </p>

<p><a class="anchor" id="gem"></a>GEM event file is a tab-delimited file (xxx_n_GEM_events.txt ) with following fields: </p>

<table style="font-family: arial,sans-serif; font-size: 10pt;" border="1">
<tr><th>Field</th><th>Description</th></tr>
<tr><td>Location</td><td>the genome coordinate of this binding event</td></tr>
<tr><td>IP binding strength</td><td>the number of IP reads associated with the event
</td></tr><tr><td>Control binding strength</td><td>the number of control reads in the corresponding region
</td></tr><tr><td>Fold</td><td>fold enrichment (IP/Control) 
</td></tr><tr><td>Expected binding strength</td><td>the number of IP read counts expected in the binding region given its local context (defined by parameter W2 or W3), this is used as the Lambda parameter for the Poisson test
</td></tr><tr><td>Q_-lg10</td><td>-log10(q-value), the q-value after multiple-testing correction, using the larger p-value of Binomial test and Poisson test
</td></tr><tr><td>P_-lg10</td><td>-log10(p-value), the p-value is computed from the Binomial test given the IP and Control read counts (when there are control data)
</td></tr><tr><td>P_poiss</td><td>-log10(p-value), the p-value is computed from the Poission test given the IP and Expected read counts (without considering control data)
</td></tr><tr><td>IPvsEMP</td><td>Shape deviation, the KL divergence of the IP reads from the empirical read distribution (log10(KL)), this is used to filter predicted events given the <code>--sd</code> cutoff (default=-0.40).
</td></tr><tr><td>Noise</td><td>the fraction of the event read count estimated to be noise               
</td></tr><tr><td>KmerGroup</td><td>the group of the k-mers associated with this binding event, only the most significant k-mer is shown, the n/n values are the total number of sequence hits of the k-mer group in the positive and negative training sequences (by default total 5000 of each), respectively
</td></tr><tr><td>KG_hgp</td><td>log10(hypergeometric p-value), the significance of enrichment of this k-mer group in the positive vs negative training sequences (by default total 5000 of each), it is the hypergeometric p-value computed using the pos/neg hit counts and total counts
</td></tr><tr><td>Strand</td><td>the sequence strand that contains the k-mer group match, the orientation of the motif is determined during the GEM motif discovery, '*' represents that no k-mer is found to associated with this event
</td></tr></table>

<p>
<p><a class="anchor" id="ksm"></a>The KSM file is a tab-delimited file (xxx_n_KSM.txt ) with following fields: </p>
First header line, e.g. "#5000/5000", shows the number of positive/negative sequences that were used for learing the motif.<br>
Second header line, e.g. "#3.01", the KSM motif score cutoff, optimized to give best motif enrichment in the training sequences.
<table style="font-family: arial,sans-serif; font-size: 10pt;" border="1">
<tr><th>Field</th><th>Description</th></tr>
<tr><td>k-mer/r.c.</td><td>The k-mer sequence and its reverse complement</td></tr>
<tr><td>Cluster</td><td>Cluster ID of the k-mer set (Cluster 0 is the primary motif, i.e. the most significant motif)</td></tr>
<tr><td>Offset</td><td>The offset of this k-mer from the seed k-mer</td></tr>
<tr><td>PosCt</td><td>Number of positive sequences that contain this k-mer</td></tr>
<tr><td>wPosCt</td><td>Weighted PosCt, calculated using the relative sequence weighting (for GEM, this is natural logarithm of the binding event strength, i.e. read count)</td></tr>
<tr><td>NegCt</td><td>Number of negative sequences that contain this k-mer</td></tr>
<tr><td>HGP_10</td><td>HyperGeometric P-value (log10)</td></tr>
<tr><td>(no name)</td><td>The IDs of positive sequences that contain this k-mer</td></tr>
<tr><td>(no name)</td><td>The IDs of negative sequences that contain this k-mer</td></tr>
</table>


<h3 style="color: RoyalBlue;"><a class="anchor" id="examples"></a>Examples:</h3>

<a href="download/gps_test.tar.gz">This data</a>
can be used to test GEM. It comes from a Ng lab publication <a href="http://www.ncbi.nlm.nih.gov/pubmed/18555785">(PMID:
18555785)</a> and consists of Bowtie alignments of mouse ES cell
CTCF ChIP-seq and GFP control reads. 
<p>Once everything is unpacked, use the following command: <br>
<code>java -Xmx10G -jar gem.jar --d Read_Distribution_default.txt --g mm8.chrom.sizes --genome your_path/mm8 --s 2000000000 --expt SRX000540_mES_CTCF.bed --ctrl SRX000543_mES_GFP.bed --f BED --out mouseCTCF --k_min 6 --k_max 13</code> </p>
<p>Note the double dashes (<code>--</code>) for GEM parameters.
<p><font color='red'>For ChIP-exo data, use ChIP-exo read distribution <code>--d Read_Distribution_ChIP-exo.txt</code>, add one more option <code>--smooth 3</code> to estimate the read distribution without too much smoothing. Depending on the quality of the data, you may want to turn off the read filtering by option <code>--nrf</code></font>.

<!--p> An
example of GEM run in multi-condition alignment mode is (Please note:
the multi-condition mode may take longer time to run) <br>
<code> java -Xmx5G -jar gem.jar --d Read_Distribution_default.txt --s 240000000 --exptCTCF-GM12878 GM12878_chr1_ip.bed --exptCTCF-HUVEC HUVEC_chr1_ip.bed --ctrlCTCF-GM12878 GM12878_chr1_ctrl.bed --ctrlCTCF-HUVEC HUVEC_chr1_ctrl.bed --f BED --g hg18_chr1.chrom.sizes --out HumanCTCF </code><br>
Note the matching <code>--expt</code> and <code>--ctrl</code> options that encode the name of different conditions. </p-->

<h3 style="color: RoyalBlue;"><a class="anchor" id="options"></a>Command-line options:</h3>

<p> The
command line parameters are in the format of <code>--flag/name</code>
pairs. Note the double dashes (<code>--</code>); GEM does not accept single dash parameters.
<p>Some parameters are required: </p>

<table style="font-family: arial,sans-serif; font-size: 10pt;" border="1">
<tbody>
<tr>
<th>Required parameters</th>
<th>Detailed information</th>
</tr>
<tr>
<td><code>--d [path]</code></td>
<td>The path to the read distribution model file</td>
</tr>
<tr>
<td><code>--exptX [path]</code></td>
<td>The path to the aligned reads file for experiment
(IP). X is condition name. In multi-condition alignment mode, X is
used to specify different conditions.</td>
</tr>
</tbody>
</table>

<p>Some parameters (those relevant to motif discovery are in bold) are optional: </p>

<table style="font-family: arial,sans-serif; font-size: 10pt;" border="1">
<tbody>
<tr>
<th>Optional parameters</th>
<th>Detailed information</th>
</tr>
<tr>
<td><code>--ctrlX [path]</code></td>
<td>The path to the aligned reads file for control. X
should match the condition name in <code>--exptX</code>.</td>
</tr>
<tr>
<td><code>--g [path]</code></td>
<td>The path to a <a href='versions.html'>genome information file (genome.chrom.sizes file)</a>. The file contains tab-delimited chromosome name/length pairs. Highly recommended, although not required. If it is not supplied, GEM will use the maximum value of read coordinate as the chomosome length.</td>
</tr>
<tr>
<td><code>--f [BED|SAM|BOWTIE|ELAND|NOVO]</code></td>
<td>Read file format: BED/SAM/BOWTIE/ELAND/NOVO. The
SAM option allows SAM or BAM file input. (default = BED)</td>
</tr>
<tr>
<td><code>--s [n]</code></td>
<td>The size of uniquely mappable genome. It depends on the genome and the read length. A good estimate is genome size * 0.8. If it is not supplied, it will be estimated using the genome information.</td>
</tr>
<tr>
<td><b><code>--genome [path]</code></b></td>
<td>the path to the genome sequence directory, which contains fasta files by chromosomes</td>
</tr>
<tr>
<td><b><code>--k [n]</code></b></td>
<td>the width of the k-mers</td>
</tr>
<tr>
<td><b><code>--k_min [n] --k_max [n]</code></b></td>
<td>minimum and maximum value of k</td>
</tr>
<tr>
<td><b><code>--seed [k-mer]</code></b></td>
<td>the seed k-mer to jump start k-mer set discovery. Exact k-mer sequence only. The width of the seed k-mer will be used to set k</td>
</tr>
<tr>
<td><b><code>--k_seqs [n]</code></b></td>
<td>the number of top ranking events to get sequences for motif discovery (default=5000)</td>
</tr>
<tr>
<td><b><code>--pp_nmotifs [n]</code></b></td>
<td>the max number of top ranking motifs to set the motif-based positional prior (default=1)</td>
</tr>

<tr>
<td><b><code>--k_win [n]</code></b></td>
<td>the sequence window size around the binding event for motif discovery (default=61bp)</td>
</tr>
<tr>
<td><b><code>--strand_type [n]</code></b></td>
<td>Double-strand or single-strand binding event calling and motif discovery, 0 for double-strand, 1 for single-strand (default=0)</td>
</tr>
<tr>
<td><code>--nd [n]</code></td>
<td>noise distribution model, 0 for no noise model, 1 for uniform noise model (default=1)</td>
</tr>
<!--tr>
<td><code>--mrc [value]</code></td>
<td>Maximum read count on a base postion.
(default=-1, GPS will estimate from the data based on a Possion model,
with a minimum value of 3)</td>
</tr-->
<tr>
<td><code>--fold [value]</code></td>
<td>Fold (IP/Control) cutoff to filter predicted events (default=3)</td>
</tr>
<tr>
<td><code>--icr [value]</code></td>
<td>IP/Control Ratio. By default, this ratio is estimated from the data using non-specific binding regions. It is important to set this value explicitly for synthetic dataset.</td>
</tr>
<tr>
<td><code>--out [name]</code></td>
<td>Output folder and file name prefix</td>
</tr>
<tr>
<td><code>--q [value]</code></td>
<td>significance level for q-value, specified as -log10(q-value). For example, to enforce a q-value threshold of 0.001,
set this value to 3. (default=2, i.e. q-value=0.01)</td>
</tr>
<tr>
<td><code>--a [value]</code></td>
<td>minimum alpha value for sparse prior (default is estimated by mean whole genome read coverage)</td>
</tr>
<tr>
<td><code>--af [value]</code></td>
<td>a constant to scale alpha value with read count
(default=3). A smaller af value will give a larger alpha value.</td>
</tr>
<tr>
<td><code>--sd [value]</code></td>
<td>Shape deviation cutoff to filter predicted events (default=-0.40).</td>
</tr>
<tr>
<td><code>--smooth [n]</code></td>
<td>The width (bp) to smooth the read distribution. If it is set to -1, there will be no smoothing (default=30).</td>
</tr>
<tr>
<td><code>--subs [region strings]</code></td>
<td>Subset of genome regions to be analyzed. The string can be in the format of "chr:start-end" or "chr", or both. For example, "1:003-1004 2 X".</td>
</tr>
<tr>
<td><code>--subf [region file]</code></td>
<td>File that contains subset of genome regions to be analyzed. Each line contains a region in "chr:start-end" format.</td>
</tr>
<tr>
<td><code>--ex [region file]</code></td>
<td>File that contains subset of genome regions to be excluded. Each line contains a region in "chr:start-end" format.</td>
</tr>
<tr>
<td><code>--t [n]</code></td>
<td>Number of threads to run GEM in paralell. It is suggested to be equal to or less than the physical CPU number on the computer. (default: physical CPU number)</td>
</tr>
<tr>
<td><code>--top [n]</code></td>
<td>Number of top ranked GEM events to be used for re-estimating the read distribution (default=2000). Note that GEM only re-estimate
when there are more than 500 significant events called.</td>
</tr>
<tr>
<td><code>--w2 [n]</code></td>
<td>Size of sliding window to estimate lambda
parameter for Possion distribution when there is no control data
(default=5,000, must be larger than 1000).</td>
</tr>
<tr>
<td><code>--w3 [n]</code></td>
<td>Size of sliding window to esitmate lambda
parameter for Possion distribution when there is no control data
(default=10,000, must be larger than w2).</td>
</tr>
</tbody>
</table>

<p>Optional
flags: </p>

<table style="font-family: arial,sans-serif; font-size: 10pt;" border="1">
<tbody>
<tr>
<th>Optional flags</th>
<th>Detailed information</th>
</tr>
<tr>
<td><code>--k_neg_dinu_shuffle</code></td>
<td>Use di-nucleotide shuffled sequences as negative sequences for motif finding</td>
</tr>
<tr>
<td><code>--bp</code></td>
<td>use Branch-seq data specific settings</td>
</tr>
<tr>
<td><code>--pp_pwm</code></td>
<td>Use PWM motif to set the motif-based positional prior (default is to use the KSM motif model)</td>
</tr>
<tr>
<td><code>--bf</code></td>
<td>Depreciated<br>
Base filtering is done by Poisson filtering by default.</td>
</tr>
<tr>
<td><code>--fa</code></td>
<td>GEM will use a fixed user-specified alpha value for all the regions</td>
</tr>
<tr>
<td><code>--help</code></td>
<td>Print help information and exit</td>
</tr>
<tr>
<td><code>--multi</code></td>
<td>Depreciated
<br>To run GEM in multi-condition mode, you only need to specify data for
conditions X and Y using --exptX and --exptY, etc.</td>
</tr>
<tr>
<td><code>--nf</code></td>
<td>Do not filter predicted events (default is to perform event filtering by shape and fold)</td>
</tr>
<tr>
<td><code>--nrf</code></td>
<td>Do not filter duplicate reads (i.e. potential PCR duplicates) (default is to apply filtering considering the read counts at its neighboring bases)
</td>
</tr>
<tr>
<td><code>--outNP</code></td>
<td>Output binding events in ENCODE NarrowPeak format (default is NO narrowPeak file output)</td>
</tr>
<tr>
<td><code>--outBED</code></td>
<td>Output binding events in BED format for UCSC Genome Browser (default is NO BED file output)</td>
</tr>
<tr>
<td><code>--outJASPAR</code></td>
<td>Output motif PFM in JASPAR format</td>
</tr>
<tr>
<td><code>--outMEME</code></td>
<td>Output motif PFM in  MEME format</td>
</tr>
<tr>
<td><code>--outHOMER</code></td>
<td>Output motif PFM in HOMER format</td>
</tr>
<tr>
<td><code>--sl</code></td>
<td>Sort GEM output by location (default is sorted by P-value)</td>
</tr>
</tbody>
</table>

<h3 style="color: RoyalBlue;"><a class="anchor" id="QA"></a>Q&amp;A</h3>
<span style="font-weight: bold;">Which round of result should I use?</span><br>
Because of the read distribution re-estimation, GEM may output event prediction and read distribution files for multiple rounds. The round numbers are coded in the file name. For example, 
<ul>

<li>xxx_0_Read_distribution.txt: The input read distribution (specified by --d).</li>
<li>xxx_0_GEM_events.txt: GPS Events used to re-estimate the read distribution of current dataset.</li>
<li>xxx_1_Read_distribution.txt: The read distribution estimated from xxx_0_GEM_events.</li>
<li>xxx_1_GEM_events.txt: GPS Events predicted using xxx_1_Read_distribution.</li>
<li>xxx_1_KSM/PFMs.txt: Motifs discovered using GPS events.</li>
<li>xxx_2_Read_distribution.txt: The read distribution estimated from GPS events.</li>
<li>xxx_2_GEM_events.txt: GEM Events predicted using xxx_2_Read_distribution and xxx_1_KSM/PFMs motifs.</li>
<li>xxx_2_KSM/PFMs.txt: Motifs discovered using GEM events.</li>
</ul>
Due to the large variability of datasets, the refined empirical spatial
distribution may become too noisy because too few events are used to
estimate the distribution. Running GEM further may make the empirical
distribution even worse. Therefore, GEM save the output files from each
round, and allow the user to check the spatial distributions to decide
which one is the best. To facilitate this process, GEM outputs an image
to plot all the read distribution curves
(xxx_All_Read_Distributions.png). Ideally, the read distribution of
later rounds should be smooth and similar to that of round 0. If so,
user can use the event predictions from these rounds. If that is not
the case, we would suggest to use the round 0 prediction results.<br>
<br>

<span style="font-weight: bold;">Some best practices?</span><br>
<ul>
<li>For each GEM result set, create a folder with the same name as <code>--out</code>. Then run GEM within this folder. This will be useful if later you need to write scripts to process multiple GEM results.
<li>After each GEM run, check the PNG image file that plot all the read distribution curves (xxx_All_Read_Distributions.png). If the estimated curves are not smooth, you may have too few binding events to estimate from.
<li>If you get something unexpected (or an error), check the command line options, make sure that each option starts with "--", and there is a space between each option/value pair.
</ul>
<span style="font-weight: bold;">Missed some events?</span><br>
Sometimes an event may be found by GEM, but it is not reported to the GEM_events.txt file. It may be reported in the _GEM_insignificant.txt file if it does not pass the statistical test. Or it could be in the _GEM_filtered.txt file if the shape of the binding event is too far away from the empirical read distribution. You can add <code>--nf</code> flag to turn off filtering.
<p>
<span style="font-weight: bold;">Can GEM/GPS process paired-end data?</span><br>
Yes. As version 2.4, GEM/GPS can process paired-end SAM/BAM format data (using the same <code>--f SAM</code> option) by treating each mate-pair as two single-end reads. This works well, but it is not optimal. We are developing a new version that explicitly models the paired-end data, which should give more accurate results.
<p>
<span style="font-weight: bold;">What if GEM finds wrong motif?</span><br>
Clearly, GEM's binding call accuracy is depending on finding the correct primary motif. Some times a co-factor motif may be more statistically significant in the data, and it is subsequently used to direct the binding calls. There are several alternative options to try:
<ul>
<li>If you know the consensus motif of the TF, use <code>--seed</code> option to set a starting k-mer for the motif discovery process.
<li>You may want to try some different <code>k</code> values, or different sequence window size <code>--k_win</code>.
<li>Try to use a different negative sequence set, for example <code>--k_neg_dinu_shuffle</code> option will use di-nucleotide shuffled sequences as negative set, instead of the default set taking from 300bp away from binding sites. This may be useful for TF that binds in CG-rich promoter regions (e.g. SP1).
</ul>
<p>


<span style="font-weight: bold;">Multi-condition v.s. Multi-replicates
</span><br>
GEM can analyze binding data from multiple conditions (time points)
simultaneously. The user need to give them different names, for
example, <code>-&#8211;exptCond1 CTCF_cond1.bed -&#8211;exptCond2
CTCF_cond2.bed</code>.<br>
<br>
For multiple replicates of same condition, you can specify multiple
replicates as separate files, for example,
<code> -&#8211;exptCond1 CTCF_cond1_rep1.bed -&#8211;exptCond1
CTCF_cond1_rep2.bed</code> (note that they need to have the same
name). GEM will combine the replicates as one large dataset for
analysis. <br>
<br>

<b>Read filtering and event filtering</b><br>
PCR amplification artifacts typically manifest as the observation of
many reads mapping to the exact same base positions. These artifacts
are quite variable and dataset-specific. Therefore, a generic approach
to exclude those regions might result in the loss of true events.
<br><br>
GEM implements an event filtering method by comparing the read
distribution of the predicted event to the expected event read
distribution. A shape deviation score (IPvsEMP field) is computed using
Kullback&#8211;Leibler divergence (see method section 2.6 of GPS paper). A
higher score means the event is more divergent from the expected read
distribution, hence more likely to be artifact or noise. A cutoff score
can be specified by user to filter out spurious events using option (<code>--sd</code>).
GEM also excludes events with less than 3 fold enrichment (IP/Control). GEM
reports the filtered events, hence allows the user to verify and adjust
cutoff threshold for a particular dataset. The shape deviation filter
is on by default, but can be turned off using option (<code>--nf</code>).
<br><br>
In addition, GEM also applies a Poisson filter for abnormal high read count at a base position. For each base, we obtain an average read count by estimating a Gaussian Kernel density (with std=20bp) on the read counts of nearby base positions (excluding the base of interest). The estimated value is used to set Lambda parameter of Poisson distribution. The actual read count value is then set to the value corresponding to p-value=0.001 if it is larger.
<!--
In addition, GEM also allows the user to set a cutoff
value for the maximum read count per base position. The cutoff value
can be estimated automatically using a Poisson model, or can be set
manually by the user (<code>--mrc</code>). -->


<h3 style="color: RoyalBlue;"><a class="anchor" id="contact"></a>Contact</h3>

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-17612319-2', 'auto');
  ga('send', 'pageview');

</script>

<p>Contact Yuchun Guo (yguo at mit dot edu) with any problems, comments, or suggestions.</p>
<p>Sign up for <a href="https://lists.csail.mit.edu/mailman/listinfo/gps">GPS
mailing list</a> to receive emails related to GEM/GPS updates,
release, etc. </p>

This software is for research use only. 


</td>

<td valign="top">


<table align="right" border="0" cellspacing="5">
<tbody>
<tr>
<td> <!--img style="width: 200px; height: 133px; float: right;" alt="GPS logo" src="http://cgs.csail.mit.edu/gps/images/gps_logo_small.png"-->
</td>
</tr>
<tr>
<td>
<table style="font-family: arial,sans-serif; font-size: 8pt; background-color: lightgrey;" border="0" frame="box" width="100%">
<tbody>
<tr>
<td align="center"> <b>Contents</b></td>
</tr>
<tr>
<td align="left" width="120">
<ul>
<li><a href="#download">Download</a></li>
<li><a href="#requirements">System
requirements</a></li>
<li><a href="#readdistrib">Read
distributions</a></li>
<li><a href="#io">Input and
output</a></li>
<li><a href="#examples">Examples</a></li>
<li><a href="#options">Command-line
options</a></li>
<li><a href="#QA">Q &amp; A</a></li>
<li><a href="#contact">Contact</a></li>
</ul>
</td>
</tr>
<tr><td>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>


<a href="images/ClustrMap_20150316/ClustrMap_20150316.htm"><img src="images/cgs.csail.mit.edu-gps--thumb.jpg" style="border:0px;" alt="Locations of visitors to this page, up to Mar 16, 2015" title="Locations of visitors to this page, up to Mar 16, 2015" id="clustrMapsImg_www4" />
</a>
<p>
<a href="https://www.google.com/analytics/web/?hl=en#report/visitors-geo/a17612319w97365702p101490300/%3F_r.tabId%3Dexplorer%26explorer-table-dataTable.sortColumnName%3Danalytics.visits%26explorer-table-dataTable.sortDescending%3Dtrue%26explorer-table.plotKeys%3D[]/">Google Analytics</a>
</td></tr></tbody>
</table>
</td>
</tr>
</tbody>
</table>

</td></tr></table>


</div>  <!-- end gem page -->
<!-- ***************  -->

      
    </body>
</html>

