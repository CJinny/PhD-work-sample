# Summary
  - **<i>Ufo1-1</i> project:**
    This is the main focus in <b>Dr. Surinder Chopra's lab</b>, we used a RNA-seq k-mer approach to identify candidate gene responsible for the natural mutation. Pacbio long-read sequencing analysis further revealed that a transposable element insertion was present in the candidate gene in <i>Ufo1-1</i> mutant. Subsequent transgenic experiment confirmed the canidate gene as <i>ufo1</i> gene. For more detail please check out our [publication](http://www.plantcell.org/content/30/12/3006.abstract).
  
  - **ITS analysis from Sorghum soil samples:** We collected soil samples from Sorghum lines with variable amount of flavonoid compound (3-deoxyanthocyanidin) and performed 16S and ITS sequencing. We hope to understand how flavonoid level affects the microbiome community. This work is done in collaboration with <b>Dr. Mary Ann Bruns's lab</b>. I learned the analysis pipeline from <b>Mara Cloutier</b>.
  
# Description of analysis
  - **RNA-seq**: Raw reads from <i>Ufo1-1</i> mutant and wild-type were generated by Illumina HiSeq2500 (single-end) and trimmed by <b>Dr. Tzuu-fen Lee</b> to remove adaptors, rRNA reads etc. In this repo, I first performed alignment using [<b>STAR-2pass</b>](https://github.com/CJinny/PhD-work-sample/blob/master/STAR_2pass.pbs), and then generated [<b>StringTie assemblies</b>](https://github.com/CJinny/PhD-work-sample/blob/master/stringtie.pbs). I used the merged StringTie assembly to perform [<b>Kallisto pseudoalignment</b>](https://github.com/CJinny/PhD-work-sample/blob/master/kallisto.pbs) to generate count data. Notice that the sequencing pipeline is different from the original paper (HiSAT2 => htseq on RefGen V3). The reason I do this is because I want to include new transcripts not found in the reference genome. In addition, the STAR-2pass pipeline is also recommended for GATK variant calling and I also used this alignment to analyze differential exon usage using DEXSeq (not included in this repo). In this [<b>Jupyter notebook</b>](https://github.com/CJinny/PhD-work-sample/blob/master/rnseq_kmers_extraction.ipynb) I also include a demo of k-mer analysis 
  
  - **PacBio**: Structural variations were identified with [<b>NGMLR](https://github.com/philres/ngmlr) => [Sniffles</b>](https://github.com/fritzsedlazeck/Sniffles) pipeline. The PBS script can be found [<b>here</b>](https://github.com/CJinny/PhD-work-sample/blob/master/pacbio_SV.pbs).
  
  - **small RNA-seq**: Again, raw reads from <i>Ufo1-1</i> mutant and wild-type were trimmed by <b>Dr. Tzuu-fen Lee</b>. I used ShortStack to conduct de-novo small RNA loci annotation as well as count accumulation at given loci (e.g. differentially methylated regions and transposable elements). I used edgeR to identify differentially expressed small RNA loci (not included in this repo) and perform a [co-occupancy] analysis to statistically measure how likely small RNA loci are co-occupied with genic/transposable element features. Some data visualization can be found [<b>here</b>](http://htmlpreview.github.io/?https://github.com/CJinny/PhD-work-sample/blob/master/sRNA_analysis_report.html).

  - **ITS**: The analysis pipeline was adopted from [<b>here</b>](http://benjjneb.github.io/dada2/tutorial.html). I'm still learning. The Markdown report can be found here [<b>here</b>](http://htmlpreview.github.io/?https://github.com/CJinny/PhD-work-sample/blob/master/sRNA_analysis_report.html).

