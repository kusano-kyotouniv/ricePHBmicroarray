Followings are the informatics methods used to study the microarray result from an line of experiment to know physiological response of rice callus in which the biosynthesis pathway genes for polyhydroxybutanoate was introduced. The order of method is as follow for each experiment refered to as Figures in the paper. 

Fig4A, hierarchical clustering
	1 -> 2 

Fig6 and Fig4B-D, scatter plotting
	1 -> 3　-> you can use mouse to use the GUI

Fig5, co-expression network
	1 -> 4 -> 5

*Note that the graphical results may looks different between every times you get it by run the script. The networks you may get mean equal although the angle and length of lines may different between the trials. This difference is caused from the random function contained in the rayout algorism. 


＊＊＊
Following 1-5 is explanation of responsible methods.


1. make_GOenrich_gene.pl

The perl script process the raw microarray data and output the result data of GO enrichment analysis. The data is text format but should be processed by following methods to get each graphic files.

input files required:
	IRGSP1.txt	annotation information. version of rice annotation project 2022.9
	raw1.txt	raw data of microarray result
	raw2.txt	raw data of microarray result
	raw3.txt	raw data of microarray result
	raw4.txt	raw data of microarray result
	* Note that it may include data out of this line of experiment.

Files may be produced:
	tamakiarray_hclust.txt
		inpur data for hierarchical clustering for hclust.r 

	tamakiarray_gene_RNAseqViewer2.txt
		gene expression data for scatterplot with P-value calculated for volcano plot, loadble by RNAseqViewer6_for_mac.pl under "RSEMdata" button.

	GOenrich_molfunc.txt
	GOenrich_celcomp.txt
	GOenrich_bioproc.txt
		formatted data to produce graphics file of network graph. Divided into three categories of GO terms, Molecular Function, Biological Process, and Cellular Component.



2. hclust.r

R script describing the setting of hierarchical clustering.

input files required:
	tamakiarray_hclust.txt	generated in method 1.

Files may be produced:
	hclust.pdf		result.


3. RNAseqViewer6_for_mac.pl

GUI software to view transcriotome, arranged for this data set. It contains function to draw scatter plot of gene expression data and output the .PNG format graphics file.

input files required:
	RNAseqViewer_settings.txt			inner data required, including the column
	tamakiarray_gene_RNAseqViewer2.txt	generated in method 1.
		You can get image file of scatter plot by push the [Load RSEM] button and set the input file above.

Files may be produced:
	RNAseqViewer_scatplot.png			result.


4. make_bigraph2.pl

Perl script to get co-occurrence data to generate the network graphics file in bigraph format required by igraph library in R.

input files required:
	transcripts.gff		annotation information. version of RAP 2022.9
	GOenrich_molfunc.txt	generated in method 1
	GOenrich_celcomp.txt	generated in method 1
	GOenrich_bioproc.txt	generated in method 1

Files may be produced:
	bigraph_molfunc.txt			co-occurrecne data in Molecular Function category
	bigraph_molfunc_rare.txt	GOs with only one co-occurrence target 
	bigraph_molfunc_orphan.txt	GOs with no co-occurrence target
	bigraph_bioproc.txt			same data set in Biological Process.
	bigraph_bioproc_rare.txt		
	bigraph_bioproc_orphan.txt
	bigraph_celcomp.txt			same data set in Biological Process.
	bigraph_celcomp_rare.txt
	bigraph_celcomp_orphan.txt

	vertex_molfunc.txt			color setting data fro each.
	vertex_molfunc_rare.txt
	vertex_molfunc_orphan.txt
	vertex_bioproc.txt
	vertex_bioproc_rare.txt
	vertex_bioproc_orphan.txt
	vertex_celcomp.txt
	vertex_celcomp_rare.txt
	vertex_celcomp_orphan.txt


5. run_igraph.sh

shell script to run the R scripts generateing network graph image from the bigraph data output from 4.

input files required:
	igraph_bioproc_orphan.r		R scripts process the bigraph into graphics file.
	igraph_bioproc_rare.r	
	igraph_bioproc.r
	igraph_celcomp_orphan.r
	igraph_celcomp_rare.r
	igraph_celcomp.r
	igraph_molfunc_orphan.r
	igraph_molfunc_rare.r
	igraph_molfunc.r
	4. で生成したファイル	

Files may be produced:
	image_bioproc_orphan.pdf	result.
	image_bioproc_rare.pdf
	image_bioproc.pdf
	image_celcomp_orphan.pdf
	image_celcomp_rare.pdf
	image_celcomp.pdf
	image_molfunc_orphan.pdf
	image_molfunc_rare.pdf
	image_molfunc.pdf

	
