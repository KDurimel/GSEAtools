GENIAL : GENe set enrIchment AnaLysis <img src="https://github.com/KDurimel/GSEAtools/blob/master/logo.png" height="40">
===================

Authors: Kevin Durimel & Victor Mataigne

Affiliation: [University of Rouen (France)](http://www.univ-rouen.fr/)

Web application now availaible (still a work in progress) : https://genial.shinyapps.io/GenialGSEA/

WARNING : this script was developed under the context of a bioinformatic project (20h) and will be not long term supported!

<i class="icon-file"></i>  Path organization
=================
	.
	│
	│--bin : windows executable files
	│
	│--dependencies : R package dependencies list (JSON)
	│
	│--log : log files
	│
	│--example_data : experiment specific raw_data (currently --> Workflow input) 
	│
	│--results : results examples (one path generated per result, Day_Month_HH_MM-SS-DD_YY format)
	│
	│--src : R scripts

     

The project 
==================================================================
From a input dataset of our choice from the literature, we had to implement R packages to realize all three GO,metabolic        pathways, and protein functional domains enrichment annotations, implementing robust statistical procedures. Outputs will be both charts (pie, bar plot, graphs, tracks) and texts results (tables). This workfow must be interfaced with Shiny, a R package developed by RStudio
    
SHINY INTERFACE OVERVIEW :
-----------------------------

Divided in 5 main sections :

	-Load and explore raw data
	-Vulcano plot
	- GO enrichment
		* GO classification (barplot)
		* Ability to set a threshold for go annotations level / ontologies
		* GO three browsing
	- KEGG enrichment
		* KEGG enrichment results overview
		* KEGG enrichment results browsing by typing target pathway IDs on a web form
	- Protein domains
		* Hypergeometric test
		* Ability to add metatada on the results :Panther_Family_ID, Panther_Family_Description,EnsemblTranscriptID, EnsemblPeptideID, PfamID, PfamSTART, PfamEND

SHINY INTERFACE OVERVIEW :
-----------------------------

Same as the Shiny interface.


STANDALONE SCRIPT :
---------------------

Rscript GENIAL-NOShiny.R

This is the GUI interfaced R script (i.e Non-Shiny), containing more features that the interfaced one.

## Usage (Linux/Mac OS X)

If necessary, set up the scripts permissions inside the src path :

	chmod +x src/*

Add the scripts to your bashrc (/home/username/.bashrc) :

	export PATH=$PATH:src/
	
Then you can run it as a shell command :

	GENIAL-NOShiny.R


Dependencies
============

GENIAL has been developped with R 3.3.2

## External dependencies

* [shiny](https://shiny.rstudio.com/) - Version 1.0.1
* [shinythemes](https://github.com/rstudio/shinythemes) - Version 1.1.1
* [ggplot2](https://github.com/tidyverse/ggplot2) - Version 2.2.1
* [clusterProfiler](https://github.com/GuangchuangYu/clusterProfiler) - Version 3.2.11
* [GO.db](https://bioconductor.org/packages/release/data/annotation/html/GO.db.html) - Version 3.2.0
* [org.Hs.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) - Version 3.4.0
* [topGO](http://bioconductor.org/packages/release/bioc/html/topGO.html) - Version 2.26.0
* [DT wrapper](https://cran.r-project.org/web/packages/DT/index.html) - Version 0.2
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) - Version 2.30.0
* [MASS](https://github.com/cran/MASS) - Version 7.3-45



BE CAUTIOUS :
---------------------------
WARNING 1 : This workflow needs some containerization (Conda, Packrat... )to work properly / ensure repdroducibility, which was not done under the context of this very short project. You need to install exactly the same R and package versions as cited in the "external dependencies" sections.
Due to this, the [Web application](https://genial.shinyapps.io/GenialGSEA/) may stop working properly at any time. This project is not long term supported and any future bugs will be not fixed.

WARNING 2 : Due to Shiny technical limitations, tasks launched in different tabs are not parallelised or processed asynchronously.

WARNING 3 : This application works only with data formated like "human_retine.txt" example file (columns names and content).

WARNING 4 : A lot of APIs are used --> The better your internet connection is the faster the application runs.


Example dataset obtained from :
-----------------------------

	Sara R. Savage, Colin A. Bretz, and John S. Penn ; RNA-Seq Reveals a Role for NFAT-Signaling in Human 
	Retinal Microvascular Endothelial Cells Treated with TNFα ; PLoS One. 2015
	
	doi:  10.1371/journal.pone.0116941
	
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4305319/

***********************************************************************************************************
Kevin Durimel, Victor Mataigne - M2.1 BIM - University of Rouen, 2017.

	Standalone script - Web app : Kevin Durimel
	Shiny interfacing : Victor Mataigne
