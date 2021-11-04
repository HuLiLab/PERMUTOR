Pipeline: PERsonalized MUTation evaluaTOR (PERMUTOR)
Version: 1.0
Author: Taylor Weiskittel
Maintainer: Taylor Weiskittel <weiskittel.taylor@mayo.edu>

R version 4.1.1 (2021-08-10) -- "Kick Things"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)
REQUIRED PACKAGES:
plyr 1.8.6
igraph 1.2.6
rlist 0.4.6.2
rapport 1.1
stringi 1.7.3
stringr 1.4.0

INPUTS:
FCRNAseq.csv= The absolute value of the log2 fold change for each gene(rows) between disease tissue and matched normal tissue per patient(columns)
patmutations.csv= A matrix annotationg for each gene(rows) and each patient(column) if a gene is mutated (mutated=1) or not (not mutated=0)
diseasecontext.csv= A single column annotating which genes (rownames) have strong evidence of disease activity (=1) or medium evidence of disease activity (=2) in the literature 
generalppi.csv= Initial protien-protien interaction network in edgelist format


OUTPUTS:
x allpatientpaths.Rdata= an R list object showing all detected shortest paths for patient x between that patient's mutated genes
x randscore n.Rdata= an R vector containing the random path scores of lenght n for patient x
x allrandpaths.Rdata= an R list object containing all generated random paths for x patient 
x allrandscores.Rdata= an R vector combining the values for patient x from the randscore files
x pvalues.Rdata= an R vector containing all empierical p-values for all detected paths, corresponds to the order of the allpatientpaths
x individualnetworkpaths.Rdata= an R list object showing all significant shortest paths included in patient x individualized disease module
x indedges.Rdata= an R data frame with patient x individualized disease module in edgelist format
x indnodes.Rdata= an R vector with all of patient x nodes in their individualized disease module 
x IDGI.Rdata= an R data frame detailing patient x IDGI scoring components and total scores
x ITGI.Rdata= an R data frame detailing patient x ITGI scoring components and total scores 
x ITGIcombinations.Rdata= an R data frame detailing patient x ITGI combinatorial scoring components and total scores 
x patientcharacteristics.Rdata= an R data frame with summary characteristics of patient x and their individualized disease module