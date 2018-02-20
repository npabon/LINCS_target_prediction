# Predicting drug targets from gene expression perturbation signatures<br>

### (Very) Brief background on drug discovery

Many medicinal drugs are small molecules that bind to specific proteins in the body and inihibit their function. This can have a variety of network-level effects depending on the target protein. In many cases the disease involves protein interactions that are happening "too much", and the inhibitor, by binding to the protein, has the effect of "slowing down" the interaction or preventing it entirely.

Historically, drugs have been discovered by chance. High throughput screens test thousands of drugs and those that produce the desired outcome in *in-vitro* assays are selected for further development. Often, drugs' specific target protein(s) are not known a-priori, and this can lead to side-effects and toxicity later down the drug developemnt pipeline when the compounds are tested in live animals. Being able to predict which specific proteins a bioactive compound will bind to and inhibit has immense potential to accelerate drug discovery and drug repurposing.

### Leveraging "big data" to predict drug-target interactions

Today we're witnessing an explosion of publicly available biological data through initiatives like the 1000 Genomes Project, The Cancer Genome Atlas, and the NIH's Library of Integrated Network-Based Cellular Signatures (LINCS) Program. The LINCS initiative generates and makes public gene expression data that indicates how different types of cells respond to various genetic and environmental perturbations, including drugs and other bioactive small molecules. 

Gene expression regulation in human cells is a complex and noisy process, making it difficult to determine the protein target of a drug by looking at which genes in the cell are up/down regulated as a result of drug treatment. However, here we hypothesize that a cell treated with an small molecule inhibitor for a given protein should produce a similar gene expression profile to that of an untreated cell in which the gene coding for that protein has been silenced. See the figure below for a visualization of this idea.

![](figures/hypothesis.jpg?raw=true)

LINCS contains gene expression profiles for over 20,000 small molecule/drug treatments and over 20,000 gene silencing, or "knockdown" experiments. We know the protein targets of many of the drugs that are tested, and some of these proteins were knocked down in separate experiments. So, we can test our hypothesis by triyng to build a machine learning classifier to predict whether a drug and a protein interact by looking at the gene expression profiles from their corresponding drug treatment and gene knockdown experiments.

If our classifier works, we can use it to predict targets for new drugs by comparing the drug's expression profiles against all the knockdown profiles in LINCS.
