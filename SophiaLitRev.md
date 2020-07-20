# Overview
The candidate will perform a search of the recent literature (2016-current) on the use of machine learning, deep learning and artificial intelligence technology in inferring genetic networks and gene functions from coexpression and other biological data [1].
The candidate will provide a 2-3 page document summarising what they have found and how it could be applied to generate knowledge from large scale gene expression data such as ARCHS4 and DEE2 [2,3]. 

# Literature Report 
## Deep learning
(Trofimov et al. 2020)
* **Factorised embeddings (FE) model** - self-supervised algorithm that simultaneously learns gene and sample representation spaces through tensor factorisation
* Despite having focused on protein-coding genes, this method may have the potential to be applied to non-coding genes and lncRNA
* Ran on RNA-Seq data
* Learned embedding spaces capture biologically meaningful information, including that regarding gene-gene co-expression and gene function.

## Machine learning 
(Du et al. 2019)
* Applies concept of word embedding in NLP deep learning research to genes
* A distributed representation of genes is produced from transcriptome-wide gene co-expression data
* Embedding matrix solely trained from gene co-expression patterns (gene co-expression provides the 'context' for gene embedding), capturing functional gene relationships.
* Embedded gene factors heavy with gene co-expression patter information was shown to be useful in tasks such as predicting gene-gene interactions.
* Aimed to represent genes as vectors through gene embedding, which grouped similar genes into spatial clusters - similar genes mapped to similar vectors

## Artificial intelligence 
* NLP (natural linguistic processing) - potential applications, as shown by Du et al.

## Other methods to explore
* Graph neural networks
* Matrix factorisation (collaborative filtering class)
* Deep Neural Networks (DNNs)

## References
* Trofimov, A, Cohen, JP, Bengio, Y, Perreault, C, Lemieux, S 2020, 'Factorized embeddings learns rich and biologically meaningful embedding spaces using factorised tensor decomposition', Bioinformatics, vol. 36, no. 1, p. i417-i426, doi:10.1093/bioinformatics/btaa488
* Du, J, Jia, P, Dai, Y, Tao, C, Zhao, Z, Zhi, D 2019 'Gene2vec: distributed representation of genes based on co-expression', BMC genomics, vol. 20, no. 82, doi: 10.1186/s12864-018-5370-x

# References

Ballouz S, Verleyen W, Gillis J. Guidance for RNA-seq co-expression network construction and analysis: safety in numbers. Bioinformatics. 2015;31(13):2123‐2130. doi:10.1093/bioinformatics/btv118

Lachmann A, Torre D, Keenan AB, et al. Massive mining of publicly available RNA-seq data from human and mouse. Nat Commun. 2018;9(1):1366. Published 2018 Apr 10. doi:10.1038/s41467-018-03751-6

Ziemann M, Kaspi A, El-Osta A. Digital expression explorer 2: a repository of uniformly processed RNA sequencing data. Gigascience. 2019;8(4):giz022. doi:10.1093/gigascience/giz022
