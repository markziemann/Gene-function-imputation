# Overview
The candidate will perform a search of the recent literature (2016-current) on the use of machine learning, deep learning and artificial intelligence technology in inferring genetic networks and gene functions from coexpression and other biological data [1].
The candidate will provide a 2-3 page document summarising what they have found and how it could be applied to generate knowledge from large scale gene expression data such as ARCHS4 and DEE2 [2,3]. 

# Literature Report 

Majority of the methods in literature that has been found and summarised below are recent studies that have been applied to protein function imputation to date. Despite this, such findings for proteins also show potential for application to the functional prediction of other molecules, such as that which is the focus of our study - lncRNA.

## Deep learning methods
Deep learning methods have been found to cooperate well with problems involving big data (Bonetta & Valentino 2019).
(Trofimov et al. 2020)
* **Factorised embeddings (FE) model** - self-supervised algorithm that simultaneously learns gene and sample representation spaces through tensor factorisation.
* Despite having focused on protein-coding genes, this method may have the potential to be applied to non-coding genes and lncRNA.
* Ran on RNA-Seq data.
* Learned embedding spaces capture biologically meaningful information, including that regarding gene-gene co-expression and gene function.

(Ioannidis, Marques & Giannakis 2019) **Graph Neural Networks**
* Semi-supervised learning (SSL) over graphs
* The inadequacy of representing the relatedness between nodal variables by a single graph has been recognised - which has been observed in many works involving SSL over graphs. 
* Association/connectivity graphs can appear different, depending on the cell type. This is an important fact to consider when modelling and predicting protein-to-protein interactions and their functions. 
* Hence, it is understandable that a single graph representation may fail in situations where nodes are involved in numerous relation types. 
* Just as how people have various relations such as familial, friendship, and professional bonds, this study shows the transition from single-relational to multi-relational graphs, and aimed to impute protein functions over multi-relational protein-to-protein networks.
* Used a Graph Residual Neural Network (GRNN) for the development of multi-relational graphs and prediction of protein functions. 
* Predicts unknown labels by mapping each node to a corresponding label, using the input feature matrix X. 
* PPI networks related two proteins via multiple cell-dependent relations, and protein classification predicted unknown functions of some proteins based on the functions of a small subset of them. With a known target function on a protein subset, known functions of proteins in X, and multi-relational protein networks, the idea was to predict whether an association existed between proteins in the unlabelled set and the target function.

(Yuan & Bar-Joseph 2019) **Convolutional neural network for co-expression (CNNC)**
* Supervised framework for gene relationship inference, which learns to distinguish between interacting, negative pairs, causal pairs, or any other gene relationship types that can be defined. 
* Each gene pair is represented as an image (histogram), with convolutional neural networks (CNNs) used to infer relationships between different expression levels encoded in the image. 
* CNNC was applied to a large amount of single-cell RNA-sequencing (scRNA-seq) data, as well as bulk RNA-seq data, to conduct various tasks, including interaction inference, causality inference, and functional assignment. 
* This model can be trained with any expression dataset, with better performance when it's provided more data.  
* In testing CNNC for functional assignments, the method was found to significantly outperform both DNN and GBA. 

## Machine learning methods
(Du et al. 2019)
* Applied concept of word embedding in NLP deep learning research to genes.
* A distributed representation of genes is produced from transcriptome-wide gene co-expression data.
* Embedding matrix solely trained from gene co-expression patterns (gene co-expression provides the 'context' for gene embedding), capturing functional gene relationships.
* Embedded gene factors heavy with gene co-expression patter information was shown to be useful in tasks such as predicting gene-gene interactions.
* Aimed to represent genes as vectors through gene embedding, which grouped similar genes into spatial clusters - similar genes mapped to similar vectors

**Random forest**
* Aggregation of multiple decision trees to produce a single final result (Liberman 2017).
* Notion of 'bagging' (Bonetta & Valentino 2019).

Novel methods of gene function imputation can also be created through the integrative effort of multiple techniques (Wekesa, Luan & Meng 2020). 
* Method: **CNPFP**
* This study aimed to identify cell cycle-specific proteins in yeast which are involved with differentially expressed proteins within the PPI network, and functionally annotate 3538 yeast proteins
* Based on and utilises differential co-expression analysis and neighbour-voting (NV) (also known as K-Nearest Neighbours) algorithm to predict protein function from yeast cell cycle gene expression and PPI data sets, with the degree of similarity between genes measured by biweight midcorrelation (BWM).
* Described a threefold process: cluster analysis through 'guilt-by-profiling', differential co-expression analysis, and the novel functional annotation of proteins through the neighbour-voting algorithm, with the exploitation of label correlations, genomic features and intrinsic information from the co-expression analysis conducted.
* Their integrative method was created in recognition that majority of protein function prediction algorithms do not make use of the diverse intrinsic information in feature and label spaces, concluding that the exploitation of intrinsic information from protein relationships improves the quality of predicting functions.

## Artificial intelligence methods
* NLP (natural linguistic processing) - potential applications, as shown by Du et al.

## Other methods to explore
* Logistic regression (machine learning)
* Matrix factorisation (collaborative filtering class)
* Deep Neural Networks (DNNs)
* SVM

## References
* Bonetta, R, Valentino, G 2019, 'Machine learning techniques for protein function prediction', Proteins: Structure, Function, and Bioinformatics, vol. 88, no. 3, p. 397, doi: 10.1002/prot.25832
* Du, J, Jia, P, Dai, Y, Tao, C, Zhao, Z, Zhi, D 2019 'Gene2vec: distributed representation of genes based on co-expression', BMC genomics, vol. 20, no. 82, doi: 10.1186/s12864-018-5370-x
* Ioannidis, VN, Marques, AG, Giannakis, GB 2019, 'Graph Neural Networks for Predicting Protein Functions', Proceedings of the IEEE 8th International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP), IEEE, Le gosier, Guadeloupe, pp. 221-225, <https://ieeexplore-ieee-org.ezproxy-b.deakin.edu.au/document/9022646#full-text-section>
* Liberman, N 2017, 'Decision Trees and Random Forests', towards data science, weblog post, 27 January, retrieved 26 July 2020, <https://towardsdatascience.com/decision-trees-and-random-forests-df0c3123f991>.
* Nguyen, G, Dlugolinsky, S, Bobak, M, Tran, V, Garcia, AL, Heredia, I, Malik, P, Hluchy, L 2019, 'Machine Learning and Deep Learning frameworks and libraries for large-scale data mining: a survey', Artificial Intelligence Review, vol. 52, p. 77-124, doi: 10.1007/s10462-018-09679-z
* Trofimov, A, Cohen, JP, Bengio, Y, Perreault, C, Lemieux, S 2020, 'Factorized embeddings learns rich and biologically meaningful embedding spaces using factorised tensor decomposition', Bioinformatics, vol. 36, no. 1, p. i417-i426, doi:10.1093/bioinformatics/btaa488
* Wekesa, JS, Luan, Y, Meng, J 2020, 'Predicting Protein Functions Based on Differential Co-expression and Neighborhood Analysis', Journal of Computational Biology, Preprint, viewed 27 July 2020, doi: 10.1089/cmb.2019.0120
* Yuan, Y, Bar-Joseph, Z 2019, 'Deep learning for inferring gene relationships from single-cell expression data', PNAS, vol. 116, no. 52, p. 27151-27158, doi: 
10.1073/pnas.1911536116


# References

Ballouz S, Verleyen W, Gillis J. Guidance for RNA-seq co-expression network construction and analysis: safety in numbers. Bioinformatics. 2015;31(13):2123‐2130. doi:10.1093/bioinformatics/btv118

Lachmann A, Torre D, Keenan AB, et al. Massive mining of publicly available RNA-seq data from human and mouse. Nat Commun. 2018;9(1):1366. Published 2018 Apr 10. doi:10.1038/s41467-018-03751-6

Ziemann M, Kaspi A, El-Osta A. Digital expression explorer 2: a repository of uniformly processed RNA sequencing data. Gigascience. 2019;8(4):giz022. doi:10.1093/gigascience/giz022
