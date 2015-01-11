# Paper08-ISMB-ProteinComplex-GraphModule
Code Rep for 2008-ISMB paper: "Protein Complex Identification by Supervised Graph Clustering , Bioinformatics 2008, 24(13), i250-i268"


<BR>

Please cite: 

@article{qi2008protein,
  title={Protein complex identification by supervised graph local clustering},
  author={Qi, Yanjun and Balem, Fernanda and Faloutsos, Christos and Klein-Seetharaman, Judith and Bar-Joseph, Ziv},
  journal={Bioinformatics},
  volume={24},
  number={13},
  pages={i250--i268},
  year={2008},
  publisher={Oxford Univ Press}
}

<BR>

--------------------

Abstract

Motivation: Protein complexes integrate multiple gene products to coordinate many biological functions. Given a graph representing pairwise protein interaction data one can search for subgraphs representing protein complexes. Previous methods for performing such search relied on the assumption that complexes form a clique in that graph. While this assumption is true for some complexes, it does not hold for many others. New algorithms are required in order to recover complexes with other types of topological structure.
 
Results: We present an algorithm for inferring protein complexes from weighted interaction graphs. By using graph topological patterns and biological properties as features, we model each complex subgraph by a probabilistic Bayesian Network (BN). We use a training set of known complexes to learn the parameters of this BN model. The log-likelihood ratio derived from the BN is then used to score subgraphs in the protein interaction graph and identify new complexes. We applied our method to protein interaction data in yeast. As we show our algorithm achieved a considerable improvement over clique based algorithms in terms of its ability to recover known complexes. We discuss some of the new complexes predicted by our algorithm and determine that they likely represent true complexes.
 

 <BR>

--------------------


Supplement web @ http://www.cs.cmu.edu/~qyj/SuperComplex/
 
<br>

Talk Slide @ http://www.cs.cmu.edu/~qyj/SuperComplex/sharefiles/complex-Talk-online.pdf

<br>
Paper PDF @ http://bioinformatics.oxfordjournals.org/content/24/13/i250.abstract?etoc


 <BR>

--------------------

Code 

·        Matlab Implementation is available @ download  
·        Plan to convert the code into one Cytoscape plug-in, please check back for updates.



<BR>

--------------------

Data
 
·        The MIPS complexes we used as references were extracted from two files:  complexCat  complexScheme  (high-throughput complexes not used)
·        The TAP complexes we used as references were extracted from the paper’s supplementary data file; we put one copy here for the reader’s convenience
 
 