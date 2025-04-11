# Testing isomorphism between tuples of subspaces

Given two tuples of subspaces, or more precisely, orthonormal bases for those subspaces, can you tell whether the tuples are isomorphic? Some Matlab code to test for isomorphism are posted here.

Code to compute an injective invariant (normalized Gram) for tuples in the Grassmannian of planes in R^d modulo O(d) given a Gram matrix of orthobases for those planes is norgPlanesGram.m.

Code to compute words which may be evaluated on a generated H^*-algebra to determine an injective invariant is canBasMatAlg.m. 

Code to compute traces of the evaluated words to determine isomorphism of n-tuples of d×d matrices generating H^∗-algebra is traceWords.m.

For more information, see Algoritm 1, Algorithm 2, and Lemma 11 in "Testing isomorphism between tuples of subspaces" by Emily J. King, Dustin G. Mixon, and Shayne Waldron. https://arxiv.org/abs/2105.03448


