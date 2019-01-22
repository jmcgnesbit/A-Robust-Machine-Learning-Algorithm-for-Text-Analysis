# Faster-than-fast NMF using random projections and Nesterov iterations
## by Farouk Yahaya, Matthieu Puigt, Gilles Delmaire, and Gilles Roussel


-----


Random projections have been recently implemented in Nonnegative Matrix Factorization (NMF) to speed-up the NMF computations, with a negligible loss of performance. In this paper, we investigate the effects of such projections when the NMF technique uses the fast Nesterov gradient descent (NeNMF). We experimentally show that structured random projections significantly speed-up NeNMF for very large data matrices. 

### Download
[Download the Matlab code of the proposed NeNMF method with random projections](https://gogs.univ-littoral.fr/puigt/Faster-than-fast_NMF/src/master/Code).
The code was written by Farouk YAHAYA and is maintained by Farouk YAHAYA, Matthieu PUIGT, Gilles DELMAIRE, and Gilles ROUSSEL (firstname.LASTNAME [at] univ-littoral.fr).

It contains several functions and **folders**:

* README.m
* add_paths.m
* Demo.m
* data_simuilation.m
* plots.m
* stop_rule.m
* **algorithms**
    * VANILLA_NeNMF.m
    * RPI_NeNMF.m
    * RSI_NeNMF.m 
* **output**
* **synthetic_data**


### Reference
If you use this code for research or educational purpose, please cite:

> F. Yahaya, M. Puigt, G. Delmaire, and G. Roussel, *"[Faster-than-fast NMF using random projections and Nesterov iterations,](https://hal.archives-ouvertes.fr/hal-01859713)"* in Proc. of  iTWIST: international Traveling Workshop on Interactions between low-complexity data models and Sensing Techniques, Marseille, France, November 21-23, 2018. 
> 

### Support

For any suggestions or questions about this code, please contact: Farouk.Yahaya [at] univ-littoral.fr and Matthieu.Puigt [at] univ-littoral.fr.