This toolbox allows the implementation of the following diffusion-based 
clustering algorithms on synthetick and real datasets included in the 
repository:

 -  Learning by Unsupervised Nonlinear Diffusion (LUND)
 -  Multiscale Learning by Unsupervised Nonlinear Diffusion (M-LUND)
 -  Spatially Regularized Diffusion Learning (SRDL)
 -  Multiscale Spatially Regularized Diffusion Learning (M-SRDL)

This package can be used to generate experiments in the following articles:

   - Murphy, James M and Polk, Sam L., 2021. A Multiscale Environment for
     Learning By Diffusion. (In Preparation)
   - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-Spatial 
     Diffusion Geometry\\ for Hyperspectral Image Clustering. (In Review)

The following scripts (in the Experiments folder) generate the relevant 
experiments:

   - M_LUND_demo.m implements M_LUND on synthetic data and Salinas A.
   - Benchmark.m generates Table 3, in which comparisons of M-LUND against 
     MMS, HSC, and SLC on benchmark datesets are provided. 
   - M_SRDL_demo.m compares M-SRDL against M-LUND on the Salinas A HSI.
   
All necessary datasets are contained in this repository, so no additional 
data downloads are necessary. All benchmark datasets were obtained from 
the UCI Machine Learning Repository: 

   - https://archive.ics.uci.edu/ml/index.php

The real HSI (Salinas A) was from the 2000 IEEE Data Fusion Contest data. 
This data is publically available via IEEE: 

   - http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes

To run comparisons against MMS clustering, the following toolboxes, written 
by Zhijing Liu and Mauricio Barahona, must be added to one's path:

   - https://github.com/barahona-research-group/GraphBasedClustering
   - https://wwwf.imperial.ac.uk/~mpbara/Partition_Stability/

Users are free to modify the Multiscale Diffusion Clustering toolbox as they 
wish. If you find it useful and use it in any publications, please cite the 
following papers:

   - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
     Diffusion. Journal of Machine Learning Research, 20(160), pp. 1-56. 
     2019.
   - Murphy, James M., and Mauro Maggioni. "Spectral-spatial diffusion 
     geometry for hyperspectral image clustering." IEEE Geoscience and 
     Remote Sensing Letters (2019).
   - Murphy, James M and Polk, Sam L., 2020. A Multiscale Environment for
     Learning By Diffusion. arXiv preprint arXiv: XXXXXXX.
   - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-Spatial 
     Diffusion Geometry for Hyperspectral Image Clustering. (In Review)

Please write with any questions: samuel.polk@tufts.edu

(c) Copyright Sam L. Polk, Tufts University, 2021. 
