This toolbox allows the implementation of the following diffusion-based 
clustering algorithms on synthetick and real datasets included in the 
repository:

 -  Learning by Unsupervised Nonlinear Diffusion (LUND)
 -  Multiscale Learning by Unsupervised Nonlinear Diffusion (M-LUND)
 -  Spatially Regularized Diffusion Learning (SRDL)
 -  Multiscale Spatially Regularized Diffusion Learning (M-SRDL)

This package can be used to generate experiments in the following articles:

   - Murphy, James M and Polk, Sam L., 2021. A Multiscale Environment for
     Learning By Diffusion. arXiv preprint arXiv: XXXXXXX.
   - Polk, Sam L. and Murphy James M., 2021. Multiscale Clustering of 
     Hyperspectral and Spectral-Spatial Diffusion Geometry. arXiv preprint 
     arXiv: XXXXXXX.

The following scripts (in the Experiments folder) generate the experiments 
in the manuscript above:

   - M_LUND_demo.m implements M_LUND on synthetic data and Salinas A.
   - Table3.m generates Table 3, in which comparisons of M-LUND against 
     MMS, HSC, and SLC on benchmark datesets are provided. 
   - M_SRDL_demo.m compares M-SRDL against M-LUND on the Salinas A HSI.
   
All benchmark datasets were obtained from the UCI Machine Learning 
Repository: 

   - https://archive.ics.uci.edu/ml/index.php

The real HSI (Salinas A) was from the 2000 IEEE Data Fusion Contest data. 
This data is publically available via IEEE: 

   - http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes

Users are free to modify this toolbox as they wish. If you find it useful
and use it in any publications, please cite the following papers:

   - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
     Diffusion. Journal of Machine Learning Research, 20(160), pp. 1-56. 
     2019.
   - Murphy, J.M., Maggioni, M. Unsupervised Clustering and Active Learning 
     of Hyperspectral Images with Nonlinear Diffusion. IEEE Transactions on 
     Geoscience and Remote Sensing, 57(3), pp. 1829-1845. 2019.
   - Murphy, James M and Polk, Sam L., 2020. A Multiscale Environment for
     Learning By Diffusion. arXiv preprint arXiv: XXXXXXX.
   - Polk, Sam L. and Murphy James M., 2021. Multiscale Clustering of 
     Hyperspectral and Spectral-Spatial Diffusion Geometry. arXiv preprint 
     arXiv: XXXXXXX.

Please write with any questions: samuel.polk@tufts.edu

(c) Copyright Sam L. Polk, Tufts University, 2021. 
