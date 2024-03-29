# MultiscaleDiffusionClustering Toolbox

This toolbox allows the implementation of the following diffusion-based clustering algorithms on synthetic and real datasets included in the repository:

   - Learning by Unsupervised Nonlinear Diffusion (LUND)
   - Multiscale Learning by Unsupervised Nonlinear Diffusion (M-LUND)
   - Spatially Regularized Diffusion Learning (SRDL)
   - Multiscale Spatially Regularized Diffusion Learning (M-SRDL)

This package can be used to generate experiments in the following articles:

   - Murphy, James M., & Polk, Sam L. (2022). A multiscale environment for learning by diffusion. Applied and Computational Harmonic Analysis, 57, 58-100.
   - Polk, Sam L. and Murphy James M. "Multiscale Spectral-Spatial Diffusion Geometry for Hyperspectral Image Clustering." To Appear In The Proceedings of IEEE IGARSS 2021 (2021).

The following scripts (in the Experiments folder) generate the relevant experiments:

   - SADataVisualization.m implements visualizes Salinas A hyperspectral image.
   - M_LUND_demo.m implements M_LUND on synthetic data and the Salinas A hyperspectral image. One can compare against related algorithms on Salinas A using this demo as well.
   - Benchmark.m compares the M-LUND algorithm against related algorithms on eleven benchmark datasets.
   - M_SRDL_demo.m compares M-SRDL against M-LUND on the Salinas A hyperspectral image.

All necessary datasets are contained in this repository, so no additional data downloads are necessary. All benchmark datasets were obtained from the UCI Machine Learning Repository:

   - https://archive.ics.uci.edu/ml/index.php

The real hyperspectral image data (Salinas A) was from the 2000 IEEE Data Fusion Contest data. This data is publically available via IEEE:

   - http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes

If MMS clustering toolboxes are not in one's path, comparisons against MMS clustering will automatically be skipped. To run comparisons against MMS clustering, the following toolboxes, written by Zhijing Liu and Mauricio Barahona, must be added to one's path:

   - https://github.com/barahona-research-group/GraphBasedClustering
   - https://wwwf.imperial.ac.uk/~mpbara/Partition_Stability/

Users are free to modify the Multiscale Diffusion Clustering toolbox as they wish. If you find it useful or use it in any publications, please cite the following papers:

   - Maggioni, Mauro, and James M. Murphy. "Learning by Unsupervised Nonlinear Diffusion." Journal of Machine Learning Research 20.160 (2019): 1-56.
   - Murphy, James M., and Mauro Maggioni. "Spectral-spatial diffusion geometry for hyperspectral image clustering." IEEE Geoscience and Remote Sensing Letters (2019).
   - Murphy, James M., & Polk, Sam L. (2022). A multiscale environment for learning by diffusion. Applied and Computational Harmonic Analysis, 57, 58-100.
   - Polk, Sam L. and Murphy James M. "Multiscale Spectral-Spatial Diffusion Geometry for Hyperspectral Image Clustering." Proceedings of IEEE IGARSS 2021 (2021): 4688-4691.

Please write with any questions: samuel.polk@tufts.edu

(c) Copyright Sam L. Polk, Tufts University, 2022.
