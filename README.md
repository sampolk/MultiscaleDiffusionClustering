# DiffusionLearning-INTEGRAL

This toolbox allows the implementation of the following diffusion-based clustering algorithms on very large datasets:

 - Learning by Unsupervised Nonlinear Diffusion (LUND)

Once this toolbox is downloaded, please do the following procedure. It only needs to be done once.

  1. Ensure that the AVIRS-NG India Forest dataset (titled "r1_reg.mat") is in your MATLAB path. 
  2. Run "preprocessing.m". This script will save nearest neighbors and a standardized version of the image locally. This file is 8.07 GB. 

After the above is done, you can analyze the AVIRIS-NG India Forest dataset using the main.m file. 

For large datasets like the AVIRS-NG India Forest dataset, we are often limited by RAM. On my MacBook Pro with 8GB of RAM, I can only store weight matrices with 1000 nearest neighbors in my workspace. 
