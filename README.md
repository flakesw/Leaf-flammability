This repository contains data and analysis scripts for the forthcoming manuscript, Linking leaf traits and litter flammability using a novel framework, tested with Brazilian Cerrado trees, by Samuel Flake, Patrick Elliott, Giselda Durigan, Davi Rossatto, and William Hoffmann. It contains all of the data necessary to generate figures 3-8 in that manuscript. The analysis was performed by Sam Flake -- please contact him at swflake@ncsu.edu with questions.

The repository contains the following files:

analysis_2025-2-18.R -- this script performs all the analysis and generates most of the figures in the paper. It uses data stored in the ./clean_data folder and outputs the tables and figures found in ./Outputs

data_prep.R -- this script takes the raw data files and creates the clean data files. It's pretty run-of-the-mill data wrangling to fix units, create some new variables, and combine and aggregate data to analysis-ready tables.

process_point_clouds.R -- this script was used to extract convex hull volumes from point clouds from 3D-scanned leaves. The original point cloud data is too large to host in this repository, but Sam is happy to share it upon request -- please use the data if you're interested in it!

./raw_data -- raw data from leaf measurements and burn experiments
./clean_data -- analysis-ready data
./Outputs -- figures and tabular outputs from analysis
