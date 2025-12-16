# HPC

This folder contains HPC (High-Performance Computing) specific scripts and configuration files.
These resources are designed to facilitate the calculation and management of large-scale computational tasks such as
calculating some properties of the compounds and reactions.

Submission scripts use SLURM and start with `sub_`. The scripts are tailored for SOL HPC at ASU. They can be modified to
suit other HPC environments as needed.

Files that start with `calc_` are used to perform computational chemistry calculations of the compounds on the HPC cluster.

Files that start with `reaction_` are used to perform calculations related to chemical reactions on the HPC cluster.

To use these scripts, ensure you have access to the HPC cluster and the necessary permissions to submit jobs.
Modify the scripts as needed to fit your specific computational requirements and HPC environment.