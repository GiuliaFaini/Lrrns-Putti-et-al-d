# Lrrns-Putti-et-al-

# Calcium imaging analysis script to separate responses to different visual stimulations in zebrafish

This repository contains a MATLAB script to analyze calcium imaging data from zebrafish and to separate neural responses to different visual stimulations.
This is related to the paper : "Lrrns define a visual circuit underlying brightness and contrast perception", Putti et al. 
https://pubmed.ncbi.nlm.nih.gov/40654605/

The script takes a mean calcium trace as input, applies OASIS-based deconvolution, and reconstructs response components restricted to predefined stimulation intervals.

For each visual stimulus interval, the script extracts peak amplitudes and integrated responses, and generates summary plots showing:
- the mean calcium trace with shaded stimulation intervals  
- the reconstructed calcium components corresponding to stimulus and non-stimulus periods  

The code is lightweight, well-documented, and designed to be easily adapted to different visual stimulation protocols by modifying the interval definitions.
