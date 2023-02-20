# POR_dinoRNAV_MVP

This repository contains the code associated with analyses of the dinoRNAV and Symbiodiniaceae community composition of Porites lobata coral holobionts in Moorea, French Polynesia. Analyses are associated with the manucript "Viruses of a key coral symbiont exhibit temperature-driven productivity across a reefscape."

#Symiodiniaceae LSU data:

LSU amplicon sequencing of Porites DNA was used to identify the Symbiodiniaceae present in each colony. Sequences were processed through the DADA2 pipeline (Callahan et al 2016) and then ASVs were collapsed using the LULU pipeline (Frøslev et al., 2017). Files and code associated with the LULU ASV curation are uploaded here.

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. https://doi.org/10.1038/nmeth.3869

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1). https ://doi.org/10.1038/s41467-017-01312-x

#dinoRNAV analyses:

All non-read files needed to reproduce read processing with vAMPirus are available here: https://zenodo.org/record/7552892#.Y_PAa-zMK3I

R scripts associated with diversity analyses and creating the phylogenetic tree have been uploaded here.

All mcp gene amplicon libraries and Symbiodiniaceae LSU gene amplicon libraries are available at the Sequence Read Archive under accession numbers PRJNA928208 and PRJNA930706, respectively.
