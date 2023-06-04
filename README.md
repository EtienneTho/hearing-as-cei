# hearing-as-cei
This repository provides the Matlab scripts useful to generate the analyses and the figures of the article: 
* Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation. Communications Biology

# Files description
Each Matlab file corresponds to one figure of sub-figure from the paper. Each script should generates figures (.eps + .fig) and source data files (.csv) in the same folder.

# Dependences
These scripts are using libraries embedded in the subfoldeer ./lib/ in particular the Empirical Modes Decomposition library developed by Patrick Flandrin and Gabriel Rilling: https://perso.ens-lyon.fr/patrick.flandrin/emd.html

Please cite us and them if you use part of these scripts :
* Flandrin, P., Rilling, G., & Goncalves, P. (2004). Empirical mode decomposition as a filter bank. IEEE signal processing letters, 11(2), 112-114.
* Rilling, G., & Flandrin, P. (2007). One or two frequencies? The empirical mode decomposition answers. IEEE transactions on signal processing, 56(1), 85-95.

Other functions picked up on different repositories are used:
* 'btqn.m': http://personal.psu.edu/drh20/code/btmatlab/

# Datasets
* The Making Sense of sounds dataset can be accessed here: https://doi.org/10.17866/rd.salford.6901475.v4
* Please cite them if re-used: Harris, Lara; Bones, Oliver Charles (2018). Making Sense Of Sounds: Data for the machine learning challenge 2018. University of Salford. Dataset. https://doi.org/10.17866/rd.salford.6901475.v4
* Please email me if you need the exact folder I reorganized, basically, all the excerpts from the development set have been pasted in the same subfolder './cmos/'

## Organisation
Each file is named with the name of a figure or a subfigure of the paper. When run, each script should generate a .eps a .fig and a .csv file.

Please let us know any bugs or questions.

Contact:
etiennethoret [at] gmail [dot] com
