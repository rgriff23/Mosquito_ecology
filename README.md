# Mosquitoes of the Field and Forest: The scale of habitat segregation in a diverse mosquito assemblage

[Michael H. Reiskind](http://www.cals.ncsu.edu/entomology/reiskind), [Randi H. Griffin](), M. Shawn Janairo, and Kristen A. Hopperstad

* MHR designed the study and wrote the manuscript. RHG and MSJ collected and identified mosquitoes. KAH collected landscape data and made the landscape figure. RHG analyzed the data, made data figures, and wrote the statistical methods and results portions of the manuscript. 

___

This repository contains the data and code for our study on the landscape ecology of mosquito communities in North Carolina. I (Randi Griffin) also provide a synopsis of the study. A manuscript is in prep for publication as of August 2015.

## Manuscript

The **Manuscript** folder contains a synopsis of our study in the file `Synopsis.md`. 

The **Manuscript/figures** folder contains the figures that are included in the synopsis.

## Analysis

The **Analysis** folder contains all the data and code needed to reproduce our results. Several R packages must be installed for the code to work:

```
install.packages("vegan")
install.packages("plyr")
install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
```

Note that `glmmADMB` is not available on CRAN. The code above worked for me, but if there are problems installing glmmADMB, see installation tips on the [glmmADMB package website](http://glmmadmb.r-forge.r-project.org/).

### Data

`FinalData.csv` contains landscape classifications and abundances of different mosquito species for all of the trap sites in our study. 

### R code

`FinalAnalysis.R` contains all of the R code to run the analyses and make the figures in our study. 

--- 
