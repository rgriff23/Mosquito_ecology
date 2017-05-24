# Mosquitoes of the Field and Forest: The scale of habitat segregation in a diverse mosquito assemblage

[Michael H. Reiskind](http://www.cals.ncsu.edu/entomology/reiskind), [Randi H. Griffin](http://rgriff23.github.io/), M. Shawn Janairo, and Kristen A. Hopperstad
___

This repository contains the data and code for [this study](http://onlinelibrary.wiley.com/doi/10.1111/mve.12193/full) on the landscape ecology of mosquito communities in North Carolina, which was published in the journal *Medical and Veternary Entomology*. I also published a `vegan` tutorial that uses the data from this study on my [blog](https://rgriff23.github.io/2017/05/23/mosquito-community-ecology-in-vegan.html).

## Analysis

The **Analysis** folder contains all the data and code needed to reproduce our results. Several R packages must be installed for the code to work:

```
install.packages("vegan")
install.packages("plyr")
install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", 
getOption("repos")), type="source")
```

If there are problems installing glmmADMB, see installation tips on the [glmmADMB package website](http://glmmadmb.r-forge.r-project.org/).

### Data

`data.csv` contains landscape classifications and abundances of different mosquito species for all of the trap sites in our study. 

### R code

`analysis.R` contains all of the R code to run the analyses and make the figures in our study. 

--- 
