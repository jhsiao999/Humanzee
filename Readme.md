### Humanzee: Gilad lab's tooks for comparative analysis of genomics data

[Gilad lab at the University of Chicago](http://giladlab.uchicago.edu/)

[Humanzee](https://github.com/jhsiao999/Humanzee) is an [R](http://www.r-project.org) package of tools for
statistical analysis, data visualization, and data preprocessing for the projects in the Gilad lab at 
the [Universit of Chicago](https://www.uchicago.edu)
interactive charts for QTL data. \[[website](http://)\]

For example usage in our projects, see the [Humanzee website](http://webaddress).

#### Installation

Install Humanzee from its
[GitHub repository](https://github.com/jhsiao999/Humanzee). You first need to
install [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install Humanzee using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package. 


```r
library(devtools)
install_github("jhsiao999/Humanzee")
```

#### Licenses

Code is released under the GPLv3 license.

Text is released under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.

Humanzee adapts some of the limma functions for linear models of different covariates across genes (lmFit) and
linaer models of different random effects across genes (mixedModel2Fit), etc. 


#### Suggestions and feedbacks

Share your Humanzee experiencs with us. The package is maintained by Joyce Hsiao (joyce.hsiao1@gmail.com).
