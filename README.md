[![Build Status](https://travis-ci.com/xmc811/Scillus.svg?branch=master)](https://travis-ci.com/xmc811/Scillus)
[![Build status](https://ci.appveyor.com/api/projects/status/dkq1xn6574kqgs0s/branch/master?svg=true)](https://ci.appveyor.com/project/xmc811/scillus/branch/master)

<img align="right" width="108" height="125" src="Scillus.png">


Please visit the Scillus documentation website:

[scillus.netlify.com](http://scillus.netlify.com)


[![Netlify Status](https://api.netlify.com/api/v1/badges/eadbcb9a-16d1-4a9a-9e50-c0e8d4104ddc/deploy-status)](https://app.netlify.com/sites/scillus/deploys)

# Scillus v0.4.0

Please use the following code to install and load the package:

```R
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)
```

### Version History

September 23, 2020

* Version 0.4.0
Additonal improvements for R Shiny applications


January 29, 2020

* Version 0.3.0


December 20, 2019

* Version 0.2.0
Optimized processing and plotting functions; simplified installation process


August 25, 2019

* Version 0.1.0
Initial release; essential functions for streamlined scRNA-seq analysis and visualization
