This is the supporting package for the projects in the papers listed below.

Please feel free to reach out with errors, critiques, or questions.

## Installation

Before proceeding, ensure you have R and RStudio installed. You will also need the devtools package, which can be installed with:

TODO: can we remove language around R studio? Does it really have to be R studio?

```
install.packages("devtools")
```

Additionally, please check that you have R version 4.4, and STRINGdb installed before moving forward.

One option for installing stringdb is using conda:

```
conda install bioconda::bioconductor-stringdb
```

Once these are both installed, you can install this package directly from GitHub using:

```
devtools::install_github("UM-Applied-Algorithms-Lab/CCCN_CFN_Tools")
```

## Usage

After installation, load the package in R with:

```
library(cccn.cfn.tools)
```

You can then use the available functions as described in the package documentation.

## Development & Contribution

If you wish to modify or contribute to the package, clone the repository locally using:

git clone https://github.com/KaylBing/CCCN_CFN_Tools.git

Then, in R, navigate to the package directory and use devtools to load the package for development:

devtools::load_all()

To check for issues before committing changes, run:

devtools::check()

