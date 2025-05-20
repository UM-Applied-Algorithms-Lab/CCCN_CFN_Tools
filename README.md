This is the supporting package for the projects in the papers listed below.

Please feel free to reach out with errors, critiques, or questions.

## Installation

Before proceeding, ensure you have R installed. You will also need the devtools package, which can be installed with:

```
install.packages("devtools")
```

TODO: understand STRINGdb and figure out how to improve this

Additionally, please check that you have R version 4.5, and STRINGdb installed before moving forward.

One option for installing stringdb is using conda:

```
conda install bioconda::bioconductor-stringdb
```

Once these are installed, you can install this package directly from GitHub using:

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

```
git clone https://github.com/UM-Applied-Algorithms-Lab/CCCN_CFN_Tools.git
```

Then, in R, navigate to the package directory. You will need to use devtools to load the package for development, so ensure you have that installed and libraried:

```
install.packages("devtools")
library(devtools)
```

You may then load the package for development:

```
load_all()
```

To check for issues before committing changes, run:

```
check()
```
