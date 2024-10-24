# GIFT

GIFT(Gene-based Integrative Fine-mapping through conditional TWAS), is an R package for efficient statistical inference of conditional TWAS fine-mapping. GIFT examines one genomic region at a time, jointly models the GReX of all genes residing in the focal region, and carries out TWAS conditional analysis in a maximum likelihood framework. In the process, GIFT explicitly models the gene expression correlation and cis-SNP LD across different genes in the region, accounts for the uncertainty in the constructed GReX through joint inference and provides calibrated p values.

## Installation
It is easy to install the development version of GIFT package using the 'devtools' package. 

```
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/GIFT")
```
## Usage
The main function in the package is GIFT, you can find the instructions by '?GIFT'.
```
library(GIFT)

?GIFT
```

## Quick Start

See [Tutorial](https://yuanzhongshang.github.io/GIFT/) for detailed documentation and examples.

## Version History

| Version | Updated      | Description                                                                                                                                                                                                                                             |
|---------|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [1.0](https://github.com/yuanzhongshang/GIFT/tree/v1.0) | 2022-10-31| Initial release.                                                                                                                                                                                                                   |
| [2.0](https://github.com/yuanzhongshang/GIFT/tree/v2.0) | 2024-10-17|Optimized some code and added a filter parameter to improve running speed.                                                                                                                                 |
| [2.1](https://github.com/yuanzhongshang/GIFT/tree/v2.1) | 2024-10-24|Optimized some code and added a split parameter to improve running speed.                                                                                                                                       |

## Development
This R package is developed by Lu Liu, Zhongshang Yuan and Xiang Zhou.
