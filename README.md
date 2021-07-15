# ccle_extractor [![Actions Status](https://github.com/okadalabipr/ccle_extractor/workflows/Tests/badge.svg)](https://github.com/okadalabipr/ccle_extractor/actions)
Extracting gene expression datasets from [CCLE](https://portals.broadinstitute.org/ccle) database.

## Manual installation of package requirements

| Language | Packages |
| ---      | ---      |
| Python >= 3.7   | [pandas](https://pandas.pydata.org), [matplotlib](https://matplotlib.org), [seaborn](https://seaborn.pydata.org) |
| R | [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html), [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) |

## CCLE Data

> - CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
> - Cell_lines_annotations_20181226.txt
> - CCLE_RNAseq_genes_counts_20180929.gct.gz

## Usage

```python
from ccle.database import CancerCellLineEncyclopedia as CCLE

# Set gene_nemes
# Set ccle_names or cell_lines

selected_CCLE_subset = CCLE(
    gene_names = ['EGFR', 'ERBB2', 'ERBB3', 'ERBB4'],
    ccle_names = ['MCF7_BREAST', 'MDAMB231_BREAST']
)

''' or
selected_CCLE_subset = CCLE(
    gene_names = ['EGFR', 'ERBB2', 'ERBB3', 'ERBB4'],
    cell_lines = ['MCF7', 'MDA-MB-231']
)
'''

# GeneCards link (https://www.genecards.org)
selected_CCLE_subset.to_gene_summary()
# TPM value
selected_CCLE_subset.to_gene_expression()
```


## Installation

    $ git clone https://github.com/okadalabipr/ccle_extractor.git

## License

[MIT](LICENSE)