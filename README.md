# ccle_processing
Generating test and/or training datasets from [CCLE](https://portals.broadinstitute.org/ccle) data

## Requirements
> - [pandas](https://pandas.pydata.org)
> - [matplotlib](https://matplotlib.org)
> - [seaborn](https://seaborn.pydata.org)

## CCLE Data
- CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
- Cell_lines_annotations_20181226.txt
- CCLE_RNAseq_genes_counts_20180929.gct.gz

## Usage
```python
from ccle_processing import CancerCellLineEncyclopedia as CCLE

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

selected_CCLE_subset.to_expression()
```

## Output
- tpm_values.csv
- expression/

## Installation
    $ git clone https://github.com/okadalabipr/ccle_processing.git

## License
[MIT](LICENSE)