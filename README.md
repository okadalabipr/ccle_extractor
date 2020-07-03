# CCLEctor
Generating test and/or training datasets from [CCLE](https://portals.broadinstitute.org/ccle) data

## Requirements
> - [pandas](https://pandas.pydata.org)
> - [matplotlib](https://matplotlib.org)
> - [seaborn](https://seaborn.pydata.org)

## Default CCLE Data

If the file paths to your CCLE data is not specified, following CCLE data will be downloaded and used automaticallly.

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
    ccle_names = ['MCF7_BREAST', 'MDAMB231_BREAST'],
    expresssion = 'Path to your CCLE_RNAseq_rsem_genes_tpm_*.txt.gz',
    annotations = 'Path to your Cell_lines_annotations_*.txt',
    counts = 'Path to your CCLE_RNAseq_genes_counts_*.gct.gz'
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
selected_CCLE_subset.to_expression()
```

## Output
- gene_summary.md
- tpm_values.csv
- expression/

## Installation
    $ git clone https://github.com/okadalabipr/CCLEctor.git

## License
[MIT](LICENSE)