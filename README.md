# ccle_extractor
Extracting gene expression/metabolomics datasets from [CCLE](https://portals.broadinstitute.org/ccle) database

## Requirements
1. ```ccle_processing.py```
    > - [pandas](https://pandas.pydata.org)
    > - [matplotlib](https://matplotlib.org)
    > - [seaborn](https://seaborn.pydata.org)
1. ```counts_normalization.R```
    > - [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
    > - [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
    > - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

## CCLE Data
- CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
- Cell_lines_annotations_20181226.txt
- CCLE_RNAseq_genes_counts_20180929.gct.gz
- CCLE_metabolomics_20190502.csv

## Usage
1. Gene expression data
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

    # GeneCards link (https://www.genecards.org)
    selected_CCLE_subset.to_gene_summary()
    # TPM value
    selected_CCLE_subset.to_gene_expression()
    ```
1. Metabolomics dataset
    ```python
    from ccle_processing import CancerCellLineEncyclopedia as CCLE

    # Set metabolite_nemes
    # Set ccle_names or cell_lines

    selected_CCLE_subset = CCLE(
        metabolite_names = ['lactate', 'glutamine', 'glutamate', 'serine', 'alanine'],
        ccle_names = ['MCF7_BREAST', 'MDAMB231_BREAST']
    )

    # Metabolite levels
    selected_CCLE_subset.to_metabolomics()
    ```

## Installation
    $ git clone https://github.com/okadalabipr/ccle_extractor.git

## License
[MIT](LICENSE)