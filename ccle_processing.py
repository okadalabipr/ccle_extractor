import os
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler


class CancerCellLineEncyclopedia(object):
    """https://portals.broadinstitute.org/ccle
    """
    def __init__(self, gene_names, cell_lines=[], ccle_names=[]):
        self.gene_names = gene_names
        self.cell_lines = cell_lines
        self.ccle_names = ccle_names

        if not cell_lines and not ccle_names:
            raise ValueError(
                'cell_lines or ccle_names must be filled in.'
            )
        
        self.gene_expression_data = pd.read_table(
            'https://data.broadinstitute.org/ccle/'\
            'CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz',
            index_col=0
        )
        self.annotations = pd.read_table(
            'https://data.broadinstitute.org/ccle/'\
            'Cell_lines_annotations_20181226.txt'
        )
        self.counts = pd.read_table(
            'https://data.broadinstitute.org/ccle/'\
            'CCLE_RNAseq_genes_counts_20180929.gct.gz',
            header=2, usecols=range(2)
        )

    def _gene2id(self, gene):
        if gene not in set(self.counts.Description):
            raise ValueError(gene)
        else:
            gene_id = self.counts.at[
                list(self.counts.Description).index(gene), 'Name'
            ]
            return gene_id
    
    def _id2gene(self, gene_id):
        gene = self.counts.at[
            list(self.counts.Name).index(gene_id), 'Description'
        ]
        return gene

    def _cell2id(self, cell):
        if cell not in set(self.annotations.Name):
            raise ValueError(cell)
        else:
            ccle_id = self.annotations.at[
                list(self.annotations.Name).index(cell), 'CCLE_ID'
            ]
            return ccle_id

    def _id2cell(self, ccle_id):
        if ccle_id not in set(self.annotations.CCLE_ID):
            raise ValueError(ccle_id)
        else:
            cell = self.annotations.at[
                list(self.annotations.CCLE_ID).index(ccle_id), 'Name'
            ]
            try:
                if math.isnan(float(cell)):
                    cell = re.findall('(.*?)_', ccle_id)[0]
            except ValueError:
                pass
            return cell

    def _get_gene_id(self):
        gene_ids = []
        for gene in self.gene_names:
            gene_ids.append(self._gene2id(gene))
        return gene_ids
    
    def _get_ccle_id(self):
        ccle_ids = []
        for cell in self.cell_lines:
            ccle_ids.append(self._cell2id(cell))
        return ccle_ids

    def extract(self):
        gene_ids = self._get_gene_id()
        ccle_ids = self._get_ccle_id() \
            if len(self.ccle_names) < len(self.cell_lines) else self.ccle_names
        data = self.gene_expression_data.loc[gene_ids, ccle_ids]
        del data.index.name
        data.rename(
            index=lambda x: self._id2gene(x),
            columns=lambda x: self._id2cell(x),
            inplace=True
        )
        data.to_csv('tpm_values.csv')
        return data

    def to_expression(self):
        os.makedirs('./expression', exist_ok=True)
        data = self.extract()
        # rcParams
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 28
        plt.rcParams['axes.linewidth'] = 2
        plt.rcParams['xtick.major.width'] = 2
        plt.rcParams['ytick.major.width'] = 2
        plt.rcParams['axes.prop_cycle'] = cycler(
            color=[
                '#4E79A7', '#F28E2B', '#E15759', '#76B7B2','#59A14E',
                '#EDC949','#B07AA2','#FF9DA7','#9C755F','#BAB0AC'
            ]
        )
        # ---
        for gene in self.gene_names:
            ax = data.loc[gene].plot.bar(
                figsize=(2*max(len(self.cell_lines), len(self.ccle_names)), 6),
                fontsize=28,
                title=gene,  # r'$\it{'+gene+'}$'
            )
            ax.set_ylabel('TPM')
            sns.despine()
            plt.savefig(
                './expression/{}.pdf'.format(gene), bbox_inches='tight'
            )
            plt.close()
