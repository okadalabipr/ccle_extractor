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
    def __init__(
            self,
            gene_names=[],
            metabolite_names=[],
            cell_lines=[],
            ccle_names=[],
    ):
        self.gene_names = gene_names
        self.metabolite_names = metabolite_names
        self.cell_lines = cell_lines
        self.ccle_names = ccle_names
        
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
            header=2
        )
        self.metabolomics_data = pd.read_csv(
            'https://data.broadinstitute.org/ccle/'\
            'CCLE_metabolomics_20190502.csv',
            index_col=0
        )

        self.annotations_name = set(self.annotations.Name)
        self.annotations_ccle_id = set(self.annotations.CCLE_ID)
        self.counts_description = set(self.counts.Description)
        self.all_metabolites = set(self.metabolomics_data.columns[1:])

    def _gene2id(self, gene):
        if gene not in self.counts_description:
            print("gene '{}' does not exist.\n".format(gene))
            return False
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
        if cell not in self.annotations_name:
            raise ValueError(cell)
        else:
            ccle_id = self.annotations.at[
                list(self.annotations.Name).index(cell), 'CCLE_ID'
            ]
            return ccle_id

    def _id2cell(self, ccle_id):
        if ccle_id not in self.annotations_ccle_id:
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
            a_gene_id = self._gene2id(gene)
            if not a_gene_id:
                pass
            else:
                gene_ids.append(a_gene_id)
        return gene_ids
    
    def _get_ccle_id(self):
        ccle_ids = []
        for cell in self.cell_lines:
            ccle_ids.append(self._cell2id(cell))
        return ccle_ids

    def _extract_tpm(self):
        gene_ids = self._get_gene_id()
        ccle_ids = self._get_ccle_id() \
            if len(self.ccle_names) < len(self.cell_lines) else self.ccle_names
        data = self.gene_expression_data.loc[gene_ids, ccle_ids]
        data.index.name = None
        data.rename(
            index=lambda x: self._id2gene(x),
            columns=lambda x: self._id2cell(x),
            inplace=True
        )
        data.to_csv('tpm_values.csv')
        return data
    
    def _extract_metabolite_level(self):
        not_exist = []
        for i, metabolite in enumerate(self.metabolite_names):
            if metabolite not in self.all_metabolites:
                print("metabolite '{}' does not exist.\n".format(metabolite))
                not_exist.append(i)
        for i in sorted(not_exist, reverse=True):
            self.metabolite_names.pop(i)
        ccle_ids = self._get_ccle_id() \
            if len(self.ccle_names) < len(self.cell_lines) else self.ccle_names
        data = self.metabolomics_data.loc[ccle_ids, self.metabolite_names]
        data.index.name = None
        data.rename(
            index=lambda x: self._id2cell(x),
            inplace=True
        )
        data.to_csv('metabolite_levels.csv')
        return data
    
    def save_CCLE_all_counts_as_csv(self):
        """
        Example
        -------
        >>> from ccle_processing import CancerCellLineEncyclopedia as CCLE
        >>> CCLE().save_CCLE_all_counts_as_csv()
        
        """
        self.counts.to_csv('CCLE_all_counts.csv')

    def to_gene_expression(self):
        if not self.cell_lines and not self.ccle_names:
            raise ValueError('cell_lines or ccle_names must be filled in.')
        os.makedirs('./gene_expression', exist_ok=True)
        data = self._extract_tpm()
        self._set_rcParams()
        for gene in self.gene_names:
            if gene in self.counts_description:
                ax = data.loc[gene].plot.bar(
                    figsize=(
                        2*max(len(self.cell_lines), len(self.ccle_names)), 6
                    ), fontsize=28, title=gene,  # r'$\it{'+gene+'}$'
                )
                ax.set_yscale('log', basey=2)
                ax.set_ylabel('log2 (TPM+1)')
                sns.despine()
                plt.savefig(
                    './gene_expression/{}.pdf'.format(gene), bbox_inches='tight'
                )
                plt.close()

        ccle_ids = self._get_ccle_id() \
            if len(self.ccle_names) < len(self.cell_lines) else self.ccle_names
        for cid in ccle_ids:
            cell = self._id2cell(cid)
            for gene in self.gene_names:
                data.loc[gene, cell] += 1
            ax = data.loc[:, cell].plot.bar(
                figsize=(
                    len(self.gene_names), 6
                ), fontsize=28, #title=cell
            )
            ax.set_yscale('log', basey=2)
            ax.set_ylabel('log2 (TPM+1)')
            sns.despine()
            plt.savefig(
                './gene_expression/{}.pdf'.format(cell), bbox_inches='tight'
            )
            plt.close()
    
    def to_gene_summary(self):
        with open('gene_summary.md', mode='w') as f:
            f.write('|gene_name|gene_id|GeneCards_URL|\n'\
                    '|---------|-------|-------------|\n')
            for gene in self.gene_names:
                if gene in self.counts_description:
                    gene_id = self._gene2id(gene)
                    gene_cards_url = (
                        'https://www.genecards.org/'\
                        'cgi-bin/carddisp.pl?gene='+gene
                    )
                    f.write(
                        '|'+gene+'|'+gene_id+'|'+gene_cards_url+'|\n'
                    )

    def to_metabolomics(self):
        if not self.cell_lines and not self.ccle_names:
            raise ValueError('cell_lines or ccle_names must be filled in.')
        os.makedirs('./metabolomics', exist_ok=True)
        data = self._extract_metabolite_level()
        self._set_rcParams()
        for metabolite in self.metabolite_names:
            ax = data.loc[:, metabolite].plot.bar(
                figsize=(
                    2*max(len(self.cell_lines), len(self.ccle_names)), 6
                ), fontsize=28, title=metabolite,
            )
            ax.set_ylabel('Metabolite levels')
            sns.despine()
            plt.savefig(
                './metabolomics/{}.pdf'.format(metabolite), bbox_inches='tight'
            )
            plt.close()

    @staticmethod
    def _set_rcParams():
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