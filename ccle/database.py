import math
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


@dataclass
class CancerCellLineEncyclopedia(object):
    """
    Extract gene expression data from https://portals.broadinstitute.org/ccle.
    """

    gene_names: List[str] = field(default_factory=list)
    cell_lines: List[str] = field(default_factory=list)
    ccle_names: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:

        self.gene_expression_data = pd.read_table(
            "https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",
            index_col=0,
        )
        self.annotations = pd.read_table(
            "https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt"
        )
        self.counts = pd.read_table(
            "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz",
            header=2,
        )

        self.annotations_name = set(self.annotations.Name)
        self.annotations_ccle_id = set(self.annotations.CCLE_ID)
        self.counts_description = set(self.counts.Description)

    def _gene2id(self, gene):
        if gene not in self.counts_description:
            print("gene '{}' does not exist.\n".format(gene))
            return False
        else:
            gene_id = self.counts.at[list(self.counts.Description).index(gene), "Name"]
            return gene_id

    def _id2gene(self, gene_id):
        gene = self.counts.at[list(self.counts.Name).index(gene_id), "Description"]
        return gene

    def _cell2id(self, cell):
        if cell not in self.annotations_name:
            raise ValueError(cell)
        else:
            ccle_id = self.annotations.at[list(self.annotations.Name).index(cell), "CCLE_ID"]
            return ccle_id

    def _id2cell(self, ccle_id):
        if ccle_id not in self.annotations_ccle_id:
            raise ValueError(ccle_id)
        else:
            cell = self.annotations.at[list(self.annotations.CCLE_ID).index(ccle_id), "Name"]
            try:
                if math.isnan(float(cell)):
                    cell = re.findall("(.*?)_", ccle_id)[0]
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
        ccle_ids = (
            self._get_ccle_id() if len(self.ccle_names) < len(self.cell_lines) else self.ccle_names
        )
        data = self.gene_expression_data.loc[gene_ids, ccle_ids]
        data.index.name = None
        data.rename(
            index=lambda x: self._id2gene(x), columns=lambda x: self._id2cell(x), inplace=True
        )
        data.to_csv("tpm_values.csv")
        return data

    def save_CCLE_all_counts_as_csv(self):
        """
        Example
        -------
        >>> from ccle.database import CancerCellLineEncyclopedia as CCLE
        >>> CCLE().save_CCLE_all_counts_as_csv()

        """
        self.counts.to_csv("CCLE_all_counts.csv")

    def to_gene_expression(
        self,
        *,
        config: Optional[dict] = None,
        comparison: Optional[str] = None,
    ):
        if not self.cell_lines and not self.ccle_names:
            raise ValueError("cell_lines or ccle_names must be filled in.")
        os.makedirs("gene_expression", exist_ok=True)
        data = self._extract_tpm()
        if config is None:
            config = {}
        config.setdefault("font.size", 28)
        config.setdefault("axes.linewidth", 2)
        config.setdefault("xtick.major.width", 2)
        config.setdefault("ytick.major.width", 2)
        config.setdefault("savefig.bbox", "tight")
        plt.rcParams.update(config)
        if comparison is None:
            pass
        elif comparison == "gene":
            for gene in self.gene_names:
                if gene in self.counts_description:
                    ax = data.loc[gene].plot.bar(
                        figsize=(2 * max(len(self.cell_lines), len(self.ccle_names)), 6),
                        fontsize=28,
                        title=gene,  # r'$\it{'+gene+'}$'
                    )
                    ax.set_yscale("log", base=2)
                    ax.set_ylabel("log2 (TPM+1)")
                    sns.despine()
                    plt.savefig(os.path.join("gene_expression", f"{gene}.pdf"))
                    plt.close()
        elif comparison == "cell":
            ccle_ids = (
                self._get_ccle_id()
                if len(self.ccle_names) < len(self.cell_lines)
                else self.ccle_names
            )
            for cid in ccle_ids:
                cell = self._id2cell(cid)
                for gene in self.gene_names:
                    data.loc[gene, cell] += 1
                ax = data.loc[:, cell].plot.bar(
                    figsize=(len(self.gene_names), 6),
                    fontsize=28,  # title=cell
                )
                ax.set_yscale("log", base=2)
                ax.set_ylabel("log2 (TPM+1)")
                sns.despine()
                plt.savefig(os.path.join("gene_expression", f"{cell}.pdf"))
                plt.close()
        else:
            raise ValueError("`comparison` must be either 'gene' or 'cell'.")

    def to_gene_summary(self):
        """
        Create table displaying gene id and annotation.
        """
        with open("gene_summary.md", mode="w") as f:
            f.write("|gene_name|gene_id|GeneCards_URL|\n|---------|-------|-------------|\n")
            for gene in self.gene_names:
                if gene in self.counts_description:
                    gene_id = self._gene2id(gene)
                    gene_cards_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene
                    f.write(f"| {gene} | {gene_id} | {gene_cards_url} |\n")
