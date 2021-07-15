import os
import random
import shutil

from ccle.database import CancerCellLineEncyclopedia as CCLE

GENE_NAMES = random.sample(CCLE().counts_description, 10)
CCLE_NAMES = random.sample(list(CCLE().gene_expression_data.columns[1:]), 10)
COMPARISON = random.choice(["gene", "cell"])


def test_subset():
    for file in ["tpm_values.csv", "gene_summary.md"]:
        if os.path.isfile(file):
            os.remove(file)
    if os.path.isdir("gene_expression"):
        shutil.rmtree("gene_expression")
    selected_CCLE_subset = CCLE(
        gene_names=GENE_NAMES,
        ccle_names=CCLE_NAMES,
    )
    selected_CCLE_subset.to_gene_summary()
    assert os.path.isfile("gene_summary.md")
    os.remove("gene_summary.md")
    selected_CCLE_subset.to_gene_expression(comparison=COMPARISON)
    assert os.path.isfile("tpm_values.csv")
    assert os.path.isdir("gene_expression")
    if COMPARISON == "gene":
        for gene in GENE_NAMES:
            assert f"{gene}.pdf" in os.listdir("gene_expression")
    else:
        for ccle_id in CCLE_NAMES:
            cell = selected_CCLE_subset._id2cell(ccle_id)
            assert f"{cell}.pdf" in os.listdir("gene_expression")
    os.remove("tpm_values.csv")
    shutil.rmtree("gene_expression")


def test_all_counts():
    if os.path.isfile("CCLE_all_counts.csv"):
        os.remove("CCLE_all_counts.csv")
    CCLE().save_CCLE_all_counts_as_csv()
    assert os.path.isfile("CCLE_all_counts.csv")
    os.remove("CCLE_all_counts.csv")
