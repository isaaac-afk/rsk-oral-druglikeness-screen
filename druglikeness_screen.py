"""
Druglikeness Screen for Oral Kinase Inhibitors

Computes molecular descriptors for a curated set of FDA-approved
oral kinase inhibitors used in oncology, and flags Lipinski compliance.
"""

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

# Curated list: FDA-approved oral kinase inhibitors used in cancer therapy.
# SMILES sourced from PubChem / DrugBank.
DRUGS = [
    ("Imatinib",     "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("Gefitinib",    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
    ("Erlotinib",    "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"),
    ("Dasatinib",    "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1"),
    ("Sunitinib",    "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ccc(F)cc32)c1C"),
    ("Sorafenib",    "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1"),
    ("Lapatinib",    "CS(=O)(=O)CCNCc1oc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)cc1"),
    ("Palbociclib",  "CC(=O)c1c(C)c2cnc(Nc3ccc(N4CCNCC4)nc3)nc2n(C3CCCC3)c1=O"),
    ("Ribociclib",   "CN(C)C(=O)c1cc2cnc(Nc3ccc(N4CCNCC4)nc3)nc2n1C1CCCC1"),
    ("Osimertinib",  "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"),
    ("Ibrutinib",    "C=CC(=O)N1CCC[C@@H](n2nc(-c3ccc(Oc4ccccc4)cc3)c3c(N)ncnc32)C1"),
    ("Pazopanib",    "Cc1ccc(Nc2nccc(N(C)c3ccc4c(C)n(C)nc4c3)n2)cc1S(N)(=O)=O"),
    ("Regorafenib",  "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)c(F)c2)ccn1"),
    ("Axitinib",     "CNC(=O)c1ccccc1Sc1ccc2c(/C=C/c3ccccn3)n[nH]c2c1"),
    ("Nilotinib",    "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1"),
    ("Cabozantinib", "COc1cc2nccc(Oc3ccc(NC(=O)C4(C(=O)Nc5ccc(F)cc5)CC4)cc3F)c2cc1OC"),
    ("Crizotinib",   "C[C@H](Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl"),
    ("Vemurafenib",  "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(-c4ccc(Cl)cc4)cc23)c1F"),
    ("Bosutinib",    "COc1cc2ncc(C#N)c(Nc3cc(Cl)c(Cl)cc3OC)c2cc1OCCCN1CCN(C)CC1"),
    ("Abemaciclib",  "CCN1CCN(Cc2ccc(Nc3ncc(F)c(-c4cc5nc(C)n(C(C)C)c5cc4F)n3)nc2)CC1"),
]


def compute_descriptors(smiles: str):
    """Return Lipinski-relevant descriptors for a SMILES, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW":   Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD":  Descriptors.NumHDonors(mol),
        "HBA":  Descriptors.NumHAcceptors(mol),
    }


def lipinski_violations(desc: dict) -> int:
    """Count how many of the four Lipinski criteria a molecule violates."""
    v = 0
    if desc["MW"]   > 500: v += 1
    if desc["LogP"] > 5:   v += 1
    if desc["HBD"]  > 5:   v += 1
    if desc["HBA"]  > 10:  v += 1
    return v


def build_dataframe(drugs):
    rows = []
    for name, smi in drugs:
        d = compute_descriptors(smi)
        if d is None:
            print(f"  [WARN] Could not parse SMILES for {name}")
            continue
        d["Name"] = name
        d["SMILES"] = smi
        d["Violations"] = lipinski_violations(d)
        d["Passes_Lipinski"] = d["Violations"] == 0
        rows.append(d)
    df = pd.DataFrame(rows)
    return df[["Name", "MW", "LogP", "HBD", "HBA",
               "Violations", "Passes_Lipinski", "SMILES"]].round({"MW": 1, "LogP": 2})


def plot_mw_vs_logp(df, out_path="lipinski_plot.png"):
    fig, ax = plt.subplots(figsize=(11, 7))
    passing = df[df["Passes_Lipinski"]]
    failing = df[~df["Passes_Lipinski"]]

    ax.scatter(passing["MW"], passing["LogP"], c="#2E7D32", s=110,
               edgecolor="black", linewidth=0.6, zorder=3,
               label=f"Pass ({len(passing)})")
    ax.scatter(failing["MW"], failing["LogP"], c="#C62828", s=110,
               edgecolor="black", linewidth=0.6, zorder=3,
               label=f"Fail ({len(failing)})")

    ax.axvline(500, ls="--", color="gray", alpha=0.6)
    ax.axhline(5,   ls="--", color="gray", alpha=0.6)

    for _, r in df.iterrows():
        ax.annotate(r["Name"], (r["MW"], r["LogP"]),
                    fontsize=8, xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("Molecular Weight (Da)")
    ax.set_ylabel("LogP")
    ax.set_title("Lipinski's Rule of Five - Oral Kinase Inhibitors (FDA-Approved)")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"  Saved plot: {out_path}")


def main():
    print("Computing descriptors...")
    df = build_dataframe(DRUGS)
    df.to_csv("lipinski_results.csv", index=False)
    print(f"  Saved table: lipinski_results.csv")

    plot_mw_vs_logp(df)

    n_pass = df["Passes_Lipinski"].sum()
    print(f"\n{n_pass}/{len(df)} drugs pass strict Lipinski (0 violations)")
    print("\nFull table:")
    print(df.drop(columns="SMILES").to_string(index=False))


if __name__ == "__main__":
    main()