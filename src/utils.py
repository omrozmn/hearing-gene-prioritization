import os
import pandas as pd
import matplotlib.pyplot as plt

def ensure_dir(path):
    """Ensure directory exists."""
    if not os.path.exists(path):
        os.makedirs(path)

def load_gene_list(path):
    """Load a list of genes from a text file."""
    if not os.path.exists(path):
        print(f"Warning: Gene list {path} not found.")
        return []
    try:
        df = pd.read_csv(path, header=None)
        return df[0].astype(str).str.upper().tolist()
    except Exception as e:
        print(f"Error loading gene list: {e}")
        return []

def save_fig(fig, filename):
    """Save matplotlib figure to results/figures."""
    ensure_dir(os.path.dirname(filename))
    fig.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close(fig)

def load_gwas_genes(path):
    """Load GWAS gene list and return as a set of uppercase gene names."""
    if not os.path.exists(path):
        print(f"Warning: GWAS gene list {path} not found. Returning empty set.")
        return set()
    try:
        df = pd.read_csv(path, header=None)
        return set(df[0].astype(str).str.upper())
    except Exception as e:
        print(f"Error loading GWAS genes: {e}")
        return set()

def load_conservation(path):
    """Load conservation scores and return as dict {GENE: score}."""
    if not os.path.exists(path):
        print(f"Warning: Conservation file {path} not found. Returning empty dict.")
        return {}
    try:
        df = pd.read_csv(path, sep='\t')
        return dict(zip(df['GENE'].str.upper(), df['CONSERVATION']))
    except Exception as e:
        print(f"Error loading conservation scores: {e}")
        return {}
