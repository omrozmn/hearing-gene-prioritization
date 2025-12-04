import pandas as pd
import os
import numpy as np

# Configuration
RESULTS_DIR = "results"
REPORT_PATH = os.path.join(RESULTS_DIR, "sanity_check_report.md")

def run_validation():
    print("--- PHASE 0: VALIDATION ---")
    
    report_lines = ["# Phase 0: Validation Report", "", "## 1. File Existence Check"]
    
    # 1. Check Files
    required_files = [
        "priority_score/gene_priority_scores.csv", # v2 (Corrected path)
        "TPS_v3_noise_integrated.csv",     # v3
        "TPS_v4_druggability.csv",         # v4
        "celltype_specific_scores.csv"
    ]
    
    all_exist = True
    for f in required_files:
        path = os.path.join(RESULTS_DIR, f)
        if os.path.exists(path):
            report_lines.append(f"- [x] Found {f}")
        else:
            report_lines.append(f"- [ ] MISSING {f}")
            all_exist = False
            
    if not all_exist:
        print("CRITICAL: Missing result files.")
        # return # Continue anyway to check what we have
        
    # 2. Load Final Scores (v4)
    print("Loading final scores...")
    try:
        df = pd.read_csv(os.path.join(RESULTS_DIR, "TPS_v4_druggability.csv"), index_col=0)
        report_lines.append("\n## 2. Score Distribution")
        report_lines.append(f"- Total Genes: {len(df)}")
        report_lines.append(f"- Score Range: {df['score_v4'].min():.4f} - {df['score_v4'].max():.4f}")
        
        if df['score_v4'].max() > 1.0001 or df['score_v4'].min() < -0.0001:
             report_lines.append("- [!] WARNING: Scores outside 0-1 range.")
        else:
             report_lines.append("- [x] Scores are properly normalized (0-1).")
             
        # 3. Positive Control Check
        report_lines.append("\n## 3. Positive Control Check (Top 100)")
        positive_controls = ['MYO7A', 'GJB2', 'POU4F3', 'TMC1', 'OTOF', 'SLC26A5']
        
        top100 = df.head(100).index.tolist()
        found_controls = []
        
        for gene in positive_controls:
            if gene in top100:
                rank = top100.index(gene) + 1
                found_controls.append(gene)
                report_lines.append(f"- [x] **{gene}** found in top 100 (Rank: {rank})")
            else:
                # Check if present in dataset but low rank
                if gene in df.index:
                    rank = df.index.get_loc(gene) + 1
                    report_lines.append(f"- [ ] {gene} present but low rank (Rank: {rank})")
                else:
                    report_lines.append(f"- [ ] {gene} NOT FOUND in dataset")
        
        if len(found_controls) >= 3:
             report_lines.append("\n**RESULT: PASS** (At least 3 positive controls in top 100)")
        else:
             report_lines.append("\n**RESULT: WARNING** (Fewer than 3 positive controls in top 100)")

    except Exception as e:
        report_lines.append(f"\n## Error during validation: {e}")

    # Write Report
    with open(REPORT_PATH, "w") as f:
        f.write("\n".join(report_lines))
        
    print(f"Validation complete. Report saved to {REPORT_PATH}")

if __name__ == "__main__":
    run_validation()
