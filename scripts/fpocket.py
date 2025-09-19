import os
import csv
import re

# --------- USER CONFIGURATION ----------
fpocket_output_base = "./PDB"  # Main folder containing *_out subfolders
output_csv = "fpocket_detailed_summary.csv"
# --------------------------------------

def extract_fpocket_info(info_file):
    data = {
        "num_pockets": 0,
        "top_pocket": {
            "druggability_score": 0.0,
            "volume": 0.0,
            "hydrophobicity_score": 0.0,
            "polarity_score": 0.0,
            "total_sasa": 0.0
        }
    }

    try:
        with open(info_file, "r") as f:
            lines = f.readlines()

        # Count number of pockets
        pocket_headers = [line for line in lines if re.match(r"Pocket \d+:", line)]
        data["num_pockets"] = len(pocket_headers)

        # Extract top pocket info only if pockets exist
        if data["num_pockets"] > 0:
            in_top_pocket = False
            for line in lines:
                if re.match(r"Pocket 1:", line):
                    in_top_pocket = True
                    continue
                elif re.match(r"Pocket \d+:", line) and not line.startswith("Pocket 1:"):
                    break  # stop after Pocket 1

                if in_top_pocket:
                    if "Druggability Score" in line:
                        data["top_pocket"]["druggability_score"] = float(line.strip().split(":")[1])
                    elif "Volume :" in line:
                        data["top_pocket"]["volume"] = float(line.strip().split(":")[1])
                    elif "Hydrophobicity score" in line:
                        data["top_pocket"]["hydrophobicity_score"] = float(line.strip().split(":")[1])
                    elif "Polarity score" in line:
                        data["top_pocket"]["polarity_score"] = float(line.strip().split(":")[1])
                    elif "Total SASA" in line:
                        data["top_pocket"]["total_sasa"] = float(line.strip().split(":")[1])
    except Exception as e:
        print(f"⚠️ Error reading {info_file}: {e}")
        return None

    return data

# Step 1: Collect data from all *_out folders
results = []
for folder in os.listdir(fpocket_output_base):
    folder_path = os.path.join(fpocket_output_base, folder)
    if os.path.isdir(folder_path) and folder.endswith("_out"):
        protein_name = folder.replace("_out", ".pdb")
        info_file = os.path.join(folder_path, folder.replace("_out", "_info.txt"))

        if os.path.exists(info_file):
            pocket_data = extract_fpocket_info(info_file)
            if pocket_data:
                top = pocket_data["top_pocket"]
                results.append({
                    "Protein": protein_name,
                    "Number of Pockets": pocket_data.get("num_pockets", 0),
                    "Druggability Score": top.get("druggability_score", 0),
                    "Volume": top.get("volume", 0),
                    "Hydrophobicity Score": top.get("hydrophobicity_score", 0),
                    "Polarity Score": top.get("polarity_score", 0),
                    "Total SASA": top.get("total_sasa", 0)
                })

# Step 2: Normalize and score
def normalize(values):
    min_val = min(values)
    max_val = max(values)
    return [(v - min_val) / (max_val - min_val) if max_val > min_val else 0 for v in values]

# Extract values to normalize
drug_scores = [r["Druggability Score"] for r in results]
volumes = [r["Volume"] for r in results]
num_pockets = [r["Number of Pockets"] for r in results]

norm_drug = normalize(drug_scores)
norm_vol = normalize(volumes)
norm_poc = normalize(num_pockets)

# Compute composite score
for i, r in enumerate(results):
    score = 0.6 * norm_drug[i] + 0.3 * norm_vol[i] + 0.1 * norm_poc[i]
    r["Composite Score"] = round(score, 4)

# Sort and save
results = sorted(results, key=lambda x: x["Composite Score"], reverse=True)

# Output file
fields = [
    "Protein", "Number of Pockets", "Druggability Score", "Volume",
    "Hydrophobicity Score", "Polarity Score", "Total SASA", "Composite Score"
]

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)

print(f"✅ Druggability summary saved to: {output_csv}")
print("🏆 Top 5 Proteins by Composite Score:")
for r in results[:5]:
    print(f"   {r['Protein']} → Score: {r['Composite Score']}")

