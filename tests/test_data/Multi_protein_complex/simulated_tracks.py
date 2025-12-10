import pandas as pd
import numpy as np
import random

# length for protein CTCF
chain_a_len = 727
chain_b_len = 727 

# --- 1. Generating Numerical Track (e.g. Disorder Score) ---
# Format: Chain, TrackName, Start, End, Value
data_numerical = []

# Chain A: Simulate segments of disordered regions with random values
for i in range(0, chain_a_len, 10):
    end = min(i + 9, chain_a_len-1)
    val = np.random.rand() # Random value between 0 and 1
    data_numerical.append(["A", "Disorder", i, end, val])

# Chain B: Simulate smoother values
for i in range(0, chain_b_len, 20):
    end = min(i + 19, chain_b_len-1)
    val = np.random.rand() * 0.5 # Lower values
    data_numerical.append(["B", "Disorder", i, end, val])

# --- 2. Generating Categorical Track (e.g. Domains) ---
data_categorical = []

# Chain A: Define several domains
domains_a = [
    (0, 100, "N-Term"),
    (200, 350, "DNA-Binding"),
    (500, 600, "Zinc-Finger"),
    (700, 726, "C-Term")
]
for start, end, name in domains_a:
    data_categorical.append(["A", "Domains", start, end, name])

# Chain B: Define a domain
data_categorical.append(["B", "Domains", 10, 90, "Linker"])

# --- Merge and Save ---
all_data = data_numerical + data_categorical
df = pd.DataFrame(all_data, columns=["chain_id", "track_name", "start", "end", "value"])

# Save as temporary BED file
bed_file_path = "/data2/AlphaFold3-SeqVisToolkit/tests/test_output/visualization/track/simulated_tracks.bed"
df.to_csv(bed_file_path, sep="\t", index=False)

print(f"Simulated data saved to: {bed_file_path}")
print(df.head())