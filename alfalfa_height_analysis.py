import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sns.set(style="whitegrid")
plt.rcParams['font.sans-serif'] = ['SimHei'] # For Chinese characters
plt.rcParams['axes.unicode_minus'] = False

# File paths
files = {
    "2001": "2005.12.7中澳苜蓿项目数据汇总.xlsx - 2001data.csv",
    "2002": "2005.12.7中澳苜蓿项目数据汇总.xlsx - 2002data .csv",
    "2003": "2005.12.7中澳苜蓿项目数据汇总.xlsx - 2003data.csv",
    "2004": "2005.12.7中澳苜蓿项目数据汇总.xlsx - 2004data.csv",
    "2005": "2005.12.7中澳苜蓿项目数据汇总.xlsx - 2005data.csv"
}

data_frames = {}

# Helper function to load and clean
def load_and_clean(year, filepath):
    # Load raw to inspect structure
    df_raw = pd.read_csv(filepath, header=None)
    
    # Identify Line/Variety column
    # Usually column 1 or 2. Let's look for "Line" or "Variety" in first few rows
    line_col_idx = -1
    for c in range(5):
        if df_raw[c].astype(str).str.contains("Line|Variety", case=False).any():
            line_col_idx = c
            break
    if line_col_idx == -1: line_col_idx = 1 # Fallback

    # Extract Data
    # 2001 is simple
    if year == "2001":
        # Look for "FH1"
        fh1_col = -1
        for r in range(5):
            for c in range(len(df_raw.columns)):
                val = str(df_raw.iloc[r, c])
                if "FH1" in val:
                    fh1_col = c
                    break
            if fh1_col != -1: break
        
        if fh1_col != -1:
            data = df_raw.iloc[3:, [line_col_idx, fh1_col]].copy()
            data.columns = ["Variety", "Fall_Height"]
            data["Cutting_Height"] = np.nan # No cutting height in 2001 usually
        else:
            return None

    # 2002-2005
    else:
        # Find "Mean" or "平均" columns
        avg_cols = []
        # We need to distinguish Cutting Avg vs Fall Avg
        # Strategy: Iterate columns, find those with "平均" or "Mean" in row 2 or 3
        # Then check the header above it (row 0 or 1) for "刈割" (Cutting) or "秋眠" (Fall)
        
        cutting_col = -1
        fall_col = -1
        
        # Search rows 0-3 for keywords
        for c in range(len(df_raw.columns)):
            col_content = df_raw.iloc[0:4, c].astype(str).str.cat(sep=" ")
            if "平均" in col_content or "Mean" in col_content:
                # Check context
                # Look at previous columns in the same block or header above
                # This is tricky with CSV. Let's look at the specific file structure from snippets.
                # In snippets: 2002 has "刈割高度" then "鲜重" then "秋眠"
                # Actually, simpler approach: The first "Average" is usually Cutting (Spring/Summer), second is Fall (if exists)
                
                # Check if "秋眠" is in the column header stack
                if "秋眠" in col_content:
                     if fall_col == -1: fall_col = c
                # Check if "刈割" is in the column header stack
                elif "刈割" in col_content:
                    if cutting_col == -1: cutting_col = c
                # Fallback based on order if keywords specific to column are missing but header spans
                # We'll refine this if needed.
        
        # Refined search for 2002-2005 based on snippets
        # 2002: Cutting Avg is ~col 19? Fall Avg ~ col ?
        # Let's try to map dynamically
        
        indices = [line_col_idx]
        names = ["Variety"]
        
        if cutting_col != -1:
            indices.append(cutting_col)
            names.append("Cutting_Height")
        if fall_col != -1:
            indices.append(fall_col)
            names.append("Fall_Height")
            
        data = df_raw.iloc[4:, indices].copy() # Assuming data starts row 4
        data.columns = names

    data["Year"] = year
    # Clean numeric
    for col in ["Cutting_Height", "Fall_Height"]:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')
    
    return data

# Process all files
all_data = []
for year, fp in files.items():
    try:
        df = load_and_clean(year, fp)
        if df is not None:
            all_data.append(df)
    except Exception as e:
        print(f"Error processing {year}: {e}")

# Combine
if all_data:
    full_df = pd.concat(all_data, ignore_index=True)
    
    # Aggregate by Variety and Year
    df_grouped = full_df.groupby(["Variety", "Year"]).mean(numeric_only=True).reset_index()
    
    # Save for user info
    print(df_grouped.head())
    print(df_grouped.info())
else:
    print("No data extracted.")

# Inspect headers more closely
for year, fp in files.items():
    print(f"--- {year} Header Inspection ---")
    df = pd.read_csv(fp, header=None, nrows=5)
    print(df.to_string())
    print("\n")

# Explicit extraction based on visual inspection
all_data = []

# 2001
df01 = pd.read_csv(files["2001"], header=None, skiprows=4)
# Column 2: Variety, Column 5: Fall Height (FH1)
# 2001 data has Fall Height but no Cutting Height usually, treating FH1 as Fall Height
d01 = df01.iloc[:, [2, 5]].copy()
d01.columns = ["Variety", "Fall_Height"]
d01["Cutting_Height"] = np.nan
d01["Year"] = "2001"
all_data.append(d01)

# 2002
df02 = pd.read_csv(files["2002"], header=None, skiprows=4)
d02 = df02.iloc[:, [2, 17, 36]].copy()
d02.columns = ["Variety", "Cutting_Height", "Fall_Height"]
d02["Year"] = "2002"
all_data.append(d02)

# 2003
df03 = pd.read_csv(files["2003"], header=None, skiprows=4)
d03 = df03.iloc[:, [2, 17, 32]].copy()
d03.columns = ["Variety", "Cutting_Height", "Fall_Height"]
d03["Year"] = "2003"
all_data.append(d03)

# 2004
df04 = pd.read_csv(files["2004"], header=None, skiprows=4)
d04 = df04.iloc[:, [2, 17, 33]].copy()
d04.columns = ["Variety", "Cutting_Height", "Fall_Height"]
d04["Year"] = "2004"
all_data.append(d04)

# 2005
df05 = pd.read_csv(files["2005"], header=None, skiprows=4)
d05 = df05.iloc[:, [2, 17, 32]].copy()
d05.columns = ["Variety", "Cutting_Height", "Fall_Height"]
d05["Year"] = "2005"
all_data.append(d05)

# Combine
full_df = pd.concat(all_data, ignore_index=True)

# Clean
full_df["Variety"] = full_df["Variety"].astype(str).str.strip()
full_df["Cutting_Height"] = pd.to_numeric(full_df["Cutting_Height"], errors='coerce')
full_df["Fall_Height"] = pd.to_numeric(full_df["Fall_Height"], errors='coerce')

# Group Analysis
# Mean height per variety across all years (ignoring NaNs)
variety_stats = full_df.groupby("Variety")[["Cutting_Height", "Fall_Height"]].mean().reset_index()

# Sort by Cutting Height (primary yield indicator)
top_cutting = variety_stats.sort_values("Cutting_Height", ascending=False).head(10)
top_fall = variety_stats.sort_values("Fall_Height", ascending=False).head(10)

print("Top 10 Cutting Height Varieties:")
print(top_cutting)
print("\nTop 10 Fall Dormancy Height Varieties:")
print(top_fall)

# Prepare for plotting
# Pivot for Heatmap or Line plot
full_pivot = full_df.groupby(["Year", "Variety"])["Cutting_Height"].mean().reset_index()

# Save processed data to csv for user download (simulated)
full_df.to_csv("processed_alfalfa_data.csv", index=False)

# Visualization

# 1. Bar Plot: Top 10 Cutting Height
plt.figure(figsize=(10, 6))
sns.barplot(x="Cutting_Height", y="Variety", data=top_cutting, palette="viridis")
plt.title("Top 10 Alfalfa Varieties by Average Cutting Height (2002-2005)")
plt.xlabel("Average Height (cm)")
plt.tight_layout()
plt.savefig("top10_cutting_height.png")
plt.close()

# 2. Line Plot: Trends for Top 5 Cutting Varieties
top5_vars = top_cutting["Variety"].head(5).tolist()
df_top5 = full_df[full_df["Variety"].isin(top5_vars) & (full_df["Year"] != "2001")] # Exclude 2001 as it has no cutting height

plt.figure(figsize=(10, 6))
sns.lineplot(x="Year", y="Cutting_Height", hue="Variety", data=df_top5, marker="o", palette="tab10")
plt.title("Yearly Cutting Height Trends for Top 5 Varieties")
plt.ylabel("Height (cm)")
plt.tight_layout()
plt.savefig("trend_top5_cutting.png")
plt.close()

# 3. Scatter: Cutting vs Fall Height
plt.figure(figsize=(8, 8))
sns.scatterplot(x="Fall_Height", y="Cutting_Height", data=variety_stats, alpha=0.6)
plt.title("Correlation: Fall Dormancy Height vs Cutting Height")
plt.xlabel("Fall Height (cm) - Proxy for Dormancy")
plt.ylabel("Cutting Height (cm)")

# Highlight top performers
for i in range(top_cutting.shape[0]):
    row = top_cutting.iloc[i]
    plt.text(row["Fall_Height"]+0.5, row["Cutting_Height"], row["Variety"], fontsize=8)

plt.tight_layout()
plt.savefig("scatter_heights.png")
plt.close()
