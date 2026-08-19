import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

# Page Configuration
st.set_page_config(
    page_title="Factor IX Structural Mutation Database",
    page_icon="🧬",
    layout="wide"
)

# Title & Description
st.title("🧬 Factor IX Missense Mutation Structural Database")
st.markdown("""
An interactive structure-function database mapping missense variants in human **Coagulation Factor IX (F9)**[cite: 1].
Explore 3D spatial metrics, active site proximity, and structural disruption predictions linked to clinical Hemophilia B phenotypes[cite: 1].
""")

# Load Dataset
@st.cache_data
def load_data():
    return pd.read_csv("factor9_variants.csv")

try:
    df = load_data()
except FileNotFoundError:
    st.error("Error: `factor9_variants.csv` not found in current directory. Please upload your exported CSV.")
    st.stop()

# Sidebar Filters
st.sidebar.header("🔍 Search & Filter Variants")

# Domain Filter
domains = ["All"] + list(df["Domain"].dropna().unique())
selected_domain = st.sidebar.selectbox("Filter by Domain", domains)

# Severity Filter
severities = ["All"] + list(df["Clinical_Severity"].dropna().unique())
selected_severity = st.sidebar.selectbox("Filter by Severity", severities)

# Filter Dataset
filtered_df = df.copy()
if selected_domain != "All":
    filtered_df = filtered_df[filtered_df["Domain"] == selected_domain]
if selected_severity != "All":
    filtered_df = filtered_df[filtered_df["Clinical_Severity"] == selected_severity]

# Variant Selection dropdown
variant_list = [v for v in filtered_df["HGVS_p"].dropna().unique() if str(v).strip().lower() != "nan"]

if variant_list:
    selected_variant_p = st.sidebar.selectbox("Select Variant (p.HGVS)", variant_list)
    matched_rows = filtered_df[filtered_df["HGVS_p"] == selected_variant_p]
    
    if not matched_rows.empty:
        variant_row = matched_rows.iloc[0]
    else:
        st.warning("No matching row found.")
        st.stop()
else:
    st.warning("No variants match the selected filters.")
    st.stop()

# Main Layout: 2 Columns (Left: Metrics & Details, Right: 3D Structural Viewer)
col1, col2 = st.columns([1, 1])
with col1:
    st.subheader(f"Variant Profile: {variant_row['HGVS_p']}")

    # Key Metric Cards
    m1, m2 = st.columns(2)
    
    res_val = pd.to_numeric(variant_row.get('Residue_ID'), errors='coerce')
    m1.metric("Residue ID", int(res_val) if pd.notnull(res_val) else "N/A")
    m2.metric("Clinical Severity", variant_row.get('Clinical_Severity', 'N/A'))

    m3, m4 = st.columns(2)
    dist_val = variant_row.get('Dist_To_Ser411', 'N/A')
    m3.metric("Dist. to Active Site (Ser411)", f"{dist_val} Å" if dist_val != 'N/A' else "N/A")

    disrupt_val = variant_row.get('Disulfide_Disrupt', 'N/A')
    m4.metric("Disulfide Disruption", str(disrupt_val))

    st.markdown("---")

    # Detailed Annotations
    st.markdown(f"**Domain Region:** {variant_row.get('Domain', 'N/A')}")
    st.markdown(f"**Secondary Structure:** {variant_row.get('Secondary_Structure', 'N/A')}")
    st.markdown(f"**Wild Type Residue:** `{variant_row.get('Wild_Type', 'N/A')}` → **Mutant Residue:** `{variant_row.get('Mutant', 'N/A')}`")

    st.info(f"**Structural Mechanism & Impact:**\n\n{variant_row.get('Structural_Impact_Summary', 'Location-specific structural perturbation in Factor IX structure.')}")

with col2:
    st.subheader("3D Spatial Visualization")
    
    # Render 3D Molecule using py3Dmol
    res_id = int(variant_row['Residue_ID'])
    
    try:
        with open("9cli.1.D.pdb", "r") as f:
            pdb_data = f.read()
            
        xyzview = py3Dmol.view(width=600, height=500)
        xyzview.addModel(pdb_data, 'pdb')
        
        # Style full protein backbone cartoon
        xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Highlight catalytic triad Ser411 in yellow sticks
        xyzview.addStyle({'resi': '411'}, {'stick': {'color': 'yellow', 'radius': 0.25}})
        
        # Highlight target mutated residue in red sticks
        xyzview.addStyle({'resi': str(res_id)}, {'stick': {'color': 'red', 'radius': 0.35}})
        
        # Zoom to target mutation site
        xyzview.zoomTo({'resi': str(res_id)})
        
        showmol(xyzview, height=500, width=600)
        st.caption("🔴 Red sticks: Target Mutation Site | 🟡 Yellow sticks: Catalytic Ser411")
        
    except FileNotFoundError:
        st.warning("`9cli.1.D.pdb` not found. Upload the PDB file to enable live 3D visualization.")

# Data Table Display below
st.markdown("---")
st.subheader("📋 Filtered Dataset Overview")
# Display only the columns that actually exist in your CSV
cols_to_show = [col for col in ['Residue_ID', 'HGVS_p', 'cDNA', 'Domain', 'Clinical_Severity'] if col in filtered_df.columns]
st.dataframe(filtered_df[cols_to_show] if cols_to_show else filtered_df)