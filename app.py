import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol

st.set_page_config(page_title="Factor IX Mutation Explorer", layout="wide")

# --- DATA LOADING ---
@st.cache_data
def load_data():
    df = pd.read_csv("factor9_variants.csv")
    return df

df = load_data()

# --- HEADER ---
st.title("🧬 Factor IX Missense Mutation Structural Database")
st.markdown("Interactive structure-function mapping of *F9* mutations in Hemophilia B.")

# --- SIDEBAR FILTERS ---
st.sidebar.header("🔍 Search & Filter Variants")

only_missense = st.sidebar.checkbox("Show Missense Variants Only", value=True)
filtered_df = df.copy()

if only_missense:
    filtered_df = filtered_df[filtered_df['HGVS_p'].astype(str).str.startswith('p.', na=False)]

# Domain Filter
domains = ["All"] + sorted([d for d in filtered_df['Domain'].dropna().unique() if str(d).strip() != ''])
selected_domain = st.sidebar.selectbox("Filter by Domain", domains)
if selected_domain != "All":
    filtered_df = filtered_df[filtered_df['Domain'] == selected_domain]

# Severity Filter
severities = ["All"] + sorted([s for s in filtered_df['Clinical_Severity'].dropna().unique() if str(s).strip() != ''])
selected_severity = st.sidebar.selectbox("Filter by Clinical Severity", severities)
if selected_severity != "All":
    filtered_df = filtered_df[filtered_df['Clinical_Severity'] == selected_severity]

# Prepare user-friendly search labels
def make_label(row):
    hgvs = str(row.get('HGVS_p', '')).strip()
    res = str(row.get('Residue_ID', '')).strip()
    cdna = str(row.get('cDNA', '')).strip()
    if hgvs and hgvs != '-' and hgvs != 'nan':
        return f"{hgvs} (Residue {res})"
    return f"{cdna} (Non-coding)"

filtered_df['Display_Label'] = filtered_df.apply(make_label, axis=1)

variant_list = filtered_df['Display_Label'].unique().tolist()

if not variant_list:
    st.warning("No variants match the selected filters.")
    st.stop()

selected_label = st.sidebar.selectbox("Type to Search Variant (e.g., p.Val92Ala or Residue 92):", variant_list)
variant_row = filtered_df[filtered_df['Display_Label'] == selected_label].iloc[0]

# --- MAIN WORKSPACE TABS ---
tab1, tab2, tab3 = st.tabs(["🔬 3D Spatial Viewer", "📊 Structural Mechanism Cards", "📋 Full Dataset Overview"])

with tab1:
    col_left, col_right = st.columns([1, 1.2])
    
    with col_left:
        st.subheader(f"Variant: {variant_row.get('HGVS_p', 'N/A')}")
        
        m1, m2 = st.columns(2)
        res_val = pd.to_numeric(variant_row.get('Residue_ID'), errors='coerce')
        m1.metric("Residue ID", int(res_val) if pd.notnull(res_val) else "N/A")
        m2.metric("Clinical Severity", str(variant_row.get('Clinical_Severity', 'N/A')))

        m3, m4 = st.columns(2)
        dist_val = variant_row.get('Dist_To_Ser411', 'N/A')
        m3.metric("Dist. to Active Site (Ser411)", f"{dist_val} Å" if pd.notnull(dist_val) and dist_val != 'N/A' else "N/A")
        m4.metric("Disulfide Disruption", str(variant_row.get('Disulfide_Disrupt', 'N/A')))

        st.markdown(f"**Domain:** {variant_row.get('Domain', 'N/A')}")
        
    with col_right:
        st.subheader("Interactive 3D Protein Structure")
        
        xyzview = py3Dmol.view(query='pdb:1CFH')
        xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
        xyzview.addStyle({'resi': '411'}, {'stick': {'color': 'yellow', 'radius': 0.3}})
        
        if pd.notnull(res_val) and res_val > 0:
            resi_str = str(int(res_val))
            xyzview.addStyle({'resi': resi_str}, {'stick': {'color': 'red', 'radius': 0.4}})
            xyzview.zoomTo({'resi': resi_str})
        else:
            xyzview.zoomTo()
            
        showmol(xyzview, height=450, width=550)
        st.caption("🔴 Red sticks = Selected Mutation Site | 🟡 Yellow sticks = Catalytic Active Site (Ser411)")

with tab2:
    st.subheader(f"Mechanism Breakdown for {variant_row.get('HGVS_p', 'N/A')}")
    
    col_a, col_b = st.columns(2)
    col_a.markdown(f"**Domain:** `{variant_row.get('Domain', 'N/A')}`")
    col_b.markdown(f"**Severity:** `{variant_row.get('Clinical_Severity', 'N/A')}`")
    
    st.divider()
    
    st.markdown("### 🔍 Observation")
    st.info(variant_row.get('Observation', 'No observation generated.'))
    
    st.markdown("### ⚡ Structural Impact")
    st.warning(variant_row.get('Structural_Impact', 'No impact generated.'))

with tab3:
    st.subheader("Filtered Dataset Overview")
    st.dataframe(filtered_df, use_container_width=True)
