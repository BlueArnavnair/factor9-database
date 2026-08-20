import streamlit as st
import pandas as pd
import py3Dmol
from stmol import showmol
import plotly.express as px

st.set_page_config(page_title="Factor IX Mutation Structural Database", layout="wide")

# --- DATA LOADING & CLEANING ---
@st.cache_data
def load_data():
    df = pd.read_csv("factor9_variants.csv")
    
    # Standardize column searching
    cols = df.columns.tolist()
    
    def clean_severity(val):
        s = str(val).strip().lower()
        if pd.isna(val) or s in ['nan', '', 'unclassified', 'n/a']:
            return 'Unclassified'
        if 'severe' in s:
            return 'Severe'
        if 'moderate' in s:
            return 'Moderate'
        if 'mild' in s:
            return 'Mild'
        try:
            num = float(val)
            if num < 1.0: return 'Severe'
            elif 1.0 <= num <= 5.0: return 'Moderate'
            elif num > 5.0: return 'Mild'
        except ValueError:
            pass
        return str(val).strip().capitalize()

    # Domain cleaning
    domain_col = next((c for c in cols if 'domain' in c.lower()), None)
    if domain_col:
        df['Domain'] = df[domain_col].fillna('Unassigned Region')
    else:
        df['Domain'] = 'Unassigned Region'

    # Severity cleaning
    sev_col = next((c for c in cols if 'sever' in c.lower() or 'clinical' in c.lower()), None)
    if sev_col:
        df['Severity_Clean'] = df[sev_col].apply(clean_severity)
    else:
        df['Severity_Clean'] = 'Unclassified'

    return df

df = load_data()

# --- TOP NAVIGATION TABS ---
nav_landing, nav_explorer = st.tabs(["🏠 Overview & Domain Insights", "🔬 Variant Explorer"])

# ==========================================
# PAGE 1: LANDING PAGE & DOMAIN INSIGHTS
# ==========================================
with nav_landing:
    st.title("🧬 Factor IX Structural Mutation Portal")
    st.markdown("### Structure-Function Relationship in Hemophilia B")
    
    st.info(
        "**Core Finding:** The clinical severity of Hemophilia B correlates not with the physical size of "
        "the amino acid substitution, but with the **structural importance of the mutation's location** "
        "(e.g., critical domain interfaces, calcium-binding coordinates, and catalytic triads)."
    )
    
    st.divider()
    
    col_chart1, col_chart2 = st.columns(2)
    
    with col_chart1:
        st.subheader("📊 Mutation Distribution across Factor IX Domains")
        domain_counts = df['Domain'].value_counts().reset_index()
        domain_counts.columns = ['Domain', 'Count']
        
        fig_pie = px.pie(
            domain_counts, 
            names='Domain', 
            values='Count',
            hole=0.4,
            color_discrete_sequence=px.colors.qualitative.Pastel
        )
        fig_pie.update_layout(margin=dict(t=20, b=20, l=10, r=10))
        st.plotly_chart(fig_pie, use_container_width=True)
        
    with col_chart2:
        st.subheader("⚡ Severity Profiles by Protein Region")
        chart_df = df[df['Severity_Clean'].isin(['Severe', 'Moderate', 'Mild'])]
        
        if not chart_df.empty:
            fig_bar = px.histogram(
                chart_df,
                x='Domain',
                color='Severity_Clean',
                barmode='stack',
                color_discrete_map={'Severe': '#EF553B', 'Moderate': '#FECB52', 'Mild': '#00CC96'},
                category_orders={'Severity_Clean': ['Severe', 'Moderate', 'Mild']}
            )
            fig_bar.update_layout(
                xaxis_title="Protein Domain", 
                yaxis_title="Number of Variants",
                legend_title="Severity",
                margin=dict(t=20, b=20, l=10, r=10)
            )
            st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.write("No severity records available for bar chart display.")

    st.divider()
    
    st.subheader("🔮 Full Structural Architecture of Human Factor IX")
    st.caption("Complete 3D tertiary fold showing domain organization and catalytic machinery.")
    
    view_main = py3Dmol.view(query='pdb:1CFH')
    view_main.setStyle({'cartoon': {'color': 'spectrum'}})
    view_main.addStyle({'resi': '411'}, {'stick': {'color': 'yellow', 'radius': 0.3}})
    view_main.zoomTo()
    showmol(view_main, height=450, width=1000)

# ==========================================
# PAGE 2: VARIANT EXPLORER
# ==========================================
with nav_explorer:
    st.sidebar.header("🔍 Filter Variants")

    only_missense = st.sidebar.checkbox("Show Missense Variants Only", value=True)
    filtered_df = df.copy()

    # Safely create display labels across variable CSV columns
    def make_label(row):
        for col_name in ['HGVS_p', 'Mutation', 'Variant', 'Change', 'p_change']:
            if col_name in row and pd.notnull(row[col_name]) and str(row[col_name]).strip() not in ['', 'nan', 'N/A']:
                res = str(row.get('Residue_ID', '')).strip()
                return f"{row[col_name]} (Residue {res})" if res else str(row[col_name])
        
        res = str(row.get('Residue_ID', '')).strip()
        return f"Residue {res}" if res else "Unknown Variant"

    filtered_df['Display_Label'] = filtered_df.apply(make_label, axis=1)

    if only_missense:
        filtered_df = filtered_df[
            filtered_df['Display_Label'].str.contains('p\.|>|Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val', regex=True, na=False)
        ]

    domains = ["All"] + sorted([d for d in filtered_df['Domain'].dropna().unique() if str(d).strip() != ''])
    selected_domain = st.sidebar.selectbox("Filter by Domain", domains)
    if selected_domain != "All":
        filtered_df = filtered_df[filtered_df['Domain'] == selected_domain]

    severities = ["All"] + sorted([s for s in filtered_df['Severity_Clean'].dropna().unique() if str(s).strip() != ''])
    selected_severity = st.sidebar.selectbox("Filter by Severity", severities)
    if selected_severity != "All":
        filtered_df = filtered_df[filtered_df['Severity_Clean'] == selected_severity]

    variant_list = filtered_df['Display_Label'].unique().tolist()

    if not variant_list:
        st.warning("No variants match the selected filters.")
        st.stop()

    selected_label = st.sidebar.selectbox("Type to Search Variant:", variant_list)
    variant_row = filtered_df[filtered_df['Display_Label'] == selected_label].iloc[0]

    st.subheader(f"Variant Profile: {selected_label}")
    
    tab_3d, tab_cards, tab_table = st.tabs(["🔬 3D Spatial Viewer", "📊 Structural Mechanism Cards", "📋 Full Dataset Overview"])

    with tab_3d:
        col_left, col_right = st.columns([1, 1.2])
        
        with col_left:
            m1, m2 = st.columns(2)
            res_val = pd.to_numeric(variant_row.get('Residue_ID'), errors='coerce')
            m1.metric("Residue ID", int(res_val) if pd.notnull(res_val) else "N/A")
            m2.metric("Clinical Severity", str(variant_row.get('Severity_Clean', 'N/A')))

            m3, m4 = st.columns(2)
            dist_val = variant_row.get('Dist_To_Ser411', 'N/A')
            m3.metric("Dist. to Active Site (Ser411)", f"{dist_val} Å" if pd.notnull(dist_val) and str(dist_val) != 'N/A' else "N/A")
            m4.metric("Disulfide Disruption", str(variant_row.get('Disulfide_Disrupt', 'N/A')))

            st.markdown(f"**Domain:** {variant_row.get('Domain', 'N/A')}")
            
        with col_right:
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
            st.caption("🔴 Red sticks = Mutation Site | 🟡 Yellow sticks = Catalytic Active Site (Ser411)")

    with tab_cards:
        st.markdown(f"### Mechanism Breakdown for `{selected_label}`")
        
        col_a, col_b = st.columns(2)
        col_a.markdown(f"**Domain:** `{variant_row.get('Domain', 'N/A')}`")
        col_b.markdown(f"**Severity:** `{variant_row.get('Severity_Clean', 'N/A')}`")
        
        st.divider()
        
        st.markdown("### 🔍 Observation")
        st.info(variant_row.get('Observation', 'No observation generated.'))
        
        st.markdown("### ⚡ Structural Impact")
        st.warning(variant_row.get('Structural_Impact', 'No impact generated.'))

    with tab_table:
        st.subheader("Filtered Dataset Overview")
        st.dataframe(filtered_df, use_container_width=True)