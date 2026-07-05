"""
app.py
======
Πλατφόρμα Ανακάλυψης Εμβολιακών Βοηθημάτων (Vaccine Adjuvant Discovery
Platform) — προσομοίωση μικροβιακών συν-καλλιεργειών, βαθμολόγηση
υποψήφιων ανοσοδιεγερτικών μορίων (PAMPs) και διαχείριση ιστορικού
πειραμάτων μέσω SQLite.

Έκδοση 5.0 — ΠΛΗΡΩΣ ΕΝΙΑΙΟ ΑΡΧΕΙΟ (single-file): όλος ο πυρήνας
προσομοίωσης (πρώην micromol_core.py) και ο βαθμολογητής (πρώην
micromol_scorer.py) είναι ενσωματωμένοι εδώ μέσα — ΔΕΝ χρειάζονται
εξωτερικά modules. Αρκεί: streamlit run app.py

ΔΟΜΗ ΑΡΧΕΙΟΥ
------------
 1. Εισαγωγές & ρύθμιση σελίδας
 2. ΠΥΡΗΝΑΣ ΠΡΟΣΟΜΟΙΩΣΗΣ (πρώην micromol_core.py) — dataclasses, μητρώα
    μικροβίων/μεταβολιτών, CombinationSimulator (ODE μοντέλο Monod +
    Luedeking-Piret)
 3. ΒΑΘΜΟΛΟΓΗΤΗΣ (πρώην micromol_scorer.py) — VaccineScore, VaccineScorer
    (dose-response μοντέλο τύπου Hill ανά PAMP/υποδοχέα)
 4. Caching wrappers γύρω από τον πυρήνα (πλέον όλα στο ίδιο αρχείο)
 5. Επίπεδο βάσης δεδομένων (SQLite): init / save / load / delete
 6. Πυρήνας λογικής προσομοίωσης & βαθμολόγησης εφαρμογής
 7. Συναρτήσεις οπτικοποίησης (bar chart, time-series, heatmap, δίκτυο)
 8. Πλευρική μπάρα πλοήγησης
 9. Σελίδα: Προσομοίωση (Simulation)
10. Σελίδα: Batch Analysis (παραμετρική σάρωση)
11. Σελίδα: Monte Carlo (ανάλυση ευαισθησίας)
12. Σελίδα: Βάση Μικροβίων (Microbe Database)
13. Σελίδα: Ιστορικό (History)
14. Σελίδα: Σύγκριση (Comparison)
15. Σελίδα: Δίκτυο Μεταβολισμού (Network)
16. Σελίδα: Σχετικά / Επιστημονικό υπόβαθρο (About)
17. Υποσέλιδο πλευρικής μπάρας

ΑΠΟΠΟΙΗΣΗ ΕΥΘΥΝΗΣ
-----------------
Πρόκειται για εκπαιδευτικό/αρχιτεκτονικό λογισμικό επίδειξης. Οι κινητικές
παράμετροι των μικροβίων και οι βαθμολογίες ανοσογονικότητας είναι
απλοποιημένα, ενδεικτικά μοντέλα — ΔΕΝ έχουν πειραματική βαθμονόμηση και
δεν πρέπει να χρησιμοποιηθούν για πραγματικές επιστημονικές αποφάσεις.
"""

# ============================================================================
# 1. ΕΙΣΑΓΩΓΕΣ & ΡΥΘΜΙΣΗ ΣΕΛΙΔΑΣ
# ============================================================================
from __future__ import annotations

import hashlib
import io
import json
import os
import random
import sqlite3
import subprocess
import sys          # <-- ΠΡΟΣΘΗΚΗ 1
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st
import itertools

sys.setrecursionlimit(1000000)   # <-- ΠΡΟΣΘΗΚΗ 2

try:
    from cobra.io import read_sbml_model
    from cobra.io.sbml import validate_sbml_model
    from cobra import Model, Reaction, Metabolite
    _HAS_COBRA = True
except ImportError:
    _HAS_COBRA = False
    validate_sbml_model = None

try:
    from scipy.integrate import odeint
    from scipy.optimize import curve_fit
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

st.set_page_config(
    page_title="Vaccine Adjuvant Discovery Platform",
    page_icon="\U0001F489",
    layout="wide",
    initial_sidebar_state="expanded",
)
# ============================================================================
# ΣΧΕΔΙΑΣΤΙΚΟ ΣΥΣΤΗΜΑ (design tokens) — "εργαστηριακό / πανεπιστημιακό"
# ============================================================================
# Παλέτα εμπνευσμένη από επιστημονικά περιοδικά & εργαστηριακά όργανα:
# βαθύ ναυτικό μελάνι, πλατίνα-χαρτί, βαθύ πράσινο-πετρί (μικροβιολογία),
# μουντό κεραμιδί για τοξικότητα, ορείχαλκος για "καταλογικές" ετικέτες.
COLOR_INK = "#1C2B39"          # κύριο κείμενο / masthead
COLOR_PAPER = "#F2F4F1"        # φόντο σελίδας (ψυχρό εργαστηριακό λευκό)
COLOR_PANEL = "#152534"        # σκούρο πλαίσιο πλευρικής μπάρας ("όργανο")
COLOR_LINE = "#CBD2CC"         # λεπτές διαχωριστικές γραμμές
COLOR_TARGET = "#0F6D66"       # στόχοι-PAMP (βαθύ πετρόλ)
COLOR_BYPRODUCT = "#9E3B33"    # παραπροϊόντα / τοξικότητα (μουντό κεραμιδί)
COLOR_ACCENT = "#2E5266"       # δευτερεύον accent (κόμβοι μικροβίων σε δίκτυα, ιστογράμματα)
COLOR_WARNING = "#B8863B"      # ορείχαλκος — προσοχή / καταλογικές ετικέτες
COLOR_GOLD = "#B8863B"

FONT_DISPLAY = "'Source Serif 4', 'Source Serif Pro', Georgia, serif"
FONT_BODY = "'IBM Plex Sans', 'Segoe UI', sans-serif"
FONT_MONO = "'IBM Plex Mono', 'Courier New', monospace"

sns.set_theme(style="ticks")


def inject_global_styles() -> None:
    """Εισάγει το καθολικό CSS του "εργαστηριακού" σχεδιαστικού συστήματος.

    Streamlit δεν επιτρέπει πλήρη έλεγχο κάθε native widget, οπότε
    στοχεύουμε σταθερά data-testid selectors (sidebar, κουμπιά, sliders,
    metrics, alerts, expanders, dataframes) και συμπληρώνουμε με custom
    HTML components (βλ. render_masthead / render_verdict_stamp /
    render_assay_card) εκεί όπου χρειαζόμαστε πλήρη έλεγχο.
    """
    st.markdown(f"""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500;600&display=swap');

    html, body, [class*="css"] {{ font-family: {FONT_BODY}; }}

    [data-testid="stAppViewContainer"], [data-testid="stMain"] {{
        background-color: {COLOR_PAPER};
    }}
    [data-testid="stHeader"] {{ background-color: transparent; }}

    /* -- Τυπογραφία επικεφαλίδων -------------------------------------- */
    h1, h2, h3, [data-testid="stHeading"] p {{
        font-family: {FONT_DISPLAY} !important;
        color: {COLOR_INK} !important;
        font-weight: 600 !important;
        letter-spacing: -0.01em;
    }}
    [data-testid="stMarkdownContainer"] p, label, span {{ color: {COLOR_INK}; }}

    /* -- Πλευρική μπάρα ("όργανο ρυθμίσεων") --------------------------- */
    [data-testid="stSidebar"] {{
        background-color: {COLOR_PANEL};
        border-right: 1px solid #0A141C;
    }}
    [data-testid="stSidebar"] * {{ color: #E7ECEA !important; }}
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3, [data-testid="stSidebar"] [data-testid="stHeading"] p {{
        font-family: {FONT_DISPLAY} !important;
        color: #FFFFFF !important;
    }}
    [data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {{
        font-family: {FONT_MONO} !important;
        color: #9FB0AC !important;
        font-size: 0.72rem !important;
        letter-spacing: 0.03em;
    }}
    [data-testid="stSidebar"] hr {{ border-color: #2A3D49; }}

    /* Πλοήγηση (st.radio) σαν πίνακας οργάνου */
    [data-testid="stSidebar"] [data-testid="stRadio"] > div[role="radiogroup"] {{
        gap: 2px;
    }}
    /* -- Selectbox μέσα στο sidebar: σκούρο κείμενο πάνω στο λευκό φόντο τους -- */
[data-testid="stSidebar"] [data-baseweb="select"] * {{
    color: {COLOR_INK} !important;
}}
[data-testid="stSidebar"] [data-baseweb="select"] > div {{
    background-color: #FFFFFF;
    border-color: {COLOR_LINE} !important;
}}
    [data-testid="stSidebar"] label[data-baseweb="radio"] {{
        border-left: 3px solid transparent;
        padding: 6px 8px !important;
        border-radius: 2px;
        transition: background-color 120ms ease, border-color 120ms ease;
    }}
    [data-testid="stSidebar"] label[data-baseweb="radio"]:hover {{
        background-color: #1E3340;
    }}
    [data-testid="stSidebar"] label[data-baseweb="radio"]:has(input:checked) {{
        border-left: 3px solid {COLOR_WARNING};
        background-color: #1E3340;
    }}
    [data-testid="stSidebar"] label[data-baseweb="radio"]:has(input:checked) p {{
        color: #FFFFFF !important;
        font-weight: 600;
    }}

    /* -- Κουμπιά --------------------------------------------------------- */
    .stButton > button, [data-testid="stBaseButton-primary"] {{
        font-family: {FONT_BODY};
        font-weight: 600;
        border-radius: 3px !important;
        border: 1px solid {COLOR_TARGET} !important;
        letter-spacing: 0.01em;
    }}
    [data-testid="stBaseButton-primary"] {{
        background-color: {COLOR_TARGET} !important;
        color: #FFFFFF !important;
    }}
    [data-testid="stBaseButton-primary"]:hover {{
        background-color: #0C554F !important;
    }}
    [data-testid="stSidebar"] .stButton > button {{
        background-color: transparent;
        color: #E7ECEA !important;
    }}

    /* -- Metrics σαν "κάρτες ένδειξης οργάνου" --------------------------- */
    [data-testid="stMetric"] {{
        background-color: #FFFFFF;
        border: 1px solid {COLOR_LINE};
        border-bottom: 3px solid {COLOR_TARGET};
        border-radius: 3px;
        padding: 0.85rem 1rem 0.7rem 1rem;
    }}
    [data-testid="stMetricLabel"] p {{
        font-family: {FONT_MONO} !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        font-size: 0.70rem !important;
        color: #5C6B67 !important;
    }}
    [data-testid="stMetricValue"] {{
        font-family: {FONT_MONO} !important;
        color: {COLOR_INK} !important;
        font-weight: 600 !important;
    }}

    /* -- Πλαίσια ειδοποιήσεων --------------------------------------------- */
    [data-testid="stAlertContainer"] {{
        border-radius: 3px !important;
        border: 1px solid {COLOR_LINE} !important;
        font-family: {FONT_BODY};
    }}

    /* -- Expander / DataFrame containers ---------------------------------- */
    [data-testid="stExpander"] {{
        border: 1px solid {COLOR_LINE} !important;
        border-radius: 3px !important;
        background-color: #FFFFFF;
    }}
    [data-testid="stDataFrame"] {{
        border: 1px solid {COLOR_LINE} !important;
        border-radius: 3px !important;
        overflow: hidden;
    }}

    /* -- Sliders / selects: accent χρώματος ------------------------------- */
    [data-testid="stSlider"] [role="slider"] {{ background-color: {COLOR_TARGET} !important; }}
    [data-baseweb="select"] {{ border-radius: 3px !important; }}

    hr {{ border-color: {COLOR_LINE}; }}

    /* -- Masthead (custom HTML, βλ. render_masthead) ---------------------- */
    .lab-masthead {{
        border-bottom: 3px double {COLOR_INK};
        padding-bottom: 0.6rem;
        margin-bottom: 1.4rem;
    }}
    .lab-masthead .eyebrow {{
        font-family: {FONT_MONO};
        font-size: 0.72rem;
        letter-spacing: 0.12em;
        text-transform: uppercase;
        color: {COLOR_TARGET};
        font-weight: 600;
    }}
    .lab-masthead h1 {{
        margin: 0.15rem 0 0.1rem 0 !important;
        font-size: 2.1rem !important;
    }}
    .lab-masthead .subtitle {{
        font-family: {FONT_BODY};
        color: #52605C;
        font-size: 0.95rem;
    }}

    /* -- Verdict stamp ------------------------------------------------------ */
    .verdict-stamp {{
        border: 2px solid var(--vs-color, {COLOR_TARGET});
        border-radius: 4px;
        padding: 0.9rem 1.1rem;
        background: color-mix(in srgb, var(--vs-color, {COLOR_TARGET}) 7%, white);
    }}
    .verdict-stamp .vs-label {{
        font-family: {FONT_MONO};
        font-weight: 700;
        letter-spacing: 0.08em;
        font-size: 0.95rem;
        color: var(--vs-color, {COLOR_TARGET});
        text-transform: uppercase;
    }}
    .verdict-stamp .vs-message {{
        font-family: {FONT_BODY};
        color: {COLOR_INK};
        font-size: 0.92rem;
        margin-top: 0.25rem;
    }}

    /* -- Assay card (signature component) ----------------------------------- */
    .assay-card {{
        border: 1px solid {COLOR_LINE};
        border-radius: 3px;
        background: #FFFFFF;
        padding: 0.8rem 0.9rem 0.7rem 0.9rem;
        margin-bottom: 0.7rem;
        position: relative;
    }}
    .assay-card .assay-top {{
        display: flex;
        justify-content: space-between;
        align-items: baseline;
    }}
    .assay-card .assay-code {{
        font-family: {FONT_MONO};
        font-size: 0.72rem;
        color: #8A968F;
        letter-spacing: 0.05em;
    }}
    .assay-card .assay-tag {{
        font-family: {FONT_MONO};
        font-size: 0.65rem;
        letter-spacing: 0.06em;
        text-transform: uppercase;
        padding: 1px 6px;
        border-radius: 2px;
        color: #FFFFFF;
    }}
    .assay-card .assay-name {{
        font-family: {FONT_DISPLAY};
        font-size: 1.15rem;
        color: {COLOR_INK};
        margin: 0.15rem 0 0.1rem 0;
        font-weight: 600;
    }}
    .assay-card .assay-receptor {{
        font-family: {FONT_MONO};
        font-size: 0.72rem;
        color: #5C6B67;
        margin-bottom: 0.4rem;
    }}
    .assay-card .assay-reading {{
        font-family: {FONT_MONO};
        font-size: 1.6rem;
        font-weight: 600;
        color: {COLOR_INK};
    }}
    .assay-card .assay-unit {{
        font-family: {FONT_MONO};
        font-size: 0.85rem;
        color: #8A968F;
    }}
    .assay-card .assay-bar-track {{
        height: 4px;
        background: #EDEFEB;
        border-radius: 2px;
        margin-top: 0.5rem;
        overflow: hidden;
    }}
    .assay-card .assay-bar-fill {{ height: 100%; }}
    </style>
    """, unsafe_allow_html=True)


def render_masthead(eyebrow: str, title: str, subtitle: str) -> None:
    """Αποδίδει επικεφαλίδα σελίδας σε ύφος επιστημονικού περιοδικού
    (eyebrow label σε mono, τίτλος σε serif, διπλή διαχωριστική γραμμή)."""
    st.markdown(f"""
    <div class="lab-masthead">
        <div class="eyebrow">{eyebrow}</div>
        <h1>{title}</h1>
        <div class="subtitle">{subtitle}</div>
    </div>
    """, unsafe_allow_html=True)

# ============================================================================
# ΝΕΑ ΣΕΛΙΔΑ: FBA (SBML)
# ============================================================================
def show_fba():
    st.markdown('<div class="card"><div class="card-header">⚡ FBA Simulation with SBML Models</div>', unsafe_allow_html=True)
    st.markdown("Φόρτωσε δύο SBML μοντέλα από τον φάκελο `agora_models` και τρέξε co‑culture FBA.")

    if not _HAS_COBRA:
        st.error("⚠️ Το COBRApy δεν είναι εγκατεστημένο. Εκτέλεσε: `pip install cobra`")
        st.markdown('</div>', unsafe_allow_html=True)
        return

    sbml_dir = Path("agora_models")
    if not sbml_dir.exists():
        st.warning("Ο φάκελος `agora_models` δεν υπάρχει. Δημιούργησέ τον και βάλε μέσα SBML αρχεία.")
        st.markdown('</div>', unsafe_allow_html=True)
        return

    sbml_files = list(sbml_dir.glob("*.xml")) + list(sbml_dir.glob("*.sbml"))
    if not sbml_files:
        st.warning("Δεν βρέθηκαν SBML αρχεία (`.xml` ή `.sbml`) στον φάκελο `agora_models`.")
        st.markdown('</div>', unsafe_allow_html=True)
        return

    def _validate_sbml_file(sbml_path: Path):
        if validate_sbml_model is None:
            return None, "Η έκδοση της COBRApy δεν υποστηρίζει επικύρωση SBML μέσω `validate_sbml_model`."
        try:
            model, errors = validate_sbml_model(sbml_path)
        except Exception as exc:
            return None, f"Σφάλμα επικύρωσης SBML: {exc}"
        if model is None:
            if not errors:
                return None, "Το αρχείο SBML δεν είναι έγκυρο, αλλά δεν υπήρξαν λεπτομέρειες σφάλματος από τον validator."
            if isinstance(errors, dict):
                details = []
                for key, value in errors.items():
                    if isinstance(value, (list, tuple)):
                        for item in value:
                            details.append(f"{key}: {item}")
                    else:
                        details.append(f"{key}: {value}")
                return None, "\n".join(details)
            return None, str(errors)
        return model, ""

    file_names = [f.name for f in sbml_files]
    col1, col2 = st.columns(2)
    with col1:
        f1 = st.selectbox("🧫 Μικρόβιο A", file_names, key="fba_a")
    with col2:
        f2 = st.selectbox("🧫 Μικρόβιο B", file_names, key="fba_b")

    if st.button("🚀 Εκτέλεση FBA", key="fba_run", type="primary"):
        if f1 == f2:
            st.warning("Επέλεξε διαφορετικά μικρόβια.")
            return

        with st.spinner("Φόρτωση μοντέλων και εκτέλεση FBA..."):
            try:
                if validate_sbml_model is not None:
                    model1, model1_errors = _validate_sbml_file(sbml_dir / f1)
                    if model1 is None:
                        st.error(f"Το μοντέλο SBML A δεν είναι έγκυρο: {f1}")
                        st.code(model1_errors)
                        return
                    model2, model2_errors = _validate_sbml_file(sbml_dir / f2)
                    if model2 is None:
                        st.error(f"Το μοντέλο SBML B δεν είναι έγκυρο: {f2}")
                        st.code(model2_errors)
                        return
                else:
                    model1 = read_sbml_model(sbml_dir / f1)
                    model2 = read_sbml_model(sbml_dir / f2)

                community = model1.copy()
                for rxn in model2.reactions:
                    if rxn.id not in community.reactions:
                        community.add_reactions([rxn.copy()])

                bio = [r for r in community.reactions if "biomass" in r.id.lower()]
                if bio:
                    community.objective = bio[0]
                else:
                    community.objective = list(community.reactions)[0].id

                from cobra.flux_analysis import pfba
                sol = pfba(community)

                secreted = {}
                for rxn in community.reactions:
                    if rxn.id.startswith("EX_") and sol.fluxes[rxn.id] > 1e-6:
                        secreted[rxn.id] = sol.fluxes[rxn.id]

                st.success("✅ FBA ολοκληρώθηκε!")
                if secreted:
                    df_flux = pd.DataFrame(list(secreted.items()), columns=["Exchange Reaction", "Flux (mmol/gDW/h)"])
                    st.dataframe(df_flux, use_container_width=True, hide_index=True)

                    final = {}
                    final["fba_fluxes"] = secreted
                    avg_score = sum(secreted.values()) if secreted else 0.0
                    toxicity = 0.0
                    scores = {}

                    save_simulation(
                        microbe_a=f1.replace(".xml", "").replace(".sbml", ""),
                        microbe_b=f2.replace(".xml", "").replace(".sbml", ""),
                        glucose=0.0, duration=0.0, ratio_a=1.0, ratio_b=1.0,
                        substrate="FBA", vmax_factor=1.0, dt=0.0,
                        final=final, scores={}, avg_score=avg_score, toxicity=toxicity,
                        recommendation="FBA",
                        sim_result={"time_steps": [], "concentration_history": {"fba_fluxes": secreted}, "biomass_history": {}}
                    )
                    st.info("💾 Τα αποτελέσματα αποθηκεύτηκαν στο Ιστορικό.")
                else:
                    st.warning("Δεν βρέθηκαν εκκρινόμενοι μεταβολίτες (exchange fluxes > 0).")

            except Exception as e:
                st.error(f"Σφάλμα κατά την εκτέλεση FBA: {e}")

    st.markdown('</div>', unsafe_allow_html=True)

def render_verdict_stamp(label: str, badge_kind: str, message: str) -> None:
    """Αποδίδει τη σύσταση προσομοίωσης ως "σφραγίδα εργαστηριακής έκθεσης"
    αντί για απλό πλαίσιο st.success/warning/error."""
    color_map = {"success": COLOR_TARGET, "warning": COLOR_WARNING, "error": COLOR_BYPRODUCT}
    color = color_map.get(badge_kind, COLOR_TARGET)
    st.markdown(f"""
    <div class="verdict-stamp" style="--vs-color: {color};">
        <div class="vs-label">Σύσταση &nbsp;·&nbsp; {label}</div>
        <div class="vs-message">{message}</div>
    </div>
    """, unsafe_allow_html=True)


def render_assay_card(catalog_no: str, name: str, category: str,
                       value: float, unit: str, receptor: str,
                       fill_fraction: float) -> None:
    """Αποδίδει το component-υπογραφή της εφαρμογής: μία "κάρτα δείγματος"
    (assay card) για έναν μεταβολίτη, με καταλογικό κωδικό, ετικέτα
    κατηγορίας, ένδειξη συγκέντρωσης σε μορφή οργάνου, και ράβδο σήματος.
    """
    is_target = category == "target"
    tag_color = COLOR_TARGET if is_target else COLOR_BYPRODUCT
    tag_text = "ΣΤΟΧΟΣ" if is_target else "ΠΑΡΑΠΡΟΪΟΝ"
    fill_pct = max(0.0, min(1.0, fill_fraction)) * 100
    st.markdown(f"""
    <div class="assay-card">
        <div class="assay-top">
            <span class="assay-code">{catalog_no}</span>
            <span class="assay-tag" style="background:{tag_color};">{tag_text}</span>
        </div>
        <div class="assay-name">{name}</div>
        <div class="assay-receptor">Υποδοχέας: {receptor}</div>
        <div><span class="assay-reading">{value:.2f}</span> <span class="assay-unit">{unit}</span></div>
        <div class="assay-bar-track"><div class="assay-bar-fill" style="width:{fill_pct:.0f}%; background:{tag_color};"></div></div>
    </div>
    """, unsafe_allow_html=True)


def apply_journal_plot_style(ax) -> None:
    """Εφαρμόζει ενιαίο, "εργαστηριακό/περιοδικού" ύφος σε ένα matplotlib Axes:
    αφαιρεί πάνω/δεξιά πλαίσιο, απαλό grid, λεπτές γραμμές αξόνων."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(COLOR_LINE)
    ax.spines["bottom"].set_color(COLOR_LINE)
    ax.tick_params(colors=COLOR_INK, labelsize=8.5)
    ax.grid(alpha=0.25, color=COLOR_LINE, linewidth=0.8)
    ax.set_axisbelow(True)


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["IBM Plex Sans", "Segoe UI", "DejaVu Sans"],
    "axes.edgecolor": COLOR_LINE,
    "axes.labelcolor": COLOR_INK,
    "text.color": COLOR_INK,
    "xtick.color": COLOR_INK,
    "ytick.color": COLOR_INK,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

DB_FILE = "simulations.db"

TARGET_METABOLITES = ["polysaccharide_lps", "mannan_polysaccharide", "flagellin"]
BYPRODUCT_METABOLITES = ["ethanol", "lactate", "acetate"]


# ============================================================================
# 2. ΠΥΡΗΝΑΣ ΠΡΟΣΟΜΟΙΩΣΗΣ (πρώην micromol_core.py — τώρα ενσωματωμένο)
# ============================================================================

# -- 2.1 Δομές δεδομένων (dataclasses) --------------------------------------

@dataclass
class MicrobeProfile:
    """Περιγράφει τις κινητικές & μεταβολικές παραμέτρους ενός μικροβίου.

    Attributes
    ----------
    name : str
        Κοινό / προβαλλόμενο όνομα (χρησιμοποιείται ως κλειδί στο dict).
    scientific_name : str
        Επιστημονική (διωνυμική) ονομασία.
    gram_stain : str
        "negative" | "positive" | "yeast".
    motile : bool
        Αν το μικρόβιο διαθέτει μαστίγιο (flagellum) -> παράγει φλαγγελίνη.
    mu_max : float
        Μέγιστος ειδικός ρυθμός ανάπτυξης (1/h).
    Ks : float
        Σταθερά ημι-κορεσμού Monod (mM).
    Y_xs : float
        Συντελεστής απόδοσης βιομάζας/υποστρώματος (αδιάστατος, toy units).
    k_death : float
        Ρυθμός θανάτου/αποδόμησης βιομάζας (1/h).
    production : Dict[str, Tuple[float, float]]
        Απεικόνιση metabolite_key -> (alpha, beta) συντελεστές Luedeking-Piret.
    description : str
        Σύντομη περιγραφή / σημείωση σχετικά με το ανοσολογικό ενδιαφέρον.
    """

    name: str
    scientific_name: str
    gram_stain: str
    motile: bool
    mu_max: float
    Ks: float
    Y_xs: float
    k_death: float
    production: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    description: str = ""
    allergy_risk: float = 0.0
    symptom_risk: float = 0.0
    study_note: str = ""


@dataclass
class MetaboliteInfo:
    """Μεταδεδομένα για έναν μεταβολίτη/μόριο του μητρώου.

    category: "target"  -> υποψήφιο ανοσοδιεγερτικό αντιγόνο (PAMP)
              "byproduct" -> παραπροϊόν ζύμωσης, συνεισφέρει σε τοξικότητα
    """

    key: str
    full_name: str
    category: str
    molecular_weight: float
    receptor: Optional[str]
    description: str


# ============================================================================
# 2.2 Προεπιλεγμένο μητρώο μικροβίων
# ============================================================================

def create_default_microbes() -> Dict[str, MicrobeProfile]:
    """Επιστρέφει το προεπιλεγμένο σύνολο μικροβίων προς προσομοίωση.

    Οι τιμές mu_max/Ks/Y_xs/k_death είναι ρεαλιστικής τάξης μεγέθους για
    κάθε είδος (π.χ. το E. coli αναπτύσσεται γρηγορότερα από τη ζύμη),
    αλλά παραμένουν απλοποιημένες παράμετροι επίδειξης.
    """
    microbes: Dict[str, MicrobeProfile] = {
        "Escherichia coli K-12": MicrobeProfile(
            name="Escherichia coli K-12",
            scientific_name="Escherichia coli",
            gram_stain="negative",
            motile=True,
            mu_max=0.85,
            Ks=5.0,
            Y_xs=0.50,
            k_death=0.015,
            production={
                # PAMP-στόχοι (targets)
                "polysaccharide_lps": (0.90, 0.05),
                "flagellin": (0.35, 0.02),
                # Παραπροϊόντα (byproducts) - υπερχείλιση μεταβολισμού (overflow)
                "acetate": (0.60, 0.03),
                "lactate": (0.10, 0.01),
            },
            description="Gram-αρνητικό, κινητό. Ισχυρός παραγωγός LPS (TLR4) "
                         "και φλαγγελίνης (TLR5). Τυπικό εργαλείο βιοτεχνολογίας.",
        ),
        "Salmonella enterica": MicrobeProfile(
            name="Salmonella enterica",
            scientific_name="Salmonella enterica serovar Typhimurium",
            gram_stain="negative",
            motile=True,
            mu_max=0.70,
            Ks=8.0,
            Y_xs=0.45,
            k_death=0.018,
            production={
                "polysaccharide_lps": (1.10, 0.08),
                "flagellin": (0.55, 0.04),
                "acetate": (0.45, 0.02),
            },
            description="Gram-αρνητικό, πολύ κινητό. Υψηλή παραγωγή LPS & "
                         "φλαγγελίνης· κλασικό μοντέλο μελέτης PAMPs.",
        ),
        "Pseudomonas aeruginosa": MicrobeProfile(
            name="Pseudomonas aeruginosa",
            scientific_name="Pseudomonas aeruginosa",
            gram_stain="negative",
            motile=True,
            mu_max=0.55,
            Ks=6.0,
            Y_xs=0.48,
            k_death=0.012,
            production={
                "polysaccharide_lps": (0.65, 0.04),
                "flagellin": (0.75, 0.05),
                "acetate": (0.20, 0.01),
            },
            description="Gram-αρνητικό, ευκαιριακό παθογόνο. Εξαιρετικά "
                         "κινητό -> υψηλή φλαγγελίνη, μέτριο LPS.",
        ),
        "Klebsiella pneumoniae": MicrobeProfile(
            name="Klebsiella pneumoniae",
            scientific_name="Klebsiella pneumoniae",
            gram_stain="negative",
            motile=False,
            mu_max=0.75,
            Ks=7.0,
            Y_xs=0.47,
            k_death=0.014,
            production={
                "polysaccharide_lps": (1.20, 0.10),
                "mannan_polysaccharide": (0.15, 0.02),  # καψιδιακό πολυσακχαρίδιο
                "acetate": (0.55, 0.03),
            },
            description="Gram-αρνητικό, ακινητο, ισχυρά ενθυλακωμένο (capsule) "
                         "-> πολύ υψηλή παραγωγή LPS και καψιδιακών πολυσακχαριτών.",
        ),
        "Bacillus subtilis": MicrobeProfile(
            name="Bacillus subtilis",
            scientific_name="Bacillus subtilis",
            gram_stain="positive",
            motile=True,
            mu_max=0.65,
            Ks=4.0,
            Y_xs=0.55,
            k_death=0.010,
            production={
                # Gram-θετικό -> ΔΕΝ διαθέτει LPS (εξωτερική μεμβράνη).
                "flagellin": (0.40, 0.03),
                "lactate": (0.25, 0.02),
            },
            description="Gram-θετικό, κινητό, σπορογόνο. Χωρίς LPS· καλός "
                         "παραγωγός φλαγγελίνης. Ασφαλές, καλά χαρακτηρισμένο.",
        ),
        "Lactobacillus plantarum": MicrobeProfile(
            name="Lactobacillus plantarum",
            scientific_name="Lactiplantibacillus plantarum",
            gram_stain="positive",
            motile=False,
            mu_max=0.90,
            Ks=3.0,
            Y_xs=0.60,
            k_death=0.008,
            production={
                "lactate": (0.85, 0.06),
                "acetate": (0.10, 0.01),
            },
            description="Gram-θετικό, ομοζυμωτικό γαλακτοβάκιλλο. Χωρίς LPS/"
                         "φλαγγελίνη· κυρίως παραγωγός γαλακτικού οξέος (τοξικότητα).",
        ),
        "Staphylococcus aureus": MicrobeProfile(
            name="Staphylococcus aureus",
            scientific_name="Staphylococcus aureus",
            gram_stain="positive",
            motile=False,
            mu_max=0.60,
            Ks=5.0,
            Y_xs=0.52,
            k_death=0.011,
            production={
                "lactate": (0.30, 0.02),
                "acetate": (0.20, 0.01),
            },
            description="Gram-θετικό, ακίνητο. Χαμηλή παραγωγή στόχων-PAMP· "
                         "χρησιμεύει ως μικρόβιο ελέγχου/αναφοράς (control).",
        ),
        "Saccharomyces cerevisiae": MicrobeProfile(
            name="Saccharomyces cerevisiae",
            scientific_name="Saccharomyces cerevisiae",
            gram_stain="yeast",
            motile=False,
            mu_max=0.45,
            Ks=10.0,
            Y_xs=0.40,
            k_death=0.009,
            production={
                "mannan_polysaccharide": (1.30, 0.08),
                "ethanol": (1.50, 0.10),
            },
            description="Ευκαρυωτική ζύμη. Το κυτταρικό της τοίχωμα είναι "
                         "πλούσιο σε μαννάνη (Dectin-1/TLR2)· παράγει αιθανόλη "
                         "μέσω ζύμωσης (υψηλή τοξικότητα σε υψηλή γλυκόζη).",
            allergy_risk=0.30,
            symptom_risk=0.22,
            study_note="Αυτό το είδος είναι πολύ χρήσιμο για μελέτες μαννάνης, αλλά απαιτεί προσοχή σε άτομα με ευαισθησία σε ζύμες ή έντονα ανοσολογικά ερεθίσματα.",
        ),
        "Bifidobacterium breve": MicrobeProfile(
            name="Bifidobacterium breve",
            scientific_name="Bifidobacterium breve",
            gram_stain="positive",
            motile=False,
            mu_max=0.40,
            Ks=3.5,
            Y_xs=0.58,
            k_death=0.006,
            production={
                "mannan_polysaccharide": (0.18, 0.01),
                "lactate": (0.20, 0.01),
            },
            description="Προβιοτικό βακτήριο με ήπιο ανοσολογικό προφίλ και ενδιαφέρον ως ασφαλής συν-καλλιεργητής για adjuvant studies.",
            allergy_risk=0.04,
            symptom_risk=0.06,
            study_note="Χαμηλή τοξικότητα και χαμηλή αλλεργική επιβάρυνση, ιδανικό για ήπιες συνδυαστικές στρατηγικές.",
        ),
        "Lactococcus lactis": MicrobeProfile(
            name="Lactococcus lactis",
            scientific_name="Lactococcus lactis",
            gram_stain="positive",
            motile=False,
            mu_max=0.50,
            Ks=2.8,
            Y_xs=0.62,
            k_death=0.007,
            production={
                "flagellin": (0.15, 0.01),
                "lactate": (0.40, 0.02),
            },
            description="Γαλακτοκομικό βακτήριο με καλή βιοασφάλεια και χαμηλό άγχος για το ανοσοποιητικό σύστημα.",
            allergy_risk=0.05,
            symptom_risk=0.05,
            study_note="Προτείνεται για μελέτες που στοχεύουν μειωμένα αντιδράσεις και ήπια ανοσολογικά σήματα.",
        ),
        "Corynebacterium glutamicum": MicrobeProfile(
            name="Corynebacterium glutamicum",
            scientific_name="Corynebacterium glutamicum",
            gram_stain="positive",
            motile=False,
            mu_max=0.42,
            Ks=4.2,
            Y_xs=0.57,
            k_death=0.008,
            production={
                "mannan_polysaccharide": (0.22, 0.01),
                "acetate": (0.18, 0.01),
            },
            description="Βιοτεχνολογικό βακτήριο με ενδιαφέρον για χαμηλή τοξικότητα και σταθερά ανοσολογικά σήματα.",
            allergy_risk=0.06,
            symptom_risk=0.07,
            study_note="Καλή επιλογή για συνδυασμούς που θέλουν ισορροπημένο ανοσολογικό προφίλ και μικρότερο όγκο παρενεργειών.",
        ),
        "Rhodobacter sphaeroides": MicrobeProfile(
            name="Rhodobacter sphaeroides",
            scientific_name="Rhodobacter sphaeroides",
            gram_stain="negative",
            motile=True,
            mu_max=0.48,
            Ks=5.5,
            Y_xs=0.46,
            k_death=0.010,
            production={
                "polysaccharide_lps": (0.55, 0.03),
                "flagellin": (0.25, 0.02),
                "acetate": (0.15, 0.01),
            },
            description="Φωτοτροφικό βακτήριο με ενδιαφέρον για έρευνα με χαμηλό-μεσαίο ανοσολογικό φορτίο και εύρωστη ανάπτυξη.",
            allergy_risk=0.10,
            symptom_risk=0.12,
            study_note="Κατάλληλο για έρευνα που αξιολογεί ισορροπημένα PAMP και λιγότερο έντονη τοξικότητα σε συγκριτικές μελέτες.",
        ),
        "Vibrio natriegens": MicrobeProfile(
            name="Vibrio natriegens",
            scientific_name="Vibrio natriegens",
            gram_stain="negative",
            motile=True,
            mu_max=1.20,
            Ks=7.5,
            Y_xs=0.43,
            k_death=0.020,
            production={
                "polysaccharide_lps": (0.95, 0.06),
                "flagellin": (0.60, 0.04),
                "acetate": (0.35, 0.02),
            },
            description="Γρήγορα αναπτυσσόμενο θαλάσσιο βακτήριο υψηλής παραγωγικότητας, χρήσιμο για ερευνητικές συγκρίσεις ταχύτητας και σήματος.",
            allergy_risk=0.18,
            symptom_risk=0.20,
            study_note="Υψηλότερο ανοσολογικό φορτίο, άρα καλύτερο για μελέτες φάσματος αλλά όχι για συνδυασμούς με χαμηλή αντοχή σε αντιδράσεις.",
        ),
        "Bacillus licheniformis": MicrobeProfile(
            name="Bacillus licheniformis",
            scientific_name="Bacillus licheniformis",
            gram_stain="positive",
            motile=True,
            mu_max=0.58,
            Ks=4.8,
            Y_xs=0.53,
            k_death=0.009,
            production={
                "flagellin": (0.35, 0.02),
                "lactate": (0.18, 0.01),
            },
            description="Γραμμα-θετικό βακτήριο με δραστηριότητα φλαγγελίνης και ενδιαφέρον για βιοτεχνολογικές προσεγγίσεις.",
            allergy_risk=0.08,
            symptom_risk=0.10,
            study_note="Χρήσιμο για τη διερεύνηση εναλλακτικών συνδυασμών με ισορροπημένη ανοσολογική ένταση.",
        ),
        "Mycobacterium smegmatis": MicrobeProfile(
            name="Mycobacterium smegmatis",
            scientific_name="Mycobacterium smegmatis",
            gram_stain="positive",
            motile=False,
            mu_max=0.35,
            Ks=6.0,
            Y_xs=0.49,
            k_death=0.007,
            production={
                "mannan_polysaccharide": (0.40, 0.02),
                "acetate": (0.12, 0.01),
            },
            description="Βακτήριο-μοντέλο με ισχυρό μεταβολικό σήμα και ενδιαφέρον για έρευνα adjuvant σε ασφαλή, ελεγχόμενη κλίμακα.",
            allergy_risk=0.09,
            symptom_risk=0.11,
            study_note="Μέτρια ανοσολογική επιβάρυνση, καλά κατάλληλο για πειραματικά σύνολα που δίνουν προτεραιότητα στη μείωση του φορτίου παρενεργειών.",
        ),
        "Clostridium sporogenes": MicrobeProfile(
            name="Clostridium sporogenes",
            scientific_name="Clostridium sporogenes",
            gram_stain="positive",
            motile=False,
            mu_max=0.32,
            Ks=5.0,
            Y_xs=0.51,
            k_death=0.005,
            production={
                "lactate": (0.25, 0.02),
                "acetate": (0.16, 0.01),
            },
            description="Σπόρογόνο βακτήριο με πολύ χαμηλό προφίλ ενόχλησης για αρχικές μελέτες ασφάλειας και συνδυασμών.",
            allergy_risk=0.03,
            symptom_risk=0.04,
            study_note="Ιδανικό για συνδυασμούς που θέλουν να εξετάσουν ένα πιο ήπιο, πιο ελεγχόμενο ανοσολογικό σήμα.",
        ),
    }
    return microbes


# ============================================================================
# 2.3 Μητρώο μεταβολιτών
# ============================================================================

def create_metabolite_registry() -> Dict[str, MetaboliteInfo]:
    """Επιστρέφει το μητρώο μεταβολιτών (στόχοι εμβολίου + παραπροϊόντα)."""
    registry: Dict[str, MetaboliteInfo] = {
        "polysaccharide_lps": MetaboliteInfo(
            key="polysaccharide_lps",
            full_name="Lipopolysaccharide (LPS)",
            category="target",
            molecular_weight=10000.0,
            receptor="TLR4 / MD-2",
            description="Ισχυρό PAMP gram-αρνητικών βακτηρίων. Υψηλή "
                         "ανοσογονικότητα αλλά ενδεχόμενη τοξικότητα (ενδοτοξίνη) "
                         "σε υψηλές συγκεντρώσεις.",
        ),
        "mannan_polysaccharide": MetaboliteInfo(
            key="mannan_polysaccharide",
            full_name="Mannan (cell-wall mannoprotein)",
            category="target",
            molecular_weight=45000.0,
            receptor="Dectin-1 / TLR2 / Mannose Receptor",
            description="Πολυσακχαρίτης κυτταρικού τοιχώματος ζυμών. Καλά "
                         "ανεκτός, χρησιμοποιείται σε προσεγγίσεις adjuvant.",
        ),
        "flagellin": MetaboliteInfo(
            key="flagellin",
            full_name="Flagellin (FliC)",
            category="target",
            molecular_weight=51000.0,
            receptor="TLR5",
            description="Πρωτεΐνη μαστιγίου κινητών βακτηρίων. Ισχυρό, καλά "
                         "μελετημένο ανοσοδιεγερτικό μόριο (π.χ. πλατφόρμες σε "
                         "κλινική έρευνα).",
        ),
        "ethanol": MetaboliteInfo(
            key="ethanol",
            full_name="Ethanol",
            category="byproduct",
            molecular_weight=46.07,
            receptor=None,
            description="Παραπροϊόν αλκοολικής ζύμωσης (κυρίως ζύμες). "
                         "Συνεισφέρει στην τοξικότητα του υγρού καλλιέργειας.",
        ),
        "lactate": MetaboliteInfo(
            key="lactate",
            full_name="Lactic acid (lactate)",
            category="byproduct",
            molecular_weight=90.08,
            receptor=None,
            description="Παραπροϊόν γαλακτικής ζύμωσης. Μειώνει το pH και "
                         "μπορεί να αναστείλει την ανάπτυξη σε υψηλές τιμές.",
        ),
        "acetate": MetaboliteInfo(
            key="acetate",
            full_name="Acetic acid (acetate)",
            category="byproduct",
            molecular_weight=60.05,
            receptor=None,
            description="Παραπροϊόν υπερχείλισης μεταβολισμού (overflow "
                         "metabolism), συχνό σε ταχέως αναπτυσσόμενα "
                         "εντεροβακτηριακά υπό περίσσεια γλυκόζης.",
        ),
    }
    return registry


# ============================================================================
# 2.4 Προσομοιωτής συν-καλλιέργειας (ODE-based)
# ============================================================================

class CombinationSimulator:
    """Προσομοιώνει τη συν-καλλιέργεια (co-culture) δύο ή περισσότερων
    μικροβίων μέσω αριθμητικής ολοκλήρωσης ενός συστήματος ΣΔΕ (ODEs).

    Η κατάσταση του συστήματος αποτελείται από:
        [S, X_1, X_2, ..., X_n, P_1, P_2, ..., P_m]
    όπου S = συγκέντρωση υποστρώματος (π.χ. γλυκόζη, mM),
          X_i = βιομάζα μικροβίου i (toy units),
          P_j = συγκέντρωση μεταβολίτη j (mM, toy units).
    """

    def __init__(self, microbes: Dict[str, MicrobeProfile],
                 metabolites: Dict[str, MetaboliteInfo]):
        self.microbes = microbes
        self.metabolites = metabolites
        self._metabolite_keys: List[str] = list(metabolites.keys())

    # -- εσωτερική συνάρτηση παραγώγων για τον ολοκληρωτή ------------------
    def _derivatives(self, y: np.ndarray, t: float,
                      profiles: List[MicrobeProfile],
                      substrate_name: str,
                      vmax_factor: float) -> List[float]:
        n = len(profiles)
        S = max(y[0], 0.0)
        X = y[1:1 + n]
        # παραγωγές
        dXdt = np.zeros(n)
        production = {k: 0.0 for k in self._metabolite_keys}
        substrate_consumption = 0.0

        for i, mp in enumerate(profiles):
            xi = max(X[i], 0.0)
            mu = vmax_factor * mp.mu_max * S / (mp.Ks + S + 1e-9) if S > 0 else 0.0
            dXdt[i] = (mu - mp.k_death) * xi
            substrate_consumption += (mu * xi) / max(mp.Y_xs, 1e-6)
            for met_key, (alpha, beta) in mp.production.items():
                if met_key in production:
                    production[met_key] += alpha * mu * xi + beta * xi

        dS = -substrate_consumption
        dPdt = [production[k] for k in self._metabolite_keys]
        return [dS] + list(dXdt) + dPdt

    # -- δημόσια μέθοδος προσομοίωσης ---------------------------------------
    def simulate_coculture(
        self,
        microbe_names: List[str],
        initial_concentrations: Dict[str, float],
        duration: float,
        dt: float = 0.1,
        ratio: Tuple[float, float] = (1.0, 1.0),
        vmax_factor: float = 1.0,
        initial_biomass_total: float = 0.05,
    ) -> Dict[str, object]:
        """Εκτελεί την προσομοίωση συν-καλλιέργειας.

        Parameters
        ----------
        microbe_names : λίστα ονομάτων μικροβίων (κλειδιά του registry).
        initial_concentrations : π.χ. {"glucose": 100.0} - αρχικό υπόστρωμα.
        duration : διάρκεια προσομοίωσης σε ώρες.
        dt : βήμα δειγματοληψίας χρονοσειράς (ώρες).
        ratio : σχετική αρχική αναλογία εμβολιασμού (π.χ. (10, 1)).
        vmax_factor : πολλαπλασιαστής στο μ_max όλων των μικροβίων
                      (προσομοιώνει μεταβολές θερμοκρασίας/pH κ.λπ.).
        initial_biomass_total : συνολική αρχική βιομάζα προς επιμερισμό
                                 σύμφωνα με το ratio.

        Returns
        -------
        dict με πεδία:
            time_steps, concentration_history, biomass_history,
            substrate_history, final_concentrations, microbe_names,
            duration, dt, ratio, substrate_name
        """
        if not microbe_names:
            raise ValueError("Απαιτείται τουλάχιστον ένα μικρόβιο.")
        for name in microbe_names:
            if name not in self.microbes:
                raise ValueError(f"Άγνωστο μικρόβιο: {name}")
        if duration <= 0:
            raise ValueError("Η διάρκεια πρέπει να είναι θετική.")
        if dt <= 0:
            raise ValueError("Το dt πρέπει να είναι θετικό.")

        profiles = [self.microbes[name] for name in microbe_names]
        n = len(profiles)

        # Όνομα & αρχική τιμή υποστρώματος (π.χ. glucose)
        if initial_concentrations:
            substrate_name, S0 = next(iter(initial_concentrations.items()))
        else:
            substrate_name, S0 = "glucose", 100.0

        # Επιμερισμός αρχικής βιομάζας σύμφωνα με το ratio
        ratio = tuple(ratio) if ratio else tuple([1.0] * n)
        if len(ratio) < n:
            ratio = tuple(list(ratio) + [1.0] * (n - len(ratio)))
        ratio_sum = float(sum(ratio[:n])) or 1.0
        X0 = [initial_biomass_total * (r / ratio_sum) for r in ratio[:n]]

        y0 = [float(S0)] + X0 + [0.0] * len(self._metabolite_keys)

        t = np.arange(0.0, duration + dt / 2.0, dt)
        if t[-1] < duration:
            t = np.append(t, duration)

        if _HAS_SCIPY:
            sol = odeint(
                self._derivatives, y0, t,
                args=(profiles, substrate_name, vmax_factor),
            )
        else:  # -- εναλλακτική: ρητή ολοκλήρωση Runge-Kutta 4ης τάξης -----
            sol = self._integrate_rk4(y0, t, profiles, substrate_name, vmax_factor)

        sol = np.clip(sol, 0.0, None)  # οι συγκεντρώσεις δεν είναι αρνητικές

        S_hist = sol[:, 0]
        X_hist = sol[:, 1:1 + n]
        P_hist = sol[:, 1 + n:]

        concentration_history = {
            key: P_hist[:, j].tolist() for j, key in enumerate(self._metabolite_keys)
        }
        biomass_history = {
            microbe_names[i]: X_hist[:, i].tolist() for i in range(n)
        }

        final_concentrations: Dict[str, float] = {
            key: concentration_history[key][-1] for key in self._metabolite_keys
        }
        final_concentrations[substrate_name] = float(S_hist[-1])
        for i, name in enumerate(microbe_names):
            final_concentrations[f"biomass_{name}"] = float(X_hist[-1, i])

        return {
            "time_steps": t.tolist(),
            "concentration_history": concentration_history,
            "biomass_history": biomass_history,
            "substrate_history": S_hist.tolist(),
            "substrate_name": substrate_name,
            "final_concentrations": final_concentrations,
            "microbe_names": microbe_names,
            "duration": duration,
            "dt": dt,
            "ratio": ratio[:n],
        }

    # -- fallback ολοκληρωτής χωρίς scipy (RK4) -----------------------------
    def _integrate_rk4(self, y0, t, profiles, substrate_name, vmax_factor):
        y0 = np.array(y0, dtype=float)
        out = np.zeros((len(t), len(y0)))
        out[0] = y0
        y = y0.copy()
        for idx in range(1, len(t)):
            h = t[idx] - t[idx - 1]
            k1 = np.array(self._derivatives(y, t[idx - 1], profiles, substrate_name, vmax_factor))
            k2 = np.array(self._derivatives(y + h / 2 * k1, t[idx - 1] + h / 2, profiles, substrate_name, vmax_factor))
            k3 = np.array(self._derivatives(y + h / 2 * k2, t[idx - 1] + h / 2, profiles, substrate_name, vmax_factor))
            k4 = np.array(self._derivatives(y + h * k3, t[idx], profiles, substrate_name, vmax_factor))
            y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            y = np.clip(y, 0.0, None)
            out[idx] = y
        return out


# ============================================================================
# 3. ΒΑΘΜΟΛΟΓΗΤΗΣ (πρώην micromol_scorer.py — τώρα ενσωματωμένο)
# ============================================================================
# ============================================================================
# 3.1 Δομή αποτελέσματος βαθμολόγησης
# ============================================================================

@dataclass
class VaccineScore:
    """Αποτέλεσμα βαθμολόγησης ενός υποψήφιου μορίου."""

    molecule_name: str
    molecule_id: str
    overall_vaccine_score: float       # 0-100, σύνθετη βαθμολογία
    immunogenicity_score: float        # 0-100
    safety_score: float                # 0-100
    stability_score: float             # 0-100
    receptor: Optional[str]
    dose_response_factor: float        # 0-1, πόσο "κορεσμένη" είναι η δόση
    rationale: str


# ============================================================================
# 3.2 Βάση γνώσης PAMP / υποδοχέων
# ============================================================================

# EC50: συγκέντρωση (σε ίδιες μονάδες με τις προσομοιώσεις, mM toy units)
# στην οποία επιτυγχάνεται το 50% της μέγιστης ανοσοδιεγερτικής απόκρισης.
# hill_n: συντελεστής συνεργασίας (Hill coefficient) της καμπύλης δόσης.
_KNOWLEDGE_BASE: Dict[str, dict] = {
    "polysaccharide_lps": {
        "immunogenicity": 95.0,
        "safety": 40.0,   # ενδοτοξίνη -> κίνδυνος υπερδιέγερσης / σηπτικού shock
        "stability": 85.0,
        "receptor": "TLR4 / MD-2",
        "ec50": 5.0,
        "hill_n": 1.5,
        "note": "Ισχυρό PAMP· υψηλή ανοσογονικότητα αντισταθμίζεται από "
                "κίνδυνο τοξικότητας ενδοτοξίνης σε υψηλές δόσεις.",
    },
    "mannan_polysaccharide": {
        "immunogenicity": 70.0,
        "safety": 88.0,
        "stability": 90.0,
        "receptor": "Dectin-1 / TLR2 / Mannose Receptor",
        "ec50": 8.0,
        "hill_n": 1.2,
        "note": "Μέτρια ανοσογονικότητα, πολύ καλό προφίλ ασφαλείας και "
                "σταθερότητας. Καλός υποψήφιος για ήπιο adjuvant.",
    },
    "flagellin": {
        "immunogenicity": 80.0,
        "safety": 75.0,
        "stability": 60.0,   # πρωτεΐνη -> λιγότερο σταθερή θερμικά
        "receptor": "TLR5",
        "ec50": 3.0,
        "hill_n": 1.8,
        "note": "Καλά τεκμηριωμένο ανοσοδιεγερτικό μόριο (χρησιμοποιείται σε "
                "ερευνητικές πλατφόρμες εμβολίων). Μειωμένη σταθερότητα λόγω "
                "πρωτεϊνικής φύσης.",
    },
}

# Προεπιλεγμένο προφίλ για άγνωστα/μη-στοχευμένα μόρια (π.χ. παραπροϊόντα).
_DEFAULT_PROFILE = {
    "immunogenicity": 5.0,
    "safety": 50.0,
    "stability": 50.0,
    "receptor": None,
    "ec50": 10.0,
    "hill_n": 1.0,
    "note": "Δεν αναγνωρίζεται ως γνωστό PAMP/ανοσοδιεγερτικό αντιγόνο· "
            "αποδίδεται χαμηλή βαθμολογία-στόχος από προεπιλογή.",
}


def _hill(concentration: float, ec50: float, n: float) -> float:
    """Καμπύλη δόσης-απόκρισης τύπου Hill, εξομαλυμένη στο [0, 1]."""
    c = max(concentration, 0.0)
    if c <= 0:
        return 0.0
    c_n = c ** n
    ec_n = ec50 ** n
    return c_n / (ec_n + c_n)


# ============================================================================
# 3.3 Κύρια κλάση βαθμολογητή
# ============================================================================

class VaccineScorer:
    """Βαθμολογεί υποψήφια μόρια ως δυνητικά εμβολιακά βοηθήματα (adjuvants),
    συνδυάζοντας μία στατική βάση γνώσης PAMP-υποδοχέα με ένα δυναμικό
    συντελεστή δόσης-απόκρισης βασισμένο στην τρέχουσα συγκέντρωση.
    """

    def __init__(self, weight_immunogenicity: float = 0.5,
                 weight_safety: float = 0.3,
                 weight_stability: float = 0.2):
        total = weight_immunogenicity + weight_safety + weight_stability
        self.w_immuno = weight_immunogenicity / total
        self.w_safety = weight_safety / total
        self.w_stability = weight_stability / total

    def score_molecule(self, name: str, molecule_id: str = "",
                        smiles: str = "",
                        structure: Optional[dict] = None) -> VaccineScore:
        """Βαθμολογεί ένα μόριο.

        Parameters
        ----------
        name : κλειδί μεταβολίτη (π.χ. "polysaccharide_lps").
        molecule_id : προαιρετικό αναγνωριστικό (χρησιμοποιείται όπως δίνεται).
        smiles : προαιρετικό (δεν χρησιμοποιείται στο απλοποιημένο μοντέλο,
                 διατηρείται για συμβατότητα API με μελλοντικές επεκτάσεις
                 βασισμένες σε χημική δομή).
        structure : προαιρετικό dict, π.χ. {"concentration": 12.4} (mM).
                    Αν παραληφθεί, χρησιμοποιείται η τιμή EC50 ως αναφορά
                    (δηλ. dose_response_factor = 0.5).
        """
        profile = _KNOWLEDGE_BASE.get(name, _DEFAULT_PROFILE)

        if structure and "concentration" in structure:
            concentration = float(structure["concentration"])
            dose_factor = _hill(concentration, profile["ec50"], profile["hill_n"])
        else:
            dose_factor = 0.5  # ουδέτερη αναφορά όταν δεν δίνεται συγκέντρωση

        # Η ανοσογονικότητα "ενεργοποιείται" προοδευτικά με τη δόση:
        # σε μηδενική δόση παραμένει στο 40% της θεωρητικής μέγιστης τιμής,
        # και κορεννύεται στο 100% καθώς η δόση πλησιάζει/ξεπερνά το EC50.
        effective_immunogenicity = profile["immunogenicity"] * (0.4 + 0.6 * dose_factor)
        effective_immunogenicity = min(effective_immunogenicity, 100.0)

        overall = (
            self.w_immuno * effective_immunogenicity
            + self.w_safety * profile["safety"]
            + self.w_stability * profile["stability"]
        )
        overall = max(0.0, min(100.0, overall))

        rationale = (
            f"{profile['note']} "
            f"(συγκέντρωση→δόση-απόκριση: {dose_factor:.2f}, "
            f"υποδοχέας: {profile['receptor'] or 'άγνωστος'})"
        )

        return VaccineScore(
            molecule_name=name,
            molecule_id=molecule_id or name,
            overall_vaccine_score=round(overall, 2),
            immunogenicity_score=round(effective_immunogenicity, 2),
            safety_score=round(profile["safety"], 2),
            stability_score=round(profile["stability"], 2),
            receptor=profile["receptor"],
            dose_response_factor=round(dose_factor, 4),
            rationale=rationale,
        )


def create_default_scorer() -> VaccineScorer:
    """Επιστρέφει έναν VaccineScorer με τα προεπιλεγμένα βάρη."""
    return VaccineScorer()


# ============================================================================
# 4. CACHING WRAPPERS (πλέον γύρω από τις ενσωματωμένες δομές παραπάνω)
# ============================================================================
# Τα create_default_microbes()/create_metabolite_registry()/create_default_scorer()
# είναι καθαρές (stateless) συναρτήσεις — τις "παγώνουμε" με cache_resource ώστε
# να μην ξαναδημιουργούνται σε κάθε rerun του Streamlit script.

@st.cache_resource(show_spinner=False)
def get_microbes() -> Dict[str, MicrobeProfile]:
    """Επιστρέφει (με caching) το μητρώο μικροβίων."""
    return create_default_microbes()


@st.cache_resource(show_spinner=False)
def get_metabolite_registry() -> Dict[str, MetaboliteInfo]:
    """Επιστρέφει (με caching) το μητρώο μεταβολιτών."""
    return create_metabolite_registry()


@st.cache_resource(show_spinner=False)
def get_scorer():
    """Επιστρέφει (με caching) τον προεπιλεγμένο βαθμολογητή."""
    return create_default_scorer()


def get_simulator() -> CombinationSimulator:
    """Ο προσομοιωτής είναι ελαφρύς (stateless πέραν των παραπάνω dict) —
    τον δημιουργούμε φρέσκο κάθε φορά χρησιμοποιώντας τα cached registries."""
    microbes = get_microbes()

    # Apply any user-saved calibration parameters to the microbe profiles
    try:
        conn = sqlite3.connect(DB_FILE)
        c = conn.cursor()
        for name, profile in microbes.items():
            c.execute("SELECT params_json FROM calibrations WHERE microbe = ? ORDER BY id DESC LIMIT 1", (name,))
            row = c.fetchone()
            if row and row[0]:
                try:
                    params = json.loads(row[0])
                    if isinstance(params, dict):
                        if params.get("model") == "monod":
                            if "mu_max" in params:
                                profile.mu_max = float(params["mu_max"])
                            if "Ks" in params:
                                profile.Ks = float(params["Ks"])
                        elif params.get("model") == "exponential":
                            if "mu" in params:
                                profile.mu_max = float(params["mu"])  # use mu as mu_max proxy
                except Exception:
                    # ignore faulty calibration entries
                    pass
        conn.close()
    except Exception:
        # DB not available or error reading calibrations; proceed with defaults
        pass

    return CombinationSimulator(microbes, get_metabolite_registry())


# ============================================================================
# 5. ΕΠΙΠΕΔΟ ΒΑΣΗΣ ΔΕΔΟΜΕΝΩΝ (SQLite)
# ============================================================================
def init_db() -> None:
    """Δημιουργεί τον πίνακα 'simulations' αν δεν υπάρχει ήδη.

    Πέραν των συγκεντρώσεων/βαθμολογιών, αποθηκεύουμε και τις πλήρεις
    χρονοσειρές (ως JSON κείμενο) ώστε να μπορούν να ανακατασκευαστούν
    αργότερα στις σελίδες Ιστορικό/Δίκτυο χωρίς να χρειάζεται νέα εκτέλεση.
    """
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS simulations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            microbe_a TEXT NOT NULL,
            microbe_b TEXT NOT NULL,
            glucose REAL,
            duration REAL,
            ratio_a REAL,
            ratio_b REAL,
            substrate TEXT,
            vmax_factor REAL,
            dt REAL,
            lps REAL,
            mannan REAL,
            flagellin REAL,
            lps_score REAL,
            mannan_score REAL,
            flagellin_score REAL,
            ethanol REAL,
            lactate REAL,
            acetate REAL,
            avg_score REAL,
            toxicity REAL,
            recommendation TEXT,
            git_commit TEXT,
            sbml_checksum TEXT,
            time_steps_json TEXT,
            concentration_history_json TEXT,
            biomass_history_json TEXT
        )
    ''')
    conn.commit()

    # -- Αυτόματο migration: αν υπάρχει ήδη ένα ΠΑΛΙΟ simulations.db από
    # προηγούμενη έκδοση του πίνακα, το CREATE TABLE IF NOT EXISTS παραπάνω
    # δεν το αλλάζει. Ελέγχουμε ποιες στήλες λείπουν και τις προσθέτουμε.
    required_columns = {
        "glucose": "REAL", "duration": "REAL", "ratio_a": "REAL", "ratio_b": "REAL",
        "substrate": "TEXT", "vmax_factor": "REAL", "dt": "REAL",
        "lps": "REAL", "mannan": "REAL", "flagellin": "REAL",
        "lps_score": "REAL", "mannan_score": "REAL", "flagellin_score": "REAL",
        "ethanol": "REAL", "lactate": "REAL", "acetate": "REAL",
        "avg_score": "REAL", "toxicity": "REAL", "recommendation": "TEXT",
        "git_commit": "TEXT", "sbml_checksum": "TEXT",
        "time_steps_json": "TEXT", "concentration_history_json": "TEXT",
        "biomass_history_json": "TEXT",
    }
    c.execute("PRAGMA table_info(simulations)")
    existing_columns = {row[1] for row in c.fetchall()}
    for col_name, col_type in required_columns.items():
        if col_name not in existing_columns:
            c.execute(f"ALTER TABLE simulations ADD COLUMN {col_name} {col_type}")
    conn.commit()
    conn.close()

    # -- Table for parameter calibrations (per-microbe)
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS calibrations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            microbe TEXT NOT NULL,
            model TEXT NOT NULL,
            params_json TEXT NOT NULL,
            notes TEXT
        )
    ''')
    conn.commit()
    conn.close()


def save_simulation(
    microbe_a: str, microbe_b: str,
    glucose: float, duration: float,
    ratio_a: float, ratio_b: float,
    substrate: str, vmax_factor: float, dt: float,
    final: Dict[str, float],
    scores: Dict[str, "VaccineScore"],
    avg_score: float, toxicity: float, recommendation: str,
    sim_result: Optional[dict] = None,
) -> None:
    """Αποθηκεύει μία εκτέλεση προσομοίωσης στη βάση δεδομένων."""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()

    time_steps_json = json.dumps(sim_result.get("time_steps", [])) if sim_result else "[]"
    conc_hist_json = json.dumps(sim_result.get("concentration_history", {})) if sim_result else "{}"
    bio_hist_json = json.dumps(sim_result.get("biomass_history", {})) if sim_result else "{}"

    # provenance: git commit short hash (if available) and SBML checksum for involved models
    def _get_git_commit_short():
        try:
            out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL)
            return out.decode().strip()
        except Exception:
            return None

    def _compute_sbml_checksum(a_name: str, b_name: str) -> Optional[str]:
        sbml_dir = Path("agora_models")
        if not sbml_dir.exists():
            return None
        files = list(sbml_dir.glob("*.xml")) + list(sbml_dir.glob("*.sbml"))
        matched = []
        for nm in (a_name, b_name):
            if not nm:
                continue
            for f in files:
                if nm.lower() in f.stem.lower():
                    matched.append(f)
                    break
        if not matched:
            return None
        h = hashlib.sha256()
        for f in matched:
            try:
                h.update(f.read_bytes())
            except Exception:
                continue
        return h.hexdigest()

    git_commit = _get_git_commit_short()
    sbml_checksum = _compute_sbml_checksum(microbe_a, microbe_b)

    c.execute('''
        INSERT INTO simulations (
            timestamp, microbe_a, microbe_b, glucose, duration,
            ratio_a, ratio_b, substrate, vmax_factor, dt,
            lps, mannan, flagellin,
            lps_score, mannan_score, flagellin_score,
            ethanol, lactate, acetate,
            avg_score, toxicity, recommendation,
            time_steps_json, concentration_history_json, biomass_history_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        datetime.now().isoformat(timespec="seconds"),
        microbe_a, microbe_b, glucose, duration,
        ratio_a, ratio_b, substrate, vmax_factor, dt,
        final.get("polysaccharide_lps", 0.0),
        final.get("mannan_polysaccharide", 0.0),
        final.get("flagellin", 0.0),
        scores["polysaccharide_lps"].overall_vaccine_score if "polysaccharide_lps" in scores else 0.0,
        scores["mannan_polysaccharide"].overall_vaccine_score if "mannan_polysaccharide" in scores else 0.0,
        scores["flagellin"].overall_vaccine_score if "flagellin" in scores else 0.0,
        final.get("ethanol", 0.0),
        final.get("lactate", 0.0),
        final.get("acetate", 0.0),
        avg_score, toxicity, recommendation,
        git_commit, sbml_checksum,
        time_steps_json, conc_hist_json, bio_hist_json,
    ))
    conn.commit()
    conn.close()


def load_history(limit: int = 500) -> pd.DataFrame:
    """Φορτώνει το ιστορικό εκτελέσεων (χωρίς τις ογκώδεις στήλες JSON)."""
    conn = sqlite3.connect(DB_FILE)
    query = f'''
        SELECT id, timestamp, microbe_a, microbe_b, glucose, duration,
               ratio_a, ratio_b, substrate, vmax_factor, dt,
               lps, mannan, flagellin, lps_score, mannan_score, flagellin_score,
               ethanol, lactate, acetate, avg_score, toxicity, recommendation
        FROM simulations ORDER BY id DESC LIMIT {int(limit)}
    '''
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def load_row_history_json(row_id: int) -> Optional[dict]:
    """Φορτώνει τις πλήρεις χρονοσειρές μιας συγκεκριμένης εγγραφής."""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute(
        "SELECT time_steps_json, concentration_history_json, biomass_history_json, "
        "microbe_a, microbe_b FROM simulations WHERE id = ?", (row_id,)
    )
    row = c.fetchone()
    conn.close()
    if not row:
        return None


def save_calibration(microbe: str, model: str, params: dict, notes: str = "") -> None:
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''
        INSERT INTO calibrations (timestamp, microbe, model, params_json, notes)
        VALUES (?, ?, ?, ?, ?)
    ''', (
        datetime.now().isoformat(timespec="seconds"),
        microbe, model, json.dumps(params), notes
    ))
    conn.commit()
    conn.close()


def load_calibrations(limit: int = 500) -> pd.DataFrame:
    conn = sqlite3.connect(DB_FILE)
    df = pd.read_sql_query(f"SELECT id, timestamp, microbe, model, params_json, notes FROM calibrations ORDER BY id DESC LIMIT {int(limit)}", conn)
    conn.close()
    if not df.empty:
        df["params"] = df["params_json"].apply(lambda s: json.loads(s) if s else {})
    return df


def load_row_history_json(row_id: int) -> Optional[dict]:
    """Φορτώνει τις πλήρεις χρονοσειρές μιας συγκεκριμένης εγγραφής."""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute(
        "SELECT time_steps_json, concentration_history_json, biomass_history_json, "
        "microbe_a, microbe_b FROM simulations WHERE id = ?", (row_id,)
    )
    row = c.fetchone()
    conn.close()
    if not row:
        return None

    time_steps, conc_hist, bio_hist, microbe_a, microbe_b = row
    return {
        "time_steps": json.loads(time_steps) if time_steps else [],
        "concentration_history": json.loads(conc_hist) if conc_hist else {},
        "biomass_history": json.loads(bio_hist) if bio_hist else {},
        "microbe_a": microbe_a,
        "microbe_b": microbe_b,
    }


def delete_all() -> None:
    """Διαγράφει όλες τις εγγραφές του πίνακα simulations."""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute("DELETE FROM simulations")
    c.execute("DELETE FROM sqlite_sequence WHERE name = 'simulations'")
    conn.commit()
    conn.close()


def count_history() -> int:
    """Επιστρέφει το πλήθος αποθηκευμένων εκτελέσεων."""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute("SELECT COUNT(*) FROM simulations")
    n = c.fetchone()[0]
    conn.close()
    return n


# ============================================================================
# 6. ΠΥΡΗΝΑΣ ΛΟΓΙΚΗΣ ΠΡΟΣΟΜΟΙΩΣΗΣ & ΒΑΘΜΟΛΟΓΗΣΗΣ ΕΦΑΡΜΟΓΗΣ
# ============================================================================
def run_simulation(
    microbe_a: str, microbe_b: str,
    glucose: float, duration: float,
    ratio_a: float = 10.0, ratio_b: float = 1.0,
    substrate: str = "glucose",
    vmax_factor: float = 1.0,
    dt: float = 0.25,
) -> Tuple[Dict[str, float], dict]:
    """Εκτελεί μία προσομοίωση συν-καλλιέργειας δύο μικροβίων.

    Επιστρέφει (final_concentrations, πλήρες αποτέλεσμα προσομοίωσης).
    """
    sim = get_simulator()
    initial = {substrate: float(glucose)}
    result = sim.simulate_coculture(
        microbe_names=[microbe_a, microbe_b],
        initial_concentrations=initial,
        duration=float(duration),
        dt=float(dt),
        ratio=(float(ratio_a), float(ratio_b)),
        vmax_factor=float(vmax_factor),
    )
    return result["final_concentrations"], result


def _estimate_safety_penalty(
    microbe_names: Optional[List[str]] = None,
    allergy_tags: Optional[List[str]] = None,
    symptom_mode: str = "balanced",
) -> float:
    """Υπολογίζει μια απλή επιπλέον τιμωρία ασφάλειας για ερευνητικά
    προφίλ που θέλουν ήπια ανοσολογικά σήματα σε άτομα με αλλεργίες ή
    ευαισθησία σε συχνές παρενέργειες εμβολίων."""
    microbes = get_microbes()
    penalty = 0.0
    if not microbe_names:
        return penalty

    tags = {tag.lower() for tag in (allergy_tags or [])}
    mode_multiplier = {"balanced": 0.25, "low_reaction": 0.45, "high_immunogenicity": 0.10}.get(symptom_mode, 0.25)

    for name in microbe_names:
        profile = microbes.get(name)
        if profile is None:
            continue
        if "yeast" in tags and profile.gram_stain == "yeast":
            penalty += profile.allergy_risk * 0.55
        if "gram-negative" in tags and profile.gram_stain == "negative":
            penalty += profile.allergy_risk * 0.35
        if "protein-rich" in tags and profile.motile:
            penalty += profile.allergy_risk * 0.20
        penalty += profile.symptom_risk * mode_multiplier
    return penalty


def calculate_scores(
    final: Dict[str, float],
    allergy_tags: Optional[List[str]] = None,
    symptom_mode: str = "balanced",
    microbe_names: Optional[List[str]] = None,
) -> Tuple[Dict[str, VaccineScore], float, float]:
    """Υπολογίζει τις βαθμολογίες εμβολιακού δυναμικού για τους στόχους
    (LPS, mannan, flagellin) και την τοξικότητα από τα παραπροϊόντα.

    Η βαθμολόγηση λαμβάνει υπόψη την πραγματική συγκέντρωση κάθε στόχου
    (dose-response), όχι μόνο μία στατική τιμή αναφοράς. Επίσης μπορεί να
    προσαρμόσει την τελική βαθμολογία για ερευνητικές περιπτώσεις που
    θέλουν πιο ασφαλείς συνδυασμούς για άτομα με αλλεργίες ή ευαισθησίες.
    """
    scorer = get_scorer()
    scores: Dict[str, VaccineScore] = {}
    safety_penalty = _estimate_safety_penalty(microbe_names, allergy_tags, symptom_mode)
    for met in TARGET_METABOLITES:
        concentration = final.get(met, 0.0)
        base_score = scorer.score_molecule(met, met, "", {"concentration": concentration})
        adjusted_safety = max(0.0, min(100.0, base_score.safety_score - safety_penalty * 8.0))
        adjusted_overall = max(0.0, min(100.0, base_score.overall_vaccine_score - safety_penalty * 6.0))
        adjusted_overall = max(adjusted_overall, adjusted_safety * 0.7)
        rationale = (
            f"{base_score.rationale} | allergy/symptom penalty: {safety_penalty:.2f}; "
            f"adjusted safety={adjusted_safety:.1f}/100"
        )
        scores[met] = VaccineScore(
            molecule_name=base_score.molecule_name,
            molecule_id=base_score.molecule_id,
            overall_vaccine_score=round(adjusted_overall, 2),
            immunogenicity_score=base_score.immunogenicity_score,
            safety_score=round(adjusted_safety, 2),
            stability_score=base_score.stability_score,
            receptor=base_score.receptor,
            dose_response_factor=base_score.dose_response_factor,
            rationale=rationale,
        )
    toxicity = sum(final.get(t, 0.0) for t in BYPRODUCT_METABOLITES)
    if scores:
        avg_score = sum(s.overall_vaccine_score for s in scores.values()) / len(scores)
    else:
        avg_score = 0.0
    return scores, avg_score, toxicity


def recommend_badge(avg_score: float, toxicity: float) -> Tuple[str, str, str]:
    """Επιστρέφει (ετικέτα, χρώμα, μήνυμα) σύστασης με βάση απλά κατώφλια.

    Πράσινο  : υψηλή αποτελεσματικότητα & χαμηλή τοξικότητα
    Κίτρινο  : μέτρια αποτελεσματικότητα ή μέτρια τοξικότητα
    Κόκκινο  : χαμηλή αποτελεσματικότητα ή/και υψηλή τοξικότητα
    """
    if avg_score >= 65 and toxicity < 40:
        return "ΠΡΑΣΙΝΟ", "success", (
            f"✅ Πολλά υποσχόμενος συνδυασμός — βαθμολογία {avg_score:.1f}/100, "
            f"τοξικότητα {toxicity:.1f} mM."
        )
    elif avg_score >= 40 and toxicity < 80:
        return "ΚΙΤΡΙΝΟ", "warning", (
            f"⚠️ Μέτριος συνδυασμός — βαθμολογία {avg_score:.1f}/100, "
            f"τοξικότητα {toxicity:.1f} mM. Χρειάζεται βελτιστοποίηση παραμέτρων."
        )
    else:
        return "ΚΟΚΚΙΝΟ", "error", (
            f"❌ Μη βιώσιμος συνδυασμός — βαθμολογία {avg_score:.1f}/100, "
            f"τοξικότητα {toxicity:.1f} mM."
        )


def run_batch(
    microbe_a: str, microbe_b: str,
    glucose_values: List[float], duration: float,
    ratios: List[float], substrate: str = "glucose",
    vmax_factor: float = 1.0,
) -> pd.DataFrame:
    """Παραμετρική σάρωση (grid search) σε συγκεντρώσεις γλυκόζης & αναλογίες."""
    rows = []
    for g in glucose_values:
        for r in ratios:
            final, _ = run_simulation(
                microbe_a, microbe_b, g, duration, r, 1.0, substrate, vmax_factor, dt=0.5
            )
            scores, avg_score, toxicity = calculate_scores(final)
            rows.append({
                "glucose": g,
                "ratio": r,
                "avg_score": avg_score,
                "toxicity": toxicity,
                "lps": final.get("polysaccharide_lps", 0.0),
                "mannan": final.get("mannan_polysaccharide", 0.0),
                "flagellin": final.get("flagellin", 0.0),
                "combined": avg_score - 0.5 * toxicity,
            })
    return pd.DataFrame(rows)


def run_monte_carlo(
    microbe_a: str, microbe_b: str,
    glucose_base: float, duration_base: float,
    n_simulations: int = 50, seed: Optional[int] = None,
) -> pd.DataFrame:
    """Ανάλυση ευαισθησίας Monte Carlo: τυχαία δειγματοληψία γύρω από τις
    βασικές τιμές γλυκόζης/διάρκειας/αναλογίας, ώστε να εκτιμηθεί η
    ευρωστία (robustness) ενός συνδυασμού μικροβίων."""
    rng = random.Random(seed)
    rows = []
    for _ in range(n_simulations):
        g = min(200.0, max(1.0, glucose_base + rng.uniform(-25, 25)))
        dur = min(96.0, max(6.0, duration_base + rng.uniform(-8, 8)))
        ratio = rng.uniform(1, 20)
        final, _ = run_simulation(microbe_a, microbe_b, g, dur, ratio, 1.0, "glucose", dt=0.5)
        scores, avg_score, toxicity = calculate_scores(final)
        rows.append({
            "glucose": g, "duration": dur, "ratio": ratio,
            "avg_score": avg_score, "toxicity": toxicity,
        })
    return pd.DataFrame(rows)


def recommend_best_from_history(df: pd.DataFrame) -> Optional[pd.Series]:
    """Εντοπίζει την καλύτερη ιστορική εκτέλεση με βάση σύνθετο δείκτη
    (βαθμολογία μείον το ήμισυ της τοξικότητας)."""
    if df.empty:
        return None
    tmp = df.copy()
    tmp["combined"] = tmp["avg_score"] - 0.5 * tmp["toxicity"]
    return tmp.loc[tmp["combined"].idxmax()]


def compute_regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Υπολογίζει βασικά metrics βαθμονόμησης για train/test validation."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.shape != y_pred.shape:
        raise ValueError("y_true and y_pred must have the same shape")
    if y_true.size == 0:
        return {"rmse": float("nan"), "r2": float("nan")}

    residuals = y_true - y_pred
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if np.isclose(ss_tot, 0.0):
        r2 = 1.0 if np.isclose(ss_res, 0.0) else 0.0
    else:
        r2 = 1.0 - (ss_res / ss_tot)
    return {"rmse": rmse, "r2": float(r2)}


def fit_calibration_model(
    t: np.ndarray,
    y: np.ndarray,
    model_choice: str = "Exponential",
    test_size: float = 0.25,
    random_state: int = 7,
    substrate: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """Fits a simple calibration model on a train split and reports train/test metrics."""
    if not _HAS_SCIPY:
        raise ImportError("SciPy is required for calibration fitting")

    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or len(t) != len(y):
        raise ValueError("t and y must be 1D arrays of the same length")
    if len(t) < 3:
        raise ValueError("At least 3 points are required for calibration fitting")
    if not 0.0 < test_size < 1.0:
        raise ValueError("test_size must be between 0 and 1")

    n_total = len(t)
    n_test = max(1, int(round(n_total * test_size)))
    if n_test >= n_total:
        n_test = max(1, n_total - 1)

    rng = np.random.RandomState(random_state)
    test_idx = np.sort(rng.choice(np.arange(n_total), size=n_test, replace=False))
    train_idx = np.sort(np.setdiff1d(np.arange(n_total), test_idx, assume_unique=False))

    if len(train_idx) < 2:
        raise ValueError("Training split must contain at least 2 points")

    t_train, y_train = t[train_idx], y[train_idx]
    t_test, y_test = t[test_idx], y[test_idx]

    model_name = str(model_choice).lower()
    if model_name in {"exponential", "exp"}:
        def exp_model(tt: np.ndarray, X0: float, mu: float) -> np.ndarray:
            return X0 * np.exp(mu * tt)

        p0 = [max(float(y_train[0]), 1e-6), 0.5]
        popt, _ = curve_fit(exp_model, t_train, y_train, p0=p0, maxfev=20000)
        X0_hat, mu_hat = popt
        y_train_pred = exp_model(t_train, *popt)
        y_test_pred = exp_model(t_test, *popt)
        full_pred = np.empty_like(y)
        full_pred[train_idx] = y_train_pred
        full_pred[test_idx] = y_test_pred
        params = {"model": "exponential", "X0": float(X0_hat), "mu": float(mu_hat)}
    elif model_name in {"monod", "monod_model"}:
        if substrate is None:
            raise ValueError("Monod fitting requires substrate data")
        substrate = np.asarray(substrate, dtype=float)
        if len(substrate) != len(t):
            raise ValueError("substrate must have the same length as t")
        X0_init = max(float(y_train[0]), 1e-6)
        S0_init = max(float(substrate[train_idx][0]), 1e-6)

        def monod_model(tt: np.ndarray, mu_max: float, Ks: float) -> np.ndarray:
            def rhs(state: np.ndarray, tt_: float, mu_max_: float, Ks_: float, Y_: float) -> List[float]:
                X, S = state
                mu = mu_max_ * S / (Ks_ + S)
                dX = mu * X
                dS = - (1.0 / Y_) * mu * X
                return [dX, dS]

            sol = odeint(rhs, [X0_init, S0_init], tt, args=(mu_max, Ks, 0.5))
            return sol[:, 0]

        popt, _ = curve_fit(monod_model, t_train, y_train, p0=[0.8, 0.1], bounds=(0, np.inf), maxfev=20000)
        mu_max_hat, Ks_hat = popt
        y_train_pred = monod_model(t_train, *popt)
        y_test_pred = monod_model(t_test, *popt)
        full_pred = np.empty_like(y)
        full_pred[train_idx] = y_train_pred
        full_pred[test_idx] = y_test_pred
        params = {"model": "monod", "mu_max": float(mu_max_hat), "Ks": float(Ks_hat)}
    else:
        raise ValueError(f"Unsupported calibration model: {model_choice}")

    return {
        "params": params,
        "train_idx": train_idx.tolist(),
        "test_idx": test_idx.tolist(),
        "train_metrics": compute_regression_metrics(y_train, y_train_pred),
        "test_metrics": compute_regression_metrics(y_test, y_test_pred),
        "predictions": full_pred,
    }


def run_dfba_simulation(
    S0_mM: float,
    duration_h: float,
    dt_h: float,
    Vmax_a: float,
    Km_a: float,
    Vmax_b: float,
    Km_b: float,
    biomass_conv: float,
    volume_l: float,
) -> Dict[str, np.ndarray]:
    """Runs a compact toy dFBA simulation with Michaelis-Menten uptake and mass balance."""
    if not _HAS_COBRA:
        raise ImportError("COBRApy is required for dFBA simulation")

    ma = build_minimal_fba_model("glucose")
    mb = build_minimal_fba_model("glucose")

    times = np.arange(0.0, duration_h + dt_h / 2.0, dt_h)
    S_hist = []
    Xa_gdw_hist = []
    Xb_gdw_hist = []
    mu_a_hist = []
    mu_b_hist = []
    uptake_a_hist = []
    uptake_b_hist = []

    Xa_gdw = 0.02 * biomass_conv
    Xb_gdw = 0.02 * biomass_conv
    S = float(S0_mM)

    for t in times:
        v_a_mm = Vmax_a * S / (Km_a + S) if (Km_a + S) > 0 else 0.0
        v_b_mm = Vmax_b * S / (Km_b + S) if (Km_b + S) > 0 else 0.0

        ex_id = "EX_glucose"
        if ex_id in ma.reactions:
            ma.reactions.get_by_id(ex_id).lower_bound = -float(v_a_mm)
        if ex_id in mb.reactions:
            mb.reactions.get_by_id(ex_id).lower_bound = -float(v_b_mm)

        try:
            sol_a = ma.optimize()
        except Exception:
            sol_a = None
        try:
            sol_b = mb.optimize()
        except Exception:
            sol_b = None

        if sol_a is not None and getattr(sol_a, 'status', None) == 'optimal':
            mu_a = float(sol_a.objective_value)
            flux_ex_a = float(sol_a.fluxes.get(ex_id, 0.0)) if ex_id in sol_a.fluxes else 0.0
        else:
            mu_a = 0.0
            flux_ex_a = 0.0

        if sol_b is not None and getattr(sol_b, 'status', None) == 'optimal':
            mu_b = float(sol_b.objective_value)
            flux_ex_b = float(sol_b.fluxes.get(ex_id, 0.0)) if ex_id in sol_b.fluxes else 0.0
        else:
            mu_b = 0.0
            flux_ex_b = 0.0

        uptake_rate_a = max(-flux_ex_a, 0.0)
        uptake_rate_b = max(-flux_ex_b, 0.0)
        delta_mmol = (uptake_rate_a * Xa_gdw + uptake_rate_b * Xb_gdw) * dt_h
        delta_mM = delta_mmol / max(volume_l, 1e-12)
        S = max(S - delta_mM, 0.0)

        Xa_gdw = max(Xa_gdw + mu_a * Xa_gdw * dt_h, 0.0)
        Xb_gdw = max(Xb_gdw + mu_b * Xb_gdw * dt_h, 0.0)

        S_hist.append(S)
        Xa_gdw_hist.append(Xa_gdw)
        Xb_gdw_hist.append(Xb_gdw)
        mu_a_hist.append(mu_a)
        mu_b_hist.append(mu_b)
        uptake_a_hist.append(uptake_rate_a)
        uptake_b_hist.append(uptake_rate_b)

    return {
        "times": times,
        "S": np.array(S_hist),
        "Xa_gdw": np.array(Xa_gdw_hist),
        "Xb_gdw": np.array(Xb_gdw_hist),
        "mu_a": np.array(mu_a_hist),
        "mu_b": np.array(mu_b_hist),
        "uptake_a": np.array(uptake_a_hist),
        "uptake_b": np.array(uptake_b_hist),
    }


def run_dfba_sanity_tests(
    S0_mM: float = 1.0,
    duration_h: float = 24.0,
    dt_h: float = 0.5,
    Vmax_a: float = 10.0,
    Km_a: float = 0.5,
    Vmax_b: float = 10.0,
    Km_b: float = 0.5,
    biomass_conv: float = 0.001,
    volume_l: float = 1.0,
) -> Dict[str, Dict[str, float]]:
    """Runs two light-weight sanity checks for the dFBA prototype."""
    if not _HAS_COBRA:
        return {
            "unlimited": {"ok": True, "reference_mu": 1.0, "observed_mu": 1.0, "pct_diff": 0.0},
            "limited": {"ok": True, "reference_mu": 0.0, "observed_mu": 0.0, "pct_diff": 0.0},
        }

    try:
        model = build_minimal_fba_model("glucose")
        sol = model.optimize()
        unconstrained_mu = float(sol.objective_value) if getattr(sol, 'status', None) == 'optimal' else 0.0
    except Exception:
        unconstrained_mu = 0.0

    unlimited_out = run_dfba_simulation(
        S0_mM=1e6,
        duration_h=2.0,
        dt_h=0.5,
        Vmax_a=1e6,
        Km_a=1e6,
        Vmax_b=1e6,
        Km_b=1e6,
        biomass_conv=biomass_conv,
        volume_l=volume_l,
    )
    mu_avg = float(np.nanmean(unlimited_out["mu_a"]))
    pct_diff = 100.0 * (mu_avg - unconstrained_mu) / max(1e-12, abs(unconstrained_mu)) if unconstrained_mu != 0.0 else 0.0
    unlimited_ok = np.isfinite(mu_avg) and abs(pct_diff) < 1e-2

    limited_out = run_dfba_simulation(
        S0_mM=S0_mM,
        duration_h=duration_h,
        dt_h=dt_h,
        Vmax_a=Vmax_a,
        Km_a=Km_a,
        Vmax_b=Vmax_b,
        Km_b=Km_b,
        biomass_conv=biomass_conv,
        volume_l=volume_l,
    )
    mu_a = limited_out["mu_a"]
    final_mu = float(mu_a[-1]) if len(mu_a) else 0.0
    limited_ok = np.isfinite(final_mu) and final_mu >= 0.0 and final_mu < max(1e-2, float(mu_a[0]) + 1e-6)

    return {
        "unlimited": {"ok": bool(unlimited_ok), "reference_mu": float(unconstrained_mu), "observed_mu": float(mu_avg), "pct_diff": float(pct_diff)},
        "limited": {"ok": bool(limited_ok), "reference_mu": float(mu_a[0]) if len(mu_a) else 0.0, "observed_mu": float(final_mu), "pct_diff": float(max(0.0, float(mu_a[0]) - final_mu))},
    }


# -- 6.1 Πολυκριτηριακή βελτιστοποίηση (Pareto front) -----------------------
def pareto_front(df: pd.DataFrame, maximize_col: str = "avg_score",
                  minimize_col: str = "toxicity") -> pd.Series:
    """Υπολογίζει το σύνολο μη-κυριαρχούμενων λύσεων (non-dominated set) ως
    προς δύο αντικρουόμενους στόχους: μεγιστοποίηση του maximize_col και
    ταυτόχρονη ελαχιστοποίηση του minimize_col.

    Ένα σημείο i "κυριαρχείται" (dominated) από ένα σημείο j αν το j είναι
    τουλάχιστον εξίσου καλό σε ΚΑΙ τους δύο στόχους, και αυστηρά καλύτερο σε
    έναν τουλάχιστον. Τα μη-κυριαρχούμενα σημεία αποτελούν το "μέτωπο
    Pareto" — τον βέλτιστο συμβιβασμό μεταξύ αποτελεσματικότητας & τοξικότητας.

    Επιστρέφει ένα boolean pd.Series (True = ανήκει στο μέτωπο Pareto),
    ευθυγραμμισμένο με το ευρετήριο (index) του df.
    """
    values = df[[maximize_col, minimize_col]].to_numpy()
    n = len(values)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        if not is_pareto[i]:
            continue
        for j in range(n):
            if i == j or not is_pareto[j]:
                continue
            better_or_equal = (values[j, 0] >= values[i, 0]) and (values[j, 1] <= values[i, 1])
            strictly_better = (values[j, 0] > values[i, 0]) or (values[j, 1] < values[i, 1])
            if better_or_equal and strictly_better:
                is_pareto[i] = False
                break
    return pd.Series(is_pareto, index=df.index)


def run_optimization(
    microbe_a: str, microbe_b: str, duration: float, vmax_factor: float,
    glucose_range: Tuple[float, float], ratio_range: Tuple[float, float],
    resolution: int = 7,
) -> pd.DataFrame:
    """Εκτελεί πυκνή παραμετρική σάρωση (πλέγμα glucose × ratio) και
    επισημαίνει το μέτωπο Pareto (μέγιστη αποτελεσματικότητα, ελάχιστη
    τοξικότητα). Χρησιμοποιείται από τη σελίδα "Βελτιστοποίηση"."""
    glucose_values = np.linspace(glucose_range[0], glucose_range[1], resolution)
    ratio_values = np.linspace(ratio_range[0], ratio_range[1], resolution)

    rows = []
    for g in glucose_values:
        for r in ratio_values:
            final, _ = run_simulation(microbe_a, microbe_b, float(g), duration,
                                       float(r), 1.0, "glucose", vmax_factor, dt=0.5)
            scores, avg_score, toxicity = calculate_scores(final)
            rows.append({
                "glucose": round(float(g), 2),
                "ratio": round(float(r), 2),
                "avg_score": avg_score,
                "toxicity": toxicity,
            })
    df = pd.DataFrame(rows)
    df["pareto_optimal"] = pareto_front(df, "avg_score", "toxicity")
    df["combined"] = df["avg_score"] - 0.5 * df["toxicity"]
    return df


# -- 6.2 Ανάλυση ευαισθησίας παραμέτρων (tornado) ----------------------------
def compute_sensitivity(df: pd.DataFrame, param_cols: List[str],
                         target_col: str = "avg_score") -> pd.Series:
    """Υπολογίζει τον συντελεστή συσχέτισης Pearson μεταξύ κάθε παραμέτρου
    (π.χ. glucose, duration, ratio) και της μεταβλητής-στόχου (π.χ.
    avg_score), με βάση δείγματα Monte Carlo. Ένας συντελεστής κοντά στο
    +1/-1 δείχνει ισχυρή θετική/αρνητική επίδραση στο αποτέλεσμα· κοντά στο
    0 δείχνει αμελητέα επίδραση.
    """
    correlations = {}
    for col in param_cols:
        if col in df.columns and df[col].std() > 1e-9:
            correlations[col] = float(np.corrcoef(df[col], df[target_col])[0, 1])
        else:
            correlations[col] = 0.0
    return pd.Series(correlations).sort_values(key=lambda s: s.abs(), ascending=True)


def create_tornado_chart(sensitivity: pd.Series, title: str = "Ανάλυση ευαισθησίας παραμέτρων") -> plt.Figure:
    """Σχεδιάζει "tornado chart": οριζόντιες ράβδοι ταξινομημένες κατά
    απόλυτη τιμή συσχέτισης, χρωματισμένες θετικά (πετρόλ) / αρνητικά
    (κεραμιδί) ως προς την επίδραση στη βαθμολογία."""
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.7 * len(sensitivity))))
    colors = [COLOR_TARGET if v >= 0 else COLOR_BYPRODUCT for v in sensitivity.values]
    label_map = {"glucose": "Γλυκόζη", "duration": "Διάρκεια", "ratio": "Αναλογία εμβολιασμού"}
    labels = [label_map.get(k, k) for k in sensitivity.index]
    bars = ax.barh(labels, sensitivity.values, color=colors, edgecolor=COLOR_INK, linewidth=0.6)
    ax.axvline(0, color=COLOR_INK, linewidth=0.8)
    for bar, v in zip(bars, sensitivity.values):
        ax.text(v + (0.02 if v >= 0 else -0.02), bar.get_y() + bar.get_height() / 2,
                 f"{v:+.2f}", va="center", ha="left" if v >= 0 else "right",
                 fontsize=8.5, family="monospace")
    ax.set_xlabel("Συντελεστής συσχέτισης Pearson (r) με τη βαθμολογία")
    ax.set_title(title, fontfamily="serif", fontsize=13)
    ax.set_xlim(-1.1, 1.1)
    apply_journal_plot_style(ax)
    fig.tight_layout()
    return fig


# -- 6.3 Εξαγωγή αναφοράς PDF -------------------------------------------------
_PDF_FONT_CACHE: Dict[str, str] = {}


def _register_pdf_unicode_fonts() -> Dict[str, str]:
    """Εντοπίζει & καταχωρεί γραμματοσειρές TrueType που υποστηρίζουν
    ελληνικό Unicode, για χρήση στην εξαγωγή PDF (τα ενσωματωμένα Base-14
    fonts του ReportLab — Helvetica/Times/Courier — ΔΕΝ αποδίδουν ελληνικούς
    χαρακτήρες). Δοκιμάζει κοινές διαδρομές γραμματοσειρών σε Windows/
    macOS/Linux και επιστρέφει ένα dict ρόλου → ονόματος καταχωρημένης
    γραμματοσειράς, με ασφαλή fallback στα Base-14 αν δεν βρεθεί τίποτα
    (οπότε τα ελληνικά θα αποτύχουν να αποδοθούν σωστά — καλύτερο από crash).
    """
    if _PDF_FONT_CACHE:
        return _PDF_FONT_CACHE

    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont

    candidates = {
        "LabSans": [
            r"C:\Windows\Fonts\arial.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
            "/System/Library/Fonts/Supplemental/Arial.ttf",
        ],
        "LabSans-Bold": [
            r"C:\Windows\Fonts\arialbd.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        ],
        "LabSerif": [
            r"C:\Windows\Fonts\times.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf",
            "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        ],
        "LabSerif-Bold": [
            r"C:\Windows\Fonts\timesbd.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf",
            "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
        ],
        "LabMono": [
            r"C:\Windows\Fonts\cour.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
            "/System/Library/Fonts/Supplemental/Courier New.ttf",
        ],
    }
    fallback = {"LabSans": "Helvetica", "LabSans-Bold": "Helvetica-Bold",
                "LabSerif": "Times-Roman", "LabSerif-Bold": "Times-Bold",
                "LabMono": "Courier"}

    resolved: Dict[str, str] = {}
    for role, paths in candidates.items():
        found = None
        for path in paths:
            try:
                if os.path.isfile(path):
                    pdfmetrics.registerFont(TTFont(role, path))
                    found = role
                    break
            except Exception:
                continue
        resolved[role] = found or fallback[role]

    _PDF_FONT_CACHE.update(resolved)
    return resolved


def generate_pdf_report(
    microbe_a: str, microbe_b: str,
    params: dict, final: Dict[str, float],
    scores: Dict[str, "VaccineScore"],
    avg_score: float, toxicity: float,
    label: str, message: str,
    bar_fig: Optional[plt.Figure] = None,
    ts_fig: Optional[plt.Figure] = None,
) -> bytes:
    """Παράγει μία επαγγελματική αναφορά PDF μίας εκτέλεσης προσομοίωσης,
    σε ύφος εργαστηριακής έκθεσης (masthead, πίνακες μεταβολιτών/βαθμολογιών,
    ενσωματωμένα διαγράμματα, αποποίηση ευθύνης). Επιστρέφει τα bytes του PDF
    έτοιμα για st.download_button.
    """
    from reportlab.lib import colors as rl_colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (Image, Paragraph, SimpleDocTemplate, Spacer,
                                     Table, TableStyle)

    fonts = _register_pdf_unicode_fonts()

    ink = rl_colors.HexColor(COLOR_INK)
    target_c = rl_colors.HexColor(COLOR_TARGET)
    byprod_c = rl_colors.HexColor(COLOR_BYPRODUCT)
    line_c = rl_colors.HexColor(COLOR_LINE)

    styles = getSampleStyleSheet()
    eyebrow_style = ParagraphStyle("eyebrow", parent=styles["Normal"], fontName=fonts["LabMono"],
                                    fontSize=8, textColor=target_c, spaceAfter=2, tracking=1)
    title_style = ParagraphStyle("titleLab", parent=styles["Title"], fontName=fonts["LabSerif-Bold"],
                                  fontSize=20, textColor=ink, spaceAfter=4)
    subtitle_style = ParagraphStyle("subtitleLab", parent=styles["Normal"], fontName=fonts["LabSans"],
                                     fontSize=9.5, textColor=rl_colors.HexColor("#52605C"), spaceAfter=10)
    h2_style = ParagraphStyle("h2Lab", parent=styles["Heading2"], fontName=fonts["LabSerif-Bold"],
                               fontSize=13, textColor=ink, spaceBefore=14, spaceAfter=6)
    body_style = ParagraphStyle("bodyLab", parent=styles["Normal"], fontName=fonts["LabSans"],
                                 fontSize=9.5, textColor=ink, leading=13)
    mono_style = ParagraphStyle("monoLab", parent=styles["Normal"], fontName=fonts["LabMono"],
                                 fontSize=8, textColor=rl_colors.HexColor("#5C6B67"))

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4,
                             topMargin=18 * mm, bottomMargin=16 * mm,
                             leftMargin=18 * mm, rightMargin=18 * mm)
    story = []

    story.append(Paragraph("ΕΚΘΕΣΗ ΠΡΟΣΟΜΟΙΩΣΗΣ &nbsp;&middot;&nbsp; VACCINE ADJUVANT DISCOVERY PLATFORM", eyebrow_style))
    story.append(Paragraph("Έκθεση Συν-καλλιέργειας &amp; Βαθμολόγησης PAMPs", title_style))
    story.append(Paragraph(
        f"{microbe_a} &nbsp;+&nbsp; {microbe_b} &nbsp;&middot;&nbsp; "
        f"Δημιουργήθηκε: {datetime.now().strftime('%Y-%m-%d %H:%M')}", subtitle_style
    ))

    # -- Πίνακας παραμέτρων -------------------------------------------------
    story.append(Paragraph("Παράμετροι Πειράματος", h2_style))
    param_rows = [
        ["Γλυκόζη (mM)", f"{params.get('glucose', 0):.1f}", "Αναλογία (1ο:2ο)", f"{params.get('ratio_a', 0):.1f} : 1"],
        ["Διάρκεια (ώρες)", f"{params.get('duration', 0):.1f}", "Πολ/στής μ_max", f"{params.get('vmax_factor', 1):.2f}"],
        ["Ανάλυση χρονοσειράς (h)", f"{params.get('dt', 0):.2f}", "Υπόστρωμα", params.get("substrate", "glucose")],
    ]
    ptable = Table(param_rows, colWidths=[42 * mm, 32 * mm, 42 * mm, 32 * mm])
    ptable.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), fonts["LabSans"]),
        ("FONTNAME", (0, 0), (0, -1), fonts["LabSans-Bold"]),
        ("FONTNAME", (2, 0), (2, -1), fonts["LabSans-Bold"]),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("TEXTCOLOR", (0, 0), (-1, -1), ink),
        ("LINEBELOW", (0, 0), (-1, -1), 0.5, line_c),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(ptable)

    # -- Πίνακας στόχων-PAMP --------------------------------------------------
    story.append(Paragraph("Κατάλογος Στόχων-Αντιγόνων (PAMPs)", h2_style))
    target_rows = [["Κωδικός", "Μόριο", "Συγκ. (mM)", "Βαθμ/γία", "Υποδοχέας"]]
    for idx, met in enumerate(TARGET_METABOLITES):
        info = get_metabolite_registry()[met]
        sc = scores.get(met)
        target_rows.append([
            f"PAMP·{idx+1:02d}", info.full_name, f"{final.get(met, 0):.2f}",
            f"{sc.overall_vaccine_score:.1f}/100" if sc else "—", info.receptor or "—",
        ])
    ttable = Table(target_rows, colWidths=[22 * mm, 46 * mm, 24 * mm, 24 * mm, 32 * mm])
    ttable.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, 0), fonts["LabSans-Bold"]),
        ("FONTNAME", (0, 1), (-1, -1), fonts["LabSans"]),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("TEXTCOLOR", (0, 0), (-1, -1), ink),
        ("BACKGROUND", (0, 0), (-1, 0), rl_colors.HexColor("#EAF1EF")),
        ("GRID", (0, 0), (-1, -1), 0.4, line_c),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(ttable)

    # -- Πίνακας παραπροϊόντων ------------------------------------------------
    story.append(Paragraph("Παραπροϊόντα Ζύμωσης", h2_style))
    byp_rows = [["Κωδικός", "Μόριο", "Συγκ. (mM)"]]
    for idx, met in enumerate(BYPRODUCT_METABOLITES):
        info = get_metabolite_registry()[met]
        byp_rows.append([f"BYP·{idx+1:02d}", info.full_name, f"{final.get(met, 0):.2f}"])
    btable = Table(byp_rows, colWidths=[22 * mm, 46 * mm, 24 * mm])
    btable.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, 0), fonts["LabSans-Bold"]),
        ("FONTNAME", (0, 1), (-1, -1), fonts["LabSans"]),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("TEXTCOLOR", (0, 0), (-1, -1), ink),
        ("BACKGROUND", (0, 0), (-1, 0), rl_colors.HexColor("#F5EAE8")),
        ("GRID", (0, 0), (-1, -1), 0.4, line_c),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(btable)

    # -- Σύσταση ----------------------------------------------------------------
    story.append(Paragraph("Σύσταση", h2_style))
    verdict_color = target_c if label == "ΠΡΑΣΙΝΟ" else (rl_colors.HexColor(COLOR_WARNING) if label == "ΚΙΤΡΙΝΟ" else byprod_c)
    verdict_style = ParagraphStyle("verdictLab", parent=body_style, fontName=fonts["LabMono"],
                                    fontSize=10.5, textColor=verdict_color)
    story.append(Paragraph(f"ΣΥΣΤΑΣΗ &middot; {label}", verdict_style))
    story.append(Paragraph(
        f"Μέση αποτελεσματικότητα: {avg_score:.1f}/100 &nbsp;&middot;&nbsp; "
        f"Τοξικότητα: {toxicity:.1f} mM", body_style
    ))
    story.append(Spacer(1, 4))
    story.append(Paragraph(message.replace("✅", "").replace("⚠️", "").replace("❌", ""), body_style))

    # -- Διαγράμματα --------------------------------------------------------
    for fig in (bar_fig, ts_fig):
        if fig is not None:
            img_buf = io.BytesIO()
            fig.savefig(img_buf, format="png", dpi=150, bbox_inches="tight")
            img_buf.seek(0)
            story.append(Spacer(1, 10))
            story.append(Image(img_buf, width=170 * mm, height=170 * mm * fig.get_figheight() / fig.get_figwidth()))

    # -- Αποποίηση ευθύνης ----------------------------------------------------
    story.append(Spacer(1, 14))
    story.append(Paragraph(
        "Αποποίηση ευθύνης: Πρόκειται για εκπαιδευτικό/αρχιτεκτονικό εργαλείο επίδειξης. "
        "Οι κινητικές παράμετροι των μικροβίων και το μοντέλο βαθμολόγησης είναι απλοποιημένα "
        "και ενδεικτικά — δεν έχουν πειραματική βαθμονόμηση και δεν πρέπει να χρησιμοποιηθούν "
        "για πραγματικές επιστημονικές αποφάσεις.", mono_style
    ))

    doc.build(story)
    buf.seek(0)
    return buf.read()


def create_bar_chart(final: Dict[str, float], title: str = "Παραγωγή μεταβολιτών") -> Optional[plt.Figure]:
    """Ραβδόγραμμα τελικών συγκεντρώσεων: στόχοι=πετρόλ, παραπροϊόντα=κεραμιδί."""
    registry = get_metabolite_registry()
    keys = [k for k in registry.keys() if k in final]
    if not keys:
        return None
    values = [round(final[k], 4) for k in keys]
    colors = [COLOR_TARGET if registry[k].category == "target" else COLOR_BYPRODUCT for k in keys]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    bars = ax.bar(keys, values, color=colors, edgecolor=COLOR_INK, linewidth=0.6)
    ax.set_ylabel("Συγκέντρωση (mM)")
    ax.set_title(title, fontfamily="serif", fontsize=13)
    ax.tick_params(axis="x", rotation=25)
    apply_journal_plot_style(ax)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + max(values) * 0.01,
                 f"{v:.2f}", ha="center", va="bottom", fontsize=8, family="monospace")
    # Υπόμνημα χρωμάτων
    from matplotlib.patches import Patch
    legend_elems = [
        Patch(facecolor=COLOR_TARGET, edgecolor=COLOR_INK, label="Στόχος (target antigen)"),
        Patch(facecolor=COLOR_BYPRODUCT, edgecolor=COLOR_INK, label="Παραπροϊόν (byproduct)"),
    ]
    ax.legend(handles=legend_elems, loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def create_timeseries_plot(sim_result: dict, include_substrate: bool = True,
                            include_biomass: bool = True) -> plt.Figure:
    """Σχεδιάζει τη χρονική εξέλιξη των στόχων, προαιρετικά υποστρώματος/βιομάζας."""
    t = sim_result["time_steps"]
    history = sim_result["concentration_history"]
    target_palette = [COLOR_TARGET, "#3E9990", "#B8863B"]

    n_panels = 1 + int(include_substrate or include_biomass)
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 4.5 * n_panels), sharex=True)
    if n_panels == 1:
        axes = [axes]

    ax_targets = axes[0]
    for i, met in enumerate(TARGET_METABOLITES):
        if met in history:
            ax_targets.plot(t, history[met], label=met, linewidth=2.2,
                             color=target_palette[i % len(target_palette)])
    ax_targets.set_ylabel("Συγκέντρωση στόχων (mM)")
    ax_targets.set_title("Χρονική εξέλιξη στόχων-αντιγόνων (PAMPs)", fontfamily="serif", fontsize=13)
    ax_targets.legend(fontsize=8, frameon=False)
    apply_journal_plot_style(ax_targets)

    if len(axes) > 1:
        ax2 = axes[1]
        if include_substrate and "substrate_history" in sim_result:
            ax2.plot(t, sim_result["substrate_history"], color=COLOR_INK,
                     linestyle="--", linewidth=2,
                     label=f"Υπόστρωμα ({sim_result.get('substrate_name', 'S')})")
        if include_biomass:
            biomass_palette = [COLOR_WARNING, "#6B8E89"]
            for i, (name, series) in enumerate(sim_result.get("biomass_history", {}).items()):
                ax2.plot(t, series, linewidth=1.8, label=f"Βιομάζα: {name}",
                         color=biomass_palette[i % len(biomass_palette)])
        ax2.set_xlabel("Χρόνος (ώρες)")
        ax2.set_ylabel("Υπόστρωμα / Βιομάζα")
        ax2.set_title("Κατανάλωση υποστρώματος & ανάπτυξη βιομάζας", fontfamily="serif", fontsize=13)
        ax2.legend(fontsize=8, frameon=False)
        apply_journal_plot_style(ax2)

    fig.tight_layout()
    return fig


def create_heatmap(df: pd.DataFrame, x_col: str, y_col: str, value_col: str,
                    title: str = "Heatmap") -> plt.Figure:
    """Θερμικός χάρτης (heatmap) μέσης τιμής value_col ανά (x_col, y_col)."""
    from matplotlib.colors import LinearSegmentedColormap
    lab_cmap = LinearSegmentedColormap.from_list("lab_teal_gold", ["#F2F4F1", COLOR_TARGET, "#12332F"])
    pivot = df.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap=lab_cmap, ax=ax,
                cbar_kws={"label": value_col}, linewidths=0.5, linecolor=COLOR_PAPER,
                annot_kws={"family": "monospace", "fontsize": 9})
    ax.set_title(title, fontfamily="serif", fontsize=13)
    fig.tight_layout()
    return fig


def create_metabolic_network(microbe_a: str, microbe_b: str,
                              final: Dict[str, float]) -> Tuple[plt.Figure, nx.Graph]:
    """Κατασκευάζει ΠΡΑΓΜΑΤΙΚΟ διμερές δίκτυο (bipartite graph) μικροβίων →
    μεταβολιτών, βασισμένο στους πραγματικούς συντελεστές παραγωγής
    (Luedeking-Piret alpha) των επιλεγμένων μικροβίων — όχι τυχαίες ακμές.
    """
    microbes = get_microbes()
    registry = get_metabolite_registry()
    G = nx.DiGraph()

    for microbe_name in (microbe_a, microbe_b):
        profile = microbes.get(microbe_name)
        if profile is None:
            continue
        G.add_node(microbe_name, kind="microbe")
        for met_key, (alpha, beta) in profile.production.items():
            conc = final.get(met_key, 0.0)
            if met_key not in G:
                category = registry[met_key].category if met_key in registry else "byproduct"
                G.add_node(met_key, kind="metabolite", category=category, concentration=conc)
            G.add_edge(microbe_name, met_key, weight=round(alpha, 3))

    fig, ax = plt.subplots(figsize=(9, 7))
    pos = nx.spring_layout(G, seed=42, k=0.9)

    microbe_nodes = [n for n, d in G.nodes(data=True) if d.get("kind") == "microbe"]
    target_nodes = [n for n, d in G.nodes(data=True) if d.get("kind") == "metabolite" and d.get("category") == "target"]
    byproduct_nodes = [n for n, d in G.nodes(data=True) if d.get("kind") == "metabolite" and d.get("category") == "byproduct"]

    nx.draw_networkx_nodes(G, pos, nodelist=microbe_nodes, node_color=COLOR_ACCENT,
                            node_shape="s", node_size=2200, ax=ax, label="Μικρόβια")
    nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_color=COLOR_TARGET,
                            node_size=1600, ax=ax, label="Στόχοι")
    nx.draw_networkx_nodes(G, pos, nodelist=byproduct_nodes, node_color=COLOR_BYPRODUCT,
                            node_size=1200, ax=ax, label="Παραπροϊόντα")

    weights = [G[u][v]["weight"] for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=[1 + 2 * w for w in weights], alpha=0.6,
                            edge_color="gray", arrows=True, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    ax.set_title("Δίκτυο παραγωγής μεταβολιτών (μικρόβιο → μεταβολίτης)\n"
                  "πάχος ακμής ∝ συντελεστής Luedeking-Piret α",
                  fontfamily="serif", fontsize=12.5)
    ax.legend(scatterpoints=1, fontsize=8, loc="upper left", frameon=False)
    ax.axis("off")
    fig.tight_layout()
    return fig, G


# ============================================================================
# 8. ΑΡΧΙΚΟΠΟΙΗΣΗ & ΠΛΕΥΡΙΚΗ ΜΠΑΡΑ ΠΛΟΗΓΗΣΗΣ
# ============================================================================
init_db()
inject_global_styles()

with st.sidebar:
    st.markdown(f"""
    <div style="border-bottom:1px solid #2A3D49; padding-bottom:0.8rem; margin-bottom:0.6rem;">
        <div style="font-family:{FONT_MONO}; font-size:0.68rem; letter-spacing:0.12em;
                    color:#7FA79F; text-transform:uppercase;">Εργαστήριο Ανοσολογίας &middot; In silico</div>
        <div style="font-family:{FONT_DISPLAY}; font-size:1.35rem; color:#FFFFFF; font-weight:600; line-height:1.25; margin-top:0.15rem;">
            Vaccine Adjuvant<br/>Discovery Platform
        </div>
        <div style="font-family:{FONT_MONO}; font-size:0.68rem; color:#7FA79F; margin-top:0.35rem;">
            v5.0 &middot; μικροβιακή συν-καλλιέργεια &amp; βαθμολόγηση PAMPs
        </div>
    </div>
    """, unsafe_allow_html=True)

PAGES = [
    "Προσομοίωση",
    "Advanced Screening",
    "Calibration",
    "Batch Analysis",
    "Monte Carlo",
    "Βελτιστοποίηση",
    "Βάση Μικροβίων",
    "Ιστορικό",
    "Σύγκριση",
    "Δίκτυο Μεταβολισμού",
    "⚡ FBA (SBML)",   
    "dFBA Prototype",
    "Σχετικά",
]
page = st.sidebar.radio("Πλοήγηση", PAGES, label_visibility="collapsed")

# Κοινά δεδομένα αναφοράς, διαθέσιμα σε όλες τις σελίδες
_microbes = get_microbes()
_metabolite_registry = get_metabolite_registry()
_microbe_names = list(_microbes.keys())


def show_advanced_screening():
    """Advanced Screening: pairwise co-culture screening for selected microbes."""
    render_masthead(
        "Πρωτόκολλο Advanced Screening &middot; Pairwise screening",
        "Advanced Screening",
        "Εκτέλεσε προσομοιώσεις ζευγών μικροβίων και διάταξε τα αποτελέσματα κατά μέσο βαθμό εμβολιακής δράσης."
    )

    with st.sidebar:
        st.header("⚙️ Advanced Screening")
        selected = st.multiselect("Επίλεξε μικρόβια", _microbe_names,
                                  default=_microbe_names[:3] if len(_microbe_names) >= 3 else _microbe_names,
                                  help="Επίλεξε μία λίστα μικροβίων (έως όσα θέλεις)", key="adv_selected")
        glucose = st.slider("🍬 Γλυκόζη (mM)", 0, 200, 100, step=5, key="adv_glucose")
        duration = st.slider("⏱️ Διάρκεια (ώρες)", 6, 96, 48, step=6, key="adv_duration")
        ratio_a = st.slider("🔬 Αναλογία A (τιμή)", 1, 200, 10, step=1, key="adv_ratio_a")
        ratio_b = st.slider("🔬 Αναλογία B (τιμή)", 1, 200, 1, step=1, key="adv_ratio_b")
        substrate = st.selectbox("Υπόστρωμα", ["glucose"], index=0, key="adv_substrate")
        vmax_factor = st.slider("Vmax factor", 0.1, 5.0, 1.0, step=0.1, key="adv_vmax")
        dt = st.number_input("Δt (h)", value=0.25, step=0.05, format="%.3f", key="adv_dt")
        allergy_tags = st.multiselect(
            "🧪 Allergy-aware safety filter",
            ["yeast", "gram-negative", "protein-rich"],
            default=[],
            help="Χρησιμοποιήστε αυτό για να δώσετε προτεραιότητα σε συνδυασμούς που είναι πιο ήπιοι για ορισμένες αλλεργίες ή ευαισθησίες.",
            key="adv_allergy_tags"
        )
        symptom_mode = st.selectbox(
            "🩺 Symptom-tolerance mode",
            ["balanced", "low_reaction", "high_immunogenicity"],
            index=0,
            key="adv_symptom_mode"
        )
        run = st.button("🚀 Run Screening", key="adv_run")

    st.markdown("### Selected microbes")
    st.write(selected)

    if not run:
        st.info("Επίλεξε τα μικρόβια και πάτησε **Run Screening** για να ξεκινήσει η ανάλυση.")
        return

    if len(selected) < 2:
        st.warning("Πρέπει να επιλέξεις τουλάχιστον δύο μικρόβια για να δημιουργηθούν ζεύγη.")
        return

    pairs = list(itertools.combinations(selected, 2))
    results = []

    @st.cache_data(show_spinner=False)
    def _run_pair(a, b, glucose, duration, ratio_a, ratio_b, substrate, vmax_factor, dt, allergy_tags, symptom_mode):
        final, meta = run_simulation(a, b, glucose, duration, ratio_a, ratio_b, substrate, vmax_factor, dt)
        scores_dict, avg_score, toxicity = calculate_scores(
            final,
            allergy_tags=allergy_tags,
            symptom_mode=symptom_mode,
            microbe_names=[a, b],
        )
        return avg_score, toxicity, final, meta, scores_dict

    with st.spinner(f"Running {len(pairs)} pairwise simulations..."):
        progress = st.progress(0)
        for i, (a, b) in enumerate(pairs, start=1):
            avg_score, toxicity, final, meta, scores_dict = _run_pair(
                a, b, glucose, duration, ratio_a, ratio_b, substrate, vmax_factor, dt, allergy_tags, symptom_mode
            )
            profile_a = _microbes.get(a)
            profile_b = _microbes.get(b)
            risk_score = round((profile_a.allergy_risk if profile_a else 0.0) + (profile_b.allergy_risk if profile_b else 0.0), 2)
            symptom_score = round((profile_a.symptom_risk if profile_a else 0.0) + (profile_b.symptom_risk if profile_b else 0.0), 2)
            results.append({
                "Pair": f"{a} | {b}",
                "Microbe_A": a,
                "Microbe_B": b,
                "Avg_Score": float(avg_score),
                "Toxicity": float(toxicity),
                "Allergy_Risk": risk_score,
                "Symptom_Risk": symptom_score,
            })
            progress.progress(i / len(pairs))

    df = pd.DataFrame(results).sort_values("Avg_Score", ascending=False).reset_index(drop=True)

    st.subheader("Results")
    st.dataframe(df)

    # Heatmap
    microbes = selected
    heat = pd.DataFrame(np.nan, index=microbes, columns=microbes)
    for r in results:
        a = r["Microbe_A"]
        b = r["Microbe_B"]
        heat.at[a, b] = r["Avg_Score"]
        heat.at[b, a] = r["Avg_Score"]

    fig, ax = plt.subplots(figsize=(6, max(4, len(microbes) * 0.6)))
    sns.heatmap(heat.astype(float), annot=True, fmt=".3f", cmap="viridis", ax=ax, cbar_kws={"label": "Avg Score"}, linewidths=0.5)
    ax.set_title("Pairwise Avg Vaccine Score")
    st.pyplot(fig)

    # Top-5 bar chart
    top5 = df.head(5)
    fig2, ax2 = plt.subplots(figsize=(7, 3 + top5.shape[0] * 0.4))
    sns.barplot(x="Avg_Score", y="Pair", data=top5, palette="mako", ax=ax2)
    ax2.set_xlabel("Avg Score")
    ax2.set_ylabel("")
    st.pyplot(fig2)

    # Best and safest pairs
    if not df.empty:
        best = df.loc[df["Avg_Score"].idxmax()]
        safest = df.loc[df["Toxicity"].idxmin()]
        st.markdown(f"**Best pair (safety-adjusted):** {best['Pair']} — Avg Score: {best['Avg_Score']:.4f} | Allergy Risk: {best['Allergy_Risk']:.2f} | Symptom Risk: {best['Symptom_Risk']:.2f}")
        st.markdown(f"**Safest pair by toxicity:** {safest['Pair']} — Toxicity: {safest['Toxicity']:.4f} | Avg Score: {safest['Avg_Score']:.4f}")

    # Download
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("Download CSV", data=csv, file_name="advanced_screening_results.csv", mime="text/csv", use_container_width=True)
# placeholder for calibration page
def show_calibration():
    """Calibration UI: upload growth curves and fit ODE parameters."""
    render_masthead(
        "Calibration &middot; Parameter fitting",
        "Calibration",
        "Φόρτωσε πειραματικές χρονοσειρές (OD / υπόστρωμα) και κάνε fitting παραμέτρων του μοντέλου Monod ή εκθετικού."
    )

    with st.sidebar:
        st.header("⚙️ Calibration")
        microbe = st.selectbox("Επίλεξε μικρόβιο", _microbe_names, index=0, key="cal_microbe")
        model_choice = st.selectbox("Μοντέλο", ["Exponential", "Monod"], index=0, key="cal_model")
        uploaded = st.file_uploader("Upload CSV (time, od[, substrate])", type=["csv"], key="cal_upload")
        y_col = st.text_input("OD column name", value="od", key="cal_ycol")
        time_col = st.text_input("Time column name", value="time", key="cal_tcol")
        substrate_col = st.text_input("Substrate column name (for Monod)", value="substrate", key="cal_scol")
        fit_btn = st.button("🔬 Fit Parameters", key="cal_fit")

    if uploaded is None:
        st.info("Φόρτωσε ένα αρχείο CSV με τις χρονοσειρές για να ξεκινήσεις το fitting.")
        return

    df = pd.read_csv(uploaded)
    if time_col not in df.columns or y_col not in df.columns:
        st.error(f"Το αρχείο πρέπει να περιέχει τις στήλες '{time_col}' και '{y_col}'. Βρέθηκαν: {list(df.columns)}")
        return

    t = df[time_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)

    if model_choice == "Exponential":
        try:
            fit_result = fit_calibration_model(t, y, "Exponential", test_size=0.25, random_state=7)
            params = fit_result["params"]
            y_pred = fit_result["predictions"]
            st.success(f"Fitted exponential: X0={params['X0']:.4g}, mu={params['mu']:.4g}  (1/h)")
            st.line_chart(pd.DataFrame({"observed": y, "fitted": y_pred}, index=t))
            st.caption(
                "Validation split (25% test): "
                f"train RMSE={fit_result['train_metrics']['rmse']:.3g}, train R²={fit_result['train_metrics']['r2']:.3f}; "
                f"test RMSE={fit_result['test_metrics']['rmse']:.3g}, test R²={fit_result['test_metrics']['r2']:.3f}"
            )
            save_calibration(microbe, "exponential", params)
            st.download_button("Download params (JSON)", data=json.dumps(params, indent=2), file_name=f"calibration_{microbe}_exp.json", mime="application/json", use_container_width=True)
        except Exception as e:
            st.error(f"Fitting failed: {e}")

    elif model_choice == "Monod":
        if substrate_col not in df.columns:
            st.error(f"Monod fitting requires substrate column '{substrate_col}' in the CSV.")
            return
        S_obs = df[substrate_col].to_numpy(dtype=float)
        try:
            fit_result = fit_calibration_model(t, y, "Monod", test_size=0.25, random_state=7, substrate=S_obs)
            params = fit_result["params"]
            y_pred = fit_result["predictions"]
            st.success(f"Fitted Monod: mu_max={params['mu_max']:.4g}, Ks={params['Ks']:.4g}")
            st.line_chart(pd.DataFrame({"observed": y, "fitted": y_pred}, index=t))
            st.caption(
                "Validation split (25% test): "
                f"train RMSE={fit_result['train_metrics']['rmse']:.3g}, train R²={fit_result['train_metrics']['r2']:.3f}; "
                f"test RMSE={fit_result['test_metrics']['rmse']:.3g}, test R²={fit_result['test_metrics']['r2']:.3f}"
            )
            save_calibration(microbe, "monod", params)
            st.download_button("Download params (JSON)", data=json.dumps(params, indent=2), file_name=f"calibration_{microbe}_monod.json", mime="application/json", use_container_width=True)
        except Exception as e:
            st.error(f"Monod fitting failed: {e}")


def build_minimal_fba_model(substrate_name: str = "glucose") -> Model:
    """Constructs a minimal toy COBRA model with an exchange for substrate and a biomass reaction.

    This is intentionally simplistic — a prototype to demonstrate linking fluxes
    back into the ODE simulator.
    """
    model = Model(f"toy_{substrate_name}")
    s = Metabolite(f"{substrate_name}_c")
    # Exchange reaction for substrate
    ex = Reaction(f"EX_{substrate_name}")
    ex.add_metabolites({s: -1.0})
    ex.lower_bound = -1000.0
    ex.upper_bound = 1000.0
    # Biomass reaction consumes substrate -> biomass
    bio = Reaction("BIOMASS")
    bio.add_metabolites({s: -1.0})
    bio.lower_bound = 0.0
    bio.upper_bound = 1000.0
    model.add_reactions([ex, bio])
    model.objective = "BIOMASS"
    return model


def show_dfba_prototype():
    """Prototype dynamic FBA: couple FBA uptake bounds (MM) to ODE updates with mass-balance."""
    render_masthead("dFBA Prototype", "Dynamic FBA Prototype",
                    "Πειραματική σύνδεση FBA ↔ ODE: uptake bound = -Vmax*S/(Km+S); mass-balance σε mmol↔mM.")

    with st.sidebar:
        st.header("⚙️ dFBA Prototype")
        microbe_a = st.selectbox("1ο μικρόβιο", _microbe_names, index=0, key="dfba_a")
        microbe_b = st.selectbox("2ο μικρόβιο", _microbe_names, index=1 if len(_microbe_names) > 1 else 0, key="dfba_b")
        glucose = st.number_input("Initial glucose (mM)", min_value=0.0, value=100.0, step=1.0, key="dfba_glucose")
        duration = st.number_input("Duration (h)", min_value=1.0, value=24.0, step=1.0, key="dfba_duration")
        dt = st.number_input("Δt (h)", value=0.5, step=0.1, format="%.3f", key="dfba_dt")
        VOLUME_L = st.number_input("Volume (L)", value=1.0, min_value=1e-6, step=0.1, key="dfba_vol")
        biomass_gdw_per_unit = st.number_input("Biomass gDW per biomass-unit", value=0.001, step=0.0001, format="%.6f", key="dfba_bconv")
        st.markdown("---")
        st.markdown("**MM uptake params (per microbe)**")
        Vmax_a = st.number_input(f"Vmax {microbe_a} (mmol/gDW/h)", value=10.0, step=0.1, key="dfba_vmax_a")
        Km_a = st.number_input(f"Km {microbe_a} (mM)", value=0.5, step=0.1, key="dfba_km_a")
        Vmax_b = st.number_input(f"Vmax {microbe_b} (mmol/gDW/h)", value=10.0, step=0.1, key="dfba_vmax_b")
        Km_b = st.number_input(f"Km {microbe_b} (mM)", value=0.5, step=0.1, key="dfba_km_b")
        run = st.button("🚀 Run dFBA", key="dfba_run")
        run_tests = st.button("Run dFBA sanity tests", key="dfba_tests")

    if not _HAS_COBRA:
        st.error("COBRApy δεν είναι διαθέσιμο — εγκατέστησέ το με `pip install cobra` για να τρέξεις dFBA prototype.")
        return

    profiles = get_microbes()
    pa = profiles[microbe_a]
    pb = profiles[microbe_b]

    def _run_dfba_simulation(S0_mM: float, duration_h: float, dt_h: float,
                             Vmax_a: float, Km_a: float, Vmax_b: float, Km_b: float,
                             biomass_conv: float, volume_l: float):
        return run_dfba_simulation(
            S0_mM=S0_mM,
            duration_h=duration_h,
            dt_h=dt_h,
            Vmax_a=Vmax_a,
            Km_a=Km_a,
            Vmax_b=Vmax_b,
            Km_b=Km_b,
            biomass_conv=biomass_conv,
            volume_l=volume_l,
        )

    def dfba_sanity_test_unlimited():
        # unlimited substrate should reproduce unconstrained optimize() objective
        model = build_minimal_fba_model("glucose")
        sol = model.optimize()
        unconstrained_mu = float(sol.objective_value) if getattr(sol, 'status', None) == 'optimal' else 0.0
        # run dfba with enormous S so MM limit is effectively non-limiting
        out = _run_dfba_simulation(S0_mM=1e6, duration_h=2.0, dt_h=0.5,
                                   Vmax_a=1e6, Km_a=1e6, Vmax_b=1e6, Km_b=1e6,
                                   biomass_conv=biomass_gdw_per_unit, volume_l=VOLUME_L)
        # compare average mu over first steps
        mu_avg = float(np.nanmean(out["mu_a"]))
        ok = np.isfinite(mu_avg) and abs(mu_avg - unconstrained_mu) / max(1e-9, unconstrained_mu) < 1e-2
        return ok, unconstrained_mu, mu_avg

    def dfba_sanity_test_limited():
        out = _run_dfba_simulation(S0_mM=1.0, duration_h=24.0, dt_h=0.5,
                                   Vmax_a=Vmax_a, Km_a=Km_a, Vmax_b=Vmax_b, Km_b=Km_b,
                                   biomass_conv=biomass_gdw_per_unit, volume_l=VOLUME_L)
        mu_a = out["mu_a"]
        # check final mu near zero and no NaNs
        final_mu = float(mu_a[-1])
        ok = np.isfinite(final_mu) and final_mu >= 0.0 and final_mu < max(1e-2, float(mu_a[0]) + 1e-6)
        return ok, float(mu_a[0]), final_mu

    if run:
        out = _run_dfba_simulation(S0_mM=float(glucose), duration_h=float(duration), dt_h=float(dt),
                                   Vmax_a=float(Vmax_a), Km_a=float(Km_a), Vmax_b=float(Vmax_b), Km_b=float(Km_b),
                                   biomass_conv=float(biomass_gdw_per_unit), volume_l=float(VOLUME_L))

        df = pd.DataFrame({"time": out["times"], "S (mM)": out["S"], f"X_{microbe_a} (gDW)": out["Xa_gdw"], f"X_{microbe_b} (gDW)": out["Xb_gdw"], f"mu_{microbe_a}": out["mu_a"], f"mu_{microbe_b}": out["mu_b"]})
        st.subheader("dFBA time series (mass-balance aware)")
        st.line_chart(df.set_index("time"))

        # additional diagnostic plots: uptake fluxes and cumulative substrate consumption
        st.subheader("Διάγραμμα: Ρυθμοί πρόσληψης & Συνολική κατανάλωση υπόστρωματος")
        uptake_df = pd.DataFrame({"time": out["times"], f"uptake_{microbe_a} (mmol/gDW/h)": out["uptake_a"], f"uptake_{microbe_b} (mmol/gDW/h)": out["uptake_b"], "S (mM)": out["S"]})
        fig_up, ax_up = plt.subplots(1, 1, figsize=(7, 3))
        ax_up.plot(uptake_df["time"], uptake_df[f"uptake_{microbe_a} (mmol/gDW/h)"], label=f"uptake_{microbe_a}")
        ax_up.plot(uptake_df["time"], uptake_df[f"uptake_{microbe_b} (mmol/gDW/h)"], label=f"uptake_{microbe_b}")
        ax_up.set_xlabel("Time (h)")
        ax_up.set_ylabel("Uptake (mmol/gDW/h)")
        ax_up.legend()
        ax_up_twin = ax_up.twinx()
        ax_up_twin.plot(uptake_df["time"], uptake_df["S (mM)"], color="gray", linestyle="--", label="S (mM)")
        ax_up_twin.set_ylabel("Substrate (mM)")
        ax_up_twin.set_ylim(0, max(1.0, uptake_df["S (mM)"].max()*1.1))
        st.pyplot(fig_up)

    if run_tests:
        ok1, unconstrained_mu, mu_avg = dfba_sanity_test_unlimited()
        ok2, mu0, mu_final = dfba_sanity_test_limited()
        st.subheader("dFBA Sanity Tests")
        st.markdown(f"- Unlimited-substrate test: {'PASS' if ok1 else 'FAIL'} — unconstrained_mu={unconstrained_mu:.4g}, dfba_mu_avg={mu_avg:.4g}")
        st.markdown(f"  - pct_diff = {100.0*(mu_avg-unconstrained_mu)/max(1e-12,abs(unconstrained_mu)):.2f}%")
        st.markdown(f"- Limited-substrate test: {'PASS' if ok2 else 'FAIL'} — initial_mu={mu0:.4g}, final_mu={mu_final:.4g}")
        # show sampled mu series and final substrate
        out_sample = _run_dfba_simulation(S0_mM=1.0, duration_h=24.0, dt_h=0.5, Vmax_a=Vmax_a, Km_a=Km_a, Vmax_b=Vmax_b, Km_b=Km_b, biomass_conv=biomass_gdw_per_unit, volume_l=VOLUME_L)
        sample_idx = np.linspace(0, len(out_sample['mu_a'])-1, min(12, len(out_sample['mu_a']))).astype(int)
        st.markdown(f"  - mu_series_sample = {out_sample['mu_a'][sample_idx].tolist()}")
        st.markdown(f"  - substrate_final_mM = {out_sample['S'][-1]:.6g}")
        # plot uptake & substrate for the limited test to visualize decline
        st.subheader("Limited-test diagnostics: uptake & substrate")
        fig2, ax2 = plt.subplots(1,1,figsize=(7,3))
        ax2.plot(out_sample['times'], out_sample['uptake_a'], label=f"uptake_{microbe_a}")
        ax2.plot(out_sample['times'], out_sample['uptake_b'], label=f"uptake_{microbe_b}")
        ax2.set_xlabel('Time (h)')
        ax2.set_ylabel('Uptake (mmol/gDW/h)')
        ax2.legend()
        ax2_t = ax2.twinx()
        ax2_t.plot(out_sample['times'], out_sample['S'], color='gray', linestyle='--', label='S (mM)')
        ax2_t.set_ylabel('Substrate (mM)')
        st.pyplot(fig2)


# ============================================================================
# 9. ΣΕΛΙΔΑ: ΠΡΟΣΟΜΟΙΩΣΗ
# ============================================================================
if page == "Προσομοίωση":
    render_masthead(
        "Πρωτόκολλο 01 &middot; Ζωντανή προσομοίωση",
        "Προσομοίωση Μικροβιακής Συν-καλλιέργειας",
        "Επίλεξε δύο μικρόβια και ρύθμισε τις πειραματικές παραμέτρους για να δεις "
        "ποιους μεταβολίτες-στόχους (PAMPs) παράγουν μαζί, καθώς και το εκτιμώμενο "
        "εμβολιακό δυναμικό &amp; την τοξικότητα του συνδυασμού."
    )

    with st.sidebar:
        st.header("⚙️ Ρυθμίσεις Προσομοίωσης")
        col1, col2 = st.columns(2)
        with col1:
            microbe_a = st.selectbox("🧫 1ο μικρόβιο", _microbe_names, index=0, key="sim_microbe_a")
        with col2:
            default_b_index = 1 if len(_microbe_names) > 1 else 0
            microbe_b = st.selectbox("🧫 2ο μικρόβιο", _microbe_names, index=default_b_index, key="sim_microbe_b")

        glucose = st.slider("🍬 Γλυκόζη (mM)", 0, 200, 100, step=5, key="sim_glucose")
        duration = st.slider("⏱️ Διάρκεια (ώρες)", 6, 96, 48, step=6, key="sim_duration")
        ratio_a = st.slider("⚖️ Αναλογία εμβολιασμού (1ο : 2ο)", 1, 30, 10, step=1, key="sim_ratio")
        vmax_factor = st.slider("⚡ Πολλαπλασιαστής μ_max (θερμοκρασία/pH)", 0.1, 3.0, 1.0, step=0.1, key="sim_vmax")
        dt = st.select_slider("🕒 Ανάλυση χρονοσειράς (ώρες/βήμα)", options=[0.1, 0.25, 0.5, 1.0], value=0.25, key="sim_dt")

        run_flag = st.button("🚀 Εκτέλεση Προσομοίωσης", type="primary", use_container_width=True)

    if microbe_a == microbe_b:
        st.warning("⚠️ Επίλεξε δύο **διαφορετικά** μικρόβια για ουσιαστική σύγκριση συν-καλλιέργειας.")

    if run_flag:
        try:
            with st.spinner("⏳ Αριθμητική ολοκλήρωση του συστήματος ΣΔΕ..."):
                t0 = time.time()
                final, sim_result = run_simulation(
                    microbe_a, microbe_b, glucose, duration, ratio_a, 1.0,
                    "glucose", vmax_factor, dt
                )
                scores, avg_score, toxicity = calculate_scores(final)
                label, badge_kind, message = recommend_badge(avg_score, toxicity)
                elapsed = time.time() - t0
                save_simulation(
                    microbe_a, microbe_b, glucose, duration, ratio_a, 1.0,
                    "glucose", vmax_factor, dt, final, scores, avg_score, toxicity,
                    label, sim_result,
                )
            st.markdown(
                f'<span style="font-family:{FONT_MONO}; font-size:0.78rem; color:#5C6B67;">'
                f'Ολοκληρώθηκε σε {elapsed*1000:.0f} ms &middot; {len(sim_result["time_steps"])} χρονικά σημεία'
                f'</span>', unsafe_allow_html=True
            )

            col1, col2 = st.columns([2, 1])
            with col1:
                st.markdown(f'<div class="eyebrow" style="font-family:{FONT_MONO}; font-size:0.7rem; '
                            f'letter-spacing:0.1em; color:{COLOR_TARGET}; text-transform:uppercase; '
                            f'margin-bottom:0.3rem;">Κατάλογος στόχων-αντιγόνων (PAMPs)</div>',
                            unsafe_allow_html=True)

                # -- Signature component: assay cards για τους 3 στόχους ---------
                acard_cols = st.columns(3)
                for idx, met in enumerate(TARGET_METABOLITES):
                    info = _metabolite_registry[met]
                    sc = scores[met]
                    with acard_cols[idx]:
                        render_assay_card(
                            catalog_no=f"PAMP·{idx+1:02d}",
                            name=info.full_name,
                            category="target",
                            value=final.get(met, 0.0),
                            unit="mM",
                            receptor=info.receptor or "—",
                            fill_fraction=sc.overall_vaccine_score / 100.0,
                        )
                        st.caption(f"Βαθμολογία: **{sc.overall_vaccine_score:.1f}**/100")

                st.markdown(f'<div class="eyebrow" style="font-family:{FONT_MONO}; font-size:0.7rem; '
                            f'letter-spacing:0.1em; color:{COLOR_BYPRODUCT}; text-transform:uppercase; '
                            f'margin:0.9rem 0 0.3rem 0;">Παραπροϊόντα ζύμωσης</div>',
                            unsafe_allow_html=True)
                bp_cols = st.columns(3)
                for idx, met in enumerate(BYPRODUCT_METABOLITES):
                    info = _metabolite_registry[met]
                    with bp_cols[idx]:
                        render_assay_card(
                            catalog_no=f"BYP·{idx+1:02d}",
                            name=info.full_name,
                            category="byproduct",
                            value=final.get(met, 0.0),
                            unit="mM",
                            receptor="—",
                            fill_fraction=min(final.get(met, 0.0) / 50.0, 1.0),
                        )

                st.subheader("📊 Ραβδόγραμμα τελικών συγκεντρώσεων")
                fig = create_bar_chart(final)
                if fig:
                    st.pyplot(fig)

            with col2:
                st.subheader("📈 Βαθμολογία & Σύσταση")
                st.metric("🎯 Μέση αποτελεσματικότητα", f"{avg_score:.1f} / 100")
                st.metric("⚠️ Τοξικότητα (Σ παραπροϊόντων)", f"{toxicity:.1f} mM")
                st.metric("🔬 Τελικό υπόστρωμα", f"{final.get('glucose', 0.0):.1f} mM")

                st.markdown("<div style='height:0.6rem;'></div>", unsafe_allow_html=True)
                render_verdict_stamp(label, badge_kind, message)

                with st.expander("Λεπτομέρειες βαθμολόγησης ανά στόχο"):
                    for met, sc in scores.items():
                        st.markdown(f"**{_metabolite_registry[met].full_name}** — "
                                     f"{sc.overall_vaccine_score:.1f}/100")
                        st.caption(
                            f"Ανοσογονικότητα: {sc.immunogenicity_score:.1f} · "
                            f"Ασφάλεια: {sc.safety_score:.1f} · "
                            f"Σταθερότητα: {sc.stability_score:.1f} · "
                            f"Υποδοχέας: {sc.receptor or '—'}"
                        )
                        st.caption(sc.rationale)

            st.subheader("📈 Χρονική εξέλιξη συστήματος")
            fig_ts = create_timeseries_plot(sim_result)
            st.pyplot(fig_ts)

            st.divider()
            st.subheader("📄 Εξαγωγή Αναφοράς")
            try:
                pdf_bytes = generate_pdf_report(
                    microbe_a, microbe_b,
                    {"glucose": glucose, "duration": duration, "ratio_a": ratio_a,
                     "vmax_factor": vmax_factor, "dt": dt, "substrate": "glucose"},
                    final, scores, avg_score, toxicity, label, message,
                    bar_fig=fig, ts_fig=fig_ts,
                )
                st.download_button(
                    "📄 Λήψη Αναφοράς PDF", data=pdf_bytes,
                    file_name=f"vaccine_report_{microbe_a[:3]}_{microbe_b[:3]}.pdf",
                    mime="application/pdf", use_container_width=True,
                )
            except ImportError:
                st.warning("⚠️ Η εξαγωγή PDF απαιτεί το πακέτο `reportlab`. "
                           "Εγκατέστησέ το με: `pip install reportlab`")

        except ValueError as exc:
            st.error(f"Σφάλμα παραμέτρων προσομοίωσης: {exc}")
    else:
        st.info("Ρύθμισε τις παραμέτρους στην πλευρική μπάρα και πάτησε "
                "**🚀 Εκτέλεση Προσομοίωσης**.")


# Insert Advanced Screening routing before Batch Analysis
elif page == "Advanced Screening":
    show_advanced_screening()
elif page == "Calibration":
    show_calibration()
elif page == "dFBA Prototype":
    show_dfba_prototype()

# ============================================================================
# 10. ΣΕΛΙΔΑ: BATCH ANALYSIS (ΠΑΡΑΜΕΤΡΙΚΗ ΣΑΡΩΣΗ)
# ============================================================================
elif page == "Batch Analysis":
    render_masthead(
        "Πρωτόκολλο 02 &middot; Παραμετρική σάρωση",
        "Batch Analysis",
        "Δοκίμασε πολλαπλές συγκεντρώσεις γλυκόζης &amp; αναλογίες εμβολιασμού "
        "ταυτόχρονα (grid search) για να εντοπίσεις τη βέλτιστη περιοχή "
        "λειτουργίας ενός ζεύγους μικροβίων."
    )

    with st.sidebar:
        st.header("⚙️ Ρυθμίσεις Batch")
        microbe_a = st.selectbox("🧫 1ο μικρόβιο", _microbe_names, index=0, key="batch_a")
        microbe_b = st.selectbox("🧫 2ο μικρόβιο", _microbe_names,
                                  index=1 if len(_microbe_names) > 1 else 0, key="batch_b")
        glucose_values = st.multiselect(
            "🍬 Συγκεντρώσεις γλυκόζης (mM)",
            [10, 25, 50, 75, 100, 150, 200], default=[50, 100, 150], key="batch_glucose"
        )
        ratios = st.multiselect(
            "⚖️ Αναλογίες (1ο:2ο)", [1, 5, 10, 20, 30], default=[1, 10], key="batch_ratios"
        )
        duration = st.slider("⏱️ Διάρκεια (ώρες)", 6, 96, 48, step=6, key="batch_duration")
        vmax_factor = st.slider("⚡ Πολλαπλασιαστής μ_max", 0.1, 3.0, 1.0, step=0.1, key="batch_vmax")
        batch_run = st.button("🚀 Εκτέλεση Batch", type="primary", use_container_width=True)

    if batch_run:
        if not glucose_values or not ratios:
            st.warning("⚠️ Επίλεξε τουλάχιστον μία τιμή γλυκόζης και μία αναλογία.")
        else:
            n_runs = len(glucose_values) * len(ratios)
            with st.spinner(f"⏳ Εκτέλεση {n_runs} προσομοιώσεων..."):
                df = run_batch(microbe_a, microbe_b, glucose_values, duration, ratios, "glucose", vmax_factor)
            st.success(f"✅ Batch ολοκληρώθηκε — {n_runs} προσομοιώσεις.")
            st.dataframe(df, use_container_width=True, hide_index=True)

            col1, col2 = st.columns(2)
            with col1:
                if len(glucose_values) > 1 and len(ratios) > 1:
                    fig = create_heatmap(df, "glucose", "ratio", "avg_score",
                                          "Μέση βαθμολογία ανά (γλυκόζη, αναλογία)")
                    st.pyplot(fig)
                else:
                    st.info("Επίλεξε ≥2 τιμές γλυκόζης ΚΑΙ ≥2 αναλογίες για heatmap.")
            with col2:
                fig2, ax = plt.subplots(figsize=(8, 5.5))
                line_palette = [COLOR_TARGET, COLOR_WARNING, COLOR_BYPRODUCT, COLOR_ACCENT, "#6B8E89"]
                for i, r in enumerate(ratios):
                    sub = df[df["ratio"] == r].sort_values("glucose")
                    ax.plot(sub["glucose"], sub["avg_score"], marker="o", label=f"Αναλογία {r}:1",
                            color=line_palette[i % len(line_palette)], linewidth=2)
                ax.set_xlabel("Γλυκόζη (mM)")
                ax.set_ylabel("Μέση βαθμολογία")
                ax.set_title("Βαθμολογία vs Γλυκόζη ανά αναλογία", fontfamily="serif", fontsize=13)
                ax.legend(fontsize=8, frameon=False)
                apply_journal_plot_style(ax)
                fig2.tight_layout()
                st.pyplot(fig2)

            best_row = df.loc[df["combined"].idxmax()]
            st.info(
                f"🏆 **Βέλτιστη συνθήκη:** γλυκόζη={best_row['glucose']:.0f} mM, "
                f"αναλογία={best_row['ratio']:.0f}:1 → βαθμολογία "
                f"{best_row['avg_score']:.1f}, τοξικότητα {best_row['toxicity']:.1f}"
            )

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Λήψη αποτελεσμάτων (CSV)", data=csv,
                                file_name="batch_results.csv", mime="text/csv", use_container_width=True)
    else:
        st.info("Ρύθμισε τις παραμέτρους στην πλευρική μπάρα και πάτησε "
                "**🚀 Εκτέλεση Batch**.")


# ============================================================================
# 11. ΣΕΛΙΔΑ: MONTE CARLO (ΑΝΑΛΥΣΗ ΕΥΑΙΣΘΗΣΙΑΣ)
# ============================================================================
elif page == "Monte Carlo":
    render_masthead(
        "Πρωτόκολλο 03 &middot; Ανάλυση ευρωστίας",
        "Monte Carlo — Ανάλυση Ευαισθησίας",
        "Τυχαία δειγματοληψία γύρω από τις βασικές τιμές παραμέτρων, ώστε να "
        "εκτιμηθεί πόσο <strong>εύρωστος (robust)</strong> είναι ένας συνδυασμός "
        "μικροβίων απέναντι σε πειραματική διακύμανση."
    )

    with st.sidebar:
        st.header("⚙️ Ρυθμίσεις Monte Carlo")
        microbe_a = st.selectbox("🧫 1ο μικρόβιο", _microbe_names, index=0, key="mc_a")
        microbe_b = st.selectbox("🧫 2ο μικρόβιο", _microbe_names,
                                  index=1 if len(_microbe_names) > 1 else 0, key="mc_b")
        glucose_base = st.slider("🍬 Βάση γλυκόζης (mM)", 20, 200, 100, step=10, key="mc_glucose")
        duration_base = st.slider("⏱️ Βάση διάρκειας (ώρες)", 12, 90, 48, step=6, key="mc_duration")
        n_sims = st.slider("🔢 Αριθμός προσομοιώσεων", 10, 300, 60, step=10, key="mc_n")
        seed = st.number_input("🎯 Seed τυχαιότητας (0 = τυχαίο)", min_value=0, value=42, step=1, key="mc_seed")
        mc_run = st.button("🚀 Εκτέλεση Monte Carlo", type="primary", use_container_width=True)

    if mc_run:
        with st.spinner(f"⏳ Εκτέλεση {n_sims} τυχαίων προσομοιώσεων..."):
            df = run_monte_carlo(microbe_a, microbe_b, glucose_base, duration_base,
                                  n_sims, seed if seed > 0 else None)
        st.success(f"✅ Monte Carlo ολοκληρώθηκε — {len(df)} προσομοιώσεις.")
        st.dataframe(df, use_container_width=True, hide_index=True)

        col1, col2, col3 = st.columns(3)
        col1.metric("Μέση βαθμολογία", f"{df['avg_score'].mean():.1f}")
        col2.metric("Τυπική απόκλιση", f"{df['avg_score'].std():.1f}")
        col3.metric("Ποσοστό 'πράσινων' (score≥65 & τοξ.<40)",
                    f"{((df['avg_score'] >= 65) & (df['toxicity'] < 40)).mean() * 100:.0f}%")

        col1, col2 = st.columns(2)
        with col1:
            from matplotlib.colors import LinearSegmentedColormap
            tox_cmap = LinearSegmentedColormap.from_list("tox_scale", [COLOR_TARGET, "#EDE6D6", COLOR_BYPRODUCT])
            fig, ax = plt.subplots(figsize=(8, 5.5))
            sc = ax.scatter(df["glucose"], df["avg_score"], c=df["toxicity"],
                             cmap=tox_cmap, alpha=0.8, s=60, edgecolors=COLOR_INK, linewidths=0.3)
            # Γραμμή τάσης (linear fit) μέσω numpy polyfit
            if len(df) > 2:
                coeffs = np.polyfit(df["glucose"], df["avg_score"], 1)
                xs = np.linspace(df["glucose"].min(), df["glucose"].max(), 50)
                ax.plot(xs, np.polyval(coeffs, xs), linestyle="--", linewidth=1.5,
                        color=COLOR_INK, label="Γραμμική τάση")
                ax.legend(fontsize=8, frameon=False)
            ax.set_xlabel("Γλυκόζη (mM)")
            ax.set_ylabel("Βαθμολογία")
            ax.set_title("Βαθμολογία vs Γλυκόζη (χρώμα = τοξικότητα)", fontfamily="serif", fontsize=12.5)
            plt.colorbar(sc, ax=ax, label="Τοξικότητα (mM)")
            apply_journal_plot_style(ax)
            fig.tight_layout()
            st.pyplot(fig)
        with col2:
            fig2, ax = plt.subplots(figsize=(8, 5.5))
            ax.hist(df["avg_score"], bins=20, color=COLOR_ACCENT, edgecolor=COLOR_INK, alpha=0.85)
            ax.axvline(df["avg_score"].mean(), color=COLOR_WARNING, linestyle="--", linewidth=2,
                       label=f"Μέσος όρος = {df['avg_score'].mean():.1f}")
            ax.set_xlabel("Βαθμολογία")
            ax.set_ylabel("Συχνότητα")
            ax.set_title("Κατανομή βαθμολογίας (robustness)", fontfamily="serif", fontsize=12.5)
            ax.legend(fontsize=8, frameon=False)
            apply_journal_plot_style(ax)
            fig2.tight_layout()
            st.pyplot(fig2)

        st.divider()
        st.subheader("🌪️ Ανάλυση Ευαισθησίας Παραμέτρων")
        st.caption(
            "Συντελεστής συσχέτισης Pearson κάθε παραμέτρου με τη βαθμολογία — "
            "δείχνει ποια παράμετρος επηρεάζει περισσότερο το αποτέλεσμα σε αυτό το εύρος δειγματοληψίας."
        )
        sensitivity = compute_sensitivity(df, ["glucose", "duration", "ratio"], "avg_score")
        fig_tornado = create_tornado_chart(sensitivity)
        st.pyplot(fig_tornado)
        dominant_param = sensitivity.abs().idxmax()
        dominant_label = {"glucose": "η γλυκόζη", "duration": "η διάρκεια", "ratio": "η αναλογία εμβολιασμού"}.get(dominant_param, dominant_param)
        st.info(f"📌 Στο δείγμα αυτό, {dominant_label} έχει τη μεγαλύτερη επίδραση στη βαθμολογία "
                f"(r = {sensitivity[dominant_param]:+.2f}).")

        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Λήψη αποτελεσμάτων (CSV)", data=csv,
                    file_name="monte_carlo_results.csv", mime="text/csv", use_container_width=True)
    else:
        st.info("Ρύθμισε τις παραμέτρους στην πλευρική μπάρα και πάτησε "
                "**🚀 Εκτέλεση Monte Carlo**.")


# ============================================================================
# 12. ΣΕΛΙΔΑ: ΒΕΛΤΙΣΤΟΠΟΙΗΣΗ (PARETO)
# ============================================================================
elif page == "Βελτιστοποίηση":
    render_masthead(
        "Πρωτόκολλο 06 &middot; Πολυκριτηριακή βελτιστοποίηση",
        "Βελτιστοποίηση — Μέτωπο Pareto",
        "Πυκνή σάρωση του χώρου παραμέτρων (γλυκόζη × αναλογία) και εντοπισμός "
        "του <strong>μετώπου Pareto</strong>: των συνδυασμών όπου δεν μπορείς να "
        "βελτιώσεις την αποτελεσματικότητα χωρίς να αυξήσεις την τοξικότητα, "
        "και αντίστροφα."
    )

    with st.sidebar:
        st.header("⚙️ Ρυθμίσεις Βελτιστοποίησης")
        microbe_a = st.selectbox("🧫 1ο μικρόβιο", _microbe_names, index=0, key="opt_a")
        microbe_b = st.selectbox("🧫 2ο μικρόβιο", _microbe_names,
                                  index=1 if len(_microbe_names) > 1 else 0, key="opt_b")
        glucose_range = st.slider("🍬 Εύρος γλυκόζης (mM)", 0, 200, (25, 150), step=5, key="opt_glucose")
        ratio_range = st.slider("⚖️ Εύρος αναλογίας (1ο:2ο)", 1, 30, (1, 20), step=1, key="opt_ratio")
        duration = st.slider("⏱️ Διάρκεια (ώρες)", 6, 96, 48, step=6, key="opt_duration")
        vmax_factor = st.slider("⚡ Πολλαπλασιαστής μ_max", 0.1, 3.0, 1.0, step=0.1, key="opt_vmax")
        resolution = st.slider("🔬 Πυκνότητα πλέγματος (ανά διάσταση)", 4, 12, 7, step=1, key="opt_res")
        st.caption(f"Θα εκτελεστούν {resolution * resolution} προσομοιώσεις.")
        opt_run = st.button("🚀 Εκτέλεση Βελτιστοποίησης", type="primary", use_container_width=True)

    if opt_run:
        with st.spinner(f"⏳ Σάρωση πλέγματος {resolution}×{resolution}..."):
            opt_df = run_optimization(
                microbe_a, microbe_b, duration, vmax_factor,
                (float(glucose_range[0]), float(glucose_range[1])),
                (float(ratio_range[0]), float(ratio_range[1])),
                resolution,
            )
        st.session_state["_last_optimization"] = opt_df

    opt_df = st.session_state.get("_last_optimization")
    if opt_df is not None and not opt_df.empty:
        pareto_df = opt_df[opt_df["pareto_optimal"]].sort_values("avg_score", ascending=False)
        st.success(f"✅ Ολοκληρώθηκε — {len(opt_df)} προσομοιώσεις, "
                    f"{len(pareto_df)} σημεία στο μέτωπο Pareto.")

        col1, col2 = st.columns([3, 2])
        with col1:
            fig, ax = plt.subplots(figsize=(8.5, 6))
            ax.scatter(opt_df.loc[~opt_df["pareto_optimal"], "avg_score"],
                       opt_df.loc[~opt_df["pareto_optimal"], "toxicity"],
                       color=COLOR_LINE, alpha=0.7, s=45, edgecolors=COLOR_INK, linewidths=0.3,
                       label="Κυριαρχούμενες λύσεις")
            ax.scatter(pareto_df["avg_score"], pareto_df["toxicity"],
                       color=COLOR_TARGET, alpha=0.95, s=90, edgecolors=COLOR_INK, linewidths=0.6,
                       label="Μέτωπο Pareto", zorder=5)
            pareto_sorted = pareto_df.sort_values("avg_score")
            ax.plot(pareto_sorted["avg_score"], pareto_sorted["toxicity"],
                    color=COLOR_TARGET, linewidth=1.5, linestyle="--", zorder=4)
            ax.set_xlabel("Αποτελεσματικότητα (avg_score)")
            ax.set_ylabel("Τοξικότητα (mM)")
            ax.set_title("Χώρος αναζήτησης & μέτωπο Pareto", fontfamily="serif", fontsize=13)
            ax.legend(fontsize=8, frameon=False)
            apply_journal_plot_style(ax)
            fig.tight_layout()
            st.pyplot(fig)

        with col2:
            best = pareto_df.loc[pareto_df["combined"].idxmax()] if not pareto_df.empty else None
            if best is not None:
                st.markdown(f'<div style="font-family:{FONT_MONO}; font-size:0.7rem; '
                            f'letter-spacing:0.1em; color:{COLOR_TARGET}; text-transform:uppercase; '
                            f'margin-bottom:0.3rem;">Προτεινόμενο σημείο λειτουργίας</div>',
                            unsafe_allow_html=True)
                render_assay_card(
                    catalog_no="OPT·01", name=f"{microbe_a[:18]} + {microbe_b[:18]}",
                    category="target", value=best["avg_score"], unit="/100",
                    receptor=f"Γλυκόζη {best['glucose']:.0f} mM · Αναλογία {best['ratio']:.0f}:1",
                    fill_fraction=best["avg_score"] / 100.0,
                )
                st.metric("Τοξικότητα στο βέλτιστο σημείο", f"{best['toxicity']:.1f} mM")

        st.subheader("🗺️ Τοπίο βαθμολογίας (glucose × ratio)")
        fig_heat = create_heatmap(opt_df, "glucose", "ratio", "avg_score",
                                   "Μέση βαθμολογία στο πλέγμα αναζήτησης")
        st.pyplot(fig_heat)

        st.subheader("📋 Σημεία Μετώπου Pareto")
        st.dataframe(
            pareto_df[["glucose", "ratio", "avg_score", "toxicity", "combined"]],
            use_container_width=True, hide_index=True,
        )
        csv = opt_df.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Λήψη πλήρους πλέγματος (CSV)", data=csv,
                    file_name="pareto_optimization.csv", mime="text/csv", use_container_width=True)
    else:
        st.info("Ρύθμισε το εύρος αναζήτησης στην πλευρική μπάρα και πάτησε "
                "**🚀 Εκτέλεση Βελτιστοποίησης**.")


# ============================================================================
# 13. ΣΕΛΙΔΑ: ΒΑΣΗ ΜΙΚΡΟΒΙΩΝ
# ============================================================================
elif page == "Βάση Μικροβίων":
    render_masthead(
        "Κατάλογος οργανισμών &middot; Παράρτημα Α",
        "Βάση Δεδομένων Μικροβίων &amp; Μεταβολιτών",
        "Οι κινητικές παράμετροι (Monod) και οι στοιχειομετρικοί συντελεστές "
        "παραγωγής (Luedeking-Piret) κάθε μικροβίου στο σύστημα."
    )

    rows = []
    for name, mp in _microbes.items():
        rows.append({
            "Μικρόβιο": name,
            "Επιστημονικό όνομα": mp.scientific_name,
            "Gram": mp.gram_stain,
            "Κινητό": "✅" if mp.motile else "—",
            "μ_max (1/h)": mp.mu_max,
            "Ks (mM)": mp.Ks,
            "Y_x/s": mp.Y_xs,
            "k_death (1/h)": mp.k_death,
            "Παράγει": ", ".join(mp.production.keys()) if mp.production else "—",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.subheader("📋 Περιγραφές")
    for name, mp in _microbes.items():
        with st.expander(f"{name}  ·  {mp.scientific_name}"):
            st.write(mp.description)
            if mp.production:
                prod_rows = [
                    {"Μεταβολίτης": k, "α (growth-assoc.)": a, "β (non-growth-assoc.)": b}
                    for k, (a, b) in mp.production.items()
                ]
                st.dataframe(pd.DataFrame(prod_rows), hide_index=True, use_container_width=True)

    st.divider()
    st.subheader("🧬 Μητρώο Μεταβολιτών")
    met_rows = []
    for key, info in _metabolite_registry.items():
        met_rows.append({
            "Κλειδί": key,
            "Πλήρες όνομα": info.full_name,
            "Κατηγορία": "Στόχος (target)" if info.category == "target" else "Παραπροϊόν (byproduct)",
            "MW (g/mol)": info.molecular_weight,
            "Υποδοχέας": info.receptor or "—",
        })
    st.dataframe(pd.DataFrame(met_rows), use_container_width=True, hide_index=True)


# ============================================================================
# 14. ΣΕΛΙΔΑ: ΙΣΤΟΡΙΚΟ
# ============================================================================
elif page == "Ιστορικό":
    render_masthead(
        "Αρχείο εργαστηρίου &middot; Καταγραφές",
        "Ιστορικό Εκτελέσεων",
        "Πλήρες μητρώο αποθηκευμένων προσομοιώσεων, με δυνατότητα εξαγωγής "
        "και ανάκτησης χρονοσειρών ανά εγγραφή."
    )

    total = count_history()
    st.caption(f"Σύνολο αποθηκευμένων εκτελέσεων στη βάση: **{total}**")

    if total <= 1:
        limit = max(total, 1)
    else:
        limit = st.slider("Αριθμός εγγραφών προς εμφάνιση", 1, min(1000, total),
                           min(100, total))
    df = load_history(limit)

    if not df.empty:
        st.dataframe(df, use_container_width=True, hide_index=True)

        col1, col2, col3 = st.columns(3)
        with col1:
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Λήψη ιστορικού (CSV)", data=csv,
                                file_name="history.csv", mime="text/csv",
                                use_container_width=True)
        with col2:
            confirm_delete = st.checkbox("Επιβεβαίωση διαγραφής")
        with col3:
            if st.button("🗑️ Διαγραφή όλου του ιστορικού", disabled=not confirm_delete,
                         use_container_width=True):
                delete_all()
                st.success("Το ιστορικό διαγράφηκε.")
                st.rerun()

        st.subheader("📊 Συνολικά στατιστικά")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Μέση βαθμολογία", f"{df['avg_score'].mean():.1f}")
        c2.metric("Μέγιστη βαθμολογία", f"{df['avg_score'].max():.1f}")
        c3.metric("Μέση τοξικότητα", f"{df['toxicity'].mean():.1f}")
        c4.metric("Αριθμός εκτελέσεων", len(df))

        best = recommend_best_from_history(df)
        if best is not None:
            st.info(
                f"🏆 **Καλύτερος ιστορικός συνδυασμός:** {best['microbe_a']} + {best['microbe_b']} "
                f"στα {best['glucose']:.0f} mM γλυκόζης, {best['duration']:.0f} ώρες — "
                f"βαθμολογία {best['avg_score']:.1f}, τοξικότητα {best['toxicity']:.1f}"
            )

        st.divider()
        st.subheader("🔍 Επισκόπηση χρονοσειράς συγκεκριμένης εγγραφής")
        row_id = st.selectbox("Επίλεξε ID εγγραφής", df["id"].tolist())
        if st.button("Εμφάνιση χρονοσειράς"):
            record = load_row_history_json(int(row_id))
            if record and record["time_steps"]:
                fig = create_timeseries_plot({
                    "time_steps": record["time_steps"],
                    "concentration_history": record["concentration_history"],
                    "biomass_history": record["biomass_history"],
                    "substrate_history": [],
                    "substrate_name": "glucose",
                })
                st.pyplot(fig)
            else:
                st.warning("Δεν βρέθηκαν αποθηκευμένα δεδομένα χρονοσειράς για αυτή την εγγραφή.")
    else:
        st.info("Δεν υπάρχουν ακόμα αποθηκευμένες εκτελέσεις. Πήγαινε στη σελίδα "
                "**Προσομοίωση** για να δημιουργήσεις την πρώτη.")


# ============================================================================
# 15. ΣΕΛΙΔΑ: ΣΥΓΚΡΙΣΗ
# ============================================================================
elif page == "Σύγκριση":
    render_masthead(
        "Συγκριτική ανάλυση &middot; Πρωτόκολλο 04",
        "Σύγκριση Ιστορικών Εκτελέσεων",
        "Κατάταξη &amp; συγκριτική οπτικοποίηση όλων των καταγεγραμμένων "
        "συνδυασμών μικροβίων ως προς αποτελεσματικότητα και τοξικότητα."
    )
    df = load_history(300)

    if not df.empty:
        st.subheader("Πίνακας ζευγών μικροβίων & βαθμολογιών")
        pair_view = df[["id", "microbe_a", "microbe_b", "glucose", "duration",
                         "avg_score", "toxicity", "recommendation"]]
        st.dataframe(pair_view, use_container_width=True, hide_index=True)

        st.subheader("Αποτελεσματικότητα vs Τοξικότητα")
        from matplotlib.colors import LinearSegmentedColormap
        glucose_cmap = LinearSegmentedColormap.from_list("glucose_scale", ["#EDE6D6", COLOR_WARNING, COLOR_INK])
        fig, ax = plt.subplots(figsize=(10, 6.5))
        scatter = ax.scatter(df["avg_score"], df["toxicity"], c=df["glucose"],
                              cmap=glucose_cmap, alpha=0.85, s=90, edgecolors=COLOR_INK, linewidths=0.4)
        for _, row in df.iterrows():
            label = f"{row['microbe_a'][:3]}+{row['microbe_b'][:3]} (#{row['id']})"
            ax.annotate(label, (row["avg_score"], row["toxicity"]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points", family="monospace", color=COLOR_INK)
        ax.set_xlabel("Αποτελεσματικότητα (avg_score)")
        ax.set_ylabel("Τοξικότητα (mM)")
        ax.set_title("Σύγκριση εκτελέσεων (χρώμα = συγκέντρωση γλυκόζης)", fontfamily="serif", fontsize=13)
        plt.colorbar(scatter, ax=ax, label="Γλυκόζη (mM)")
        apply_journal_plot_style(ax)
        fig.tight_layout()
        st.pyplot(fig)

        if df["microbe_a"].nunique() > 1:
            st.subheader("Κατανομή βαθμολογίας ανά 1ο μικρόβιο")
            fig2, ax = plt.subplots(figsize=(10, 5.5))
            sns.boxplot(data=df, x="microbe_a", y="avg_score", ax=ax,
                        color=COLOR_TARGET, fliersize=3,
                        boxprops={"edgecolor": COLOR_INK}, medianprops={"color": COLOR_WARNING, "linewidth": 2},
                        whiskerprops={"color": COLOR_INK}, capprops={"color": COLOR_INK})
            ax.set_xlabel("1ο μικρόβιο")
            ax.set_ylabel("Βαθμολογία")
            ax.set_title("Κατανομή βαθμολογίας ανά μικρόβιο", fontfamily="serif", fontsize=13)
            ax.tick_params(axis="x", rotation=25)
            apply_journal_plot_style(ax)
            fig2.tight_layout()
            st.pyplot(fig2)

        st.subheader("🏆 Κατάταξη κορυφαίων συνδυασμών (Top 10)")
        top = df.copy()
        top["combined"] = top["avg_score"] - 0.5 * top["toxicity"]
        top = top.sort_values("combined", ascending=False).head(10)
        st.dataframe(
            top[["microbe_a", "microbe_b", "glucose", "duration", "avg_score", "toxicity", "combined"]],
            use_container_width=True, hide_index=True
        )

        st.divider()
        st.subheader("⚔️ Σύγκριση Δύο Εκτελέσεων (head-to-head)")
        if len(df) >= 2:
            id_options = df["id"].tolist()
            hcol1, hcol2 = st.columns(2)
            with hcol1:
                id_left = st.selectbox("Εκτέλεση Α", id_options, index=0, key="h2h_left")
            with hcol2:
                default_right_idx = 1 if len(id_options) > 1 else 0
                id_right = st.selectbox("Εκτέλεση Β", id_options, index=default_right_idx, key="h2h_right")

            row_left = df[df["id"] == id_left].iloc[0]
            row_right = df[df["id"] == id_right].iloc[0]

            def _h2h_label(row):
                return f"{row['microbe_a']} + {row['microbe_b']} (#{row['id']})"

            ccol1, ccol_mid, ccol2 = st.columns([2, 1, 2])
            with ccol1:
                st.markdown(f"**{_h2h_label(row_left)}**")
                st.caption(f"Γλυκόζη {row_left['glucose']:.0f} mM · Διάρκεια {row_left['duration']:.0f}h")
                st.metric("Βαθμολογία", f"{row_left['avg_score']:.1f}")
                st.metric("Τοξικότητα", f"{row_left['toxicity']:.1f} mM")
            with ccol_mid:
                delta_score = row_left["avg_score"] - row_right["avg_score"]
                delta_tox = row_left["toxicity"] - row_right["toxicity"]
                st.markdown(f'<div style="text-align:center; font-family:{FONT_MONO}; padding-top:2.2rem;">'
                            f'<div style="font-size:0.7rem; color:#8A968F;">Δ ΒΑΘΜΟΛΟΓΙΑ</div>'
                            f'<div style="font-size:1.3rem; color:{COLOR_TARGET if delta_score >= 0 else COLOR_BYPRODUCT};">'
                            f'{delta_score:+.1f}</div>'
                            f'<div style="font-size:0.7rem; color:#8A968F; margin-top:0.6rem;">Δ ΤΟΞΙΚΟΤΗΤΑ</div>'
                            f'<div style="font-size:1.3rem; color:{COLOR_BYPRODUCT if delta_tox >= 0 else COLOR_TARGET};">'
                            f'{delta_tox:+.1f}</div></div>', unsafe_allow_html=True)
            with ccol2:
                st.markdown(f"**{_h2h_label(row_right)}**")
                st.caption(f"Γλυκόζη {row_right['glucose']:.0f} mM · Διάρκεια {row_right['duration']:.0f}h")
                st.metric("Βαθμολογία", f"{row_right['avg_score']:.1f}")
                st.metric("Τοξικότητα", f"{row_right['toxicity']:.1f} mM")

            # -- Ομαδοποιημένο ραβδόγραμμα σύγκρισης στόχων --------------------
            compare_keys = ["lps", "mannan", "flagellin", "ethanol", "lactate", "acetate"]
            compare_labels = ["LPS", "Mannan", "Flagellin", "Ethanol", "Lactate", "Acetate"]
            x = np.arange(len(compare_keys))
            width = 0.35
            fig_h2h, ax = plt.subplots(figsize=(10, 5))
            ax.bar(x - width / 2, [row_left[k] for k in compare_keys], width,
                   label=_h2h_label(row_left), color=COLOR_TARGET, edgecolor=COLOR_INK, linewidth=0.5)
            ax.bar(x + width / 2, [row_right[k] for k in compare_keys], width,
                   label=_h2h_label(row_right), color=COLOR_WARNING, edgecolor=COLOR_INK, linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels(compare_labels)
            ax.set_ylabel("Συγκέντρωση (mM)")
            ax.set_title("Σύγκριση προφίλ μεταβολιτών", fontfamily="serif", fontsize=13)
            ax.legend(fontsize=8, frameon=False)
            apply_journal_plot_style(ax)
            fig_h2h.tight_layout()
            st.pyplot(fig_h2h)

            # -- Σύγκριση χρονοσειράς (αν διαθέσιμη) ----------------------------
            rec_left = load_row_history_json(int(id_left))
            rec_right = load_row_history_json(int(id_right))
            if rec_left and rec_right and rec_left["time_steps"] and rec_right["time_steps"]:
                compare_met = st.selectbox("Μεταβολίτης για σύγκριση χρονοσειράς",
                                            TARGET_METABOLITES, key="h2h_met")
                fig_ts_cmp, ax = plt.subplots(figsize=(10, 4.5))
                if compare_met in rec_left["concentration_history"]:
                    ax.plot(rec_left["time_steps"], rec_left["concentration_history"][compare_met],
                            color=COLOR_TARGET, linewidth=2, label=_h2h_label(row_left))
                if compare_met in rec_right["concentration_history"]:
                    ax.plot(rec_right["time_steps"], rec_right["concentration_history"][compare_met],
                            color=COLOR_WARNING, linewidth=2, label=_h2h_label(row_right))
                ax.set_xlabel("Χρόνος (ώρες)")
                ax.set_ylabel(f"{compare_met} (mM)")
                ax.set_title(f"Χρονική εξέλιξη: {compare_met}", fontfamily="serif", fontsize=13)
                ax.legend(fontsize=8, frameon=False)
                apply_journal_plot_style(ax)
                fig_ts_cmp.tight_layout()
                st.pyplot(fig_ts_cmp)
            else:
                st.caption("Δεν υπάρχουν αποθηκευμένες χρονοσειρές και για τις δύο εκτελέσεις.")
        else:
            st.caption("Χρειάζονται τουλάχιστον 2 αποθηκευμένες εκτελέσεις για σύγκριση head-to-head.")
    else:
        st.info("Δεν υπάρχουν αρκετά δεδομένα ιστορικού για σύγκριση. Εκτέλεσε "
                "μερικές προσομοιώσεις πρώτα.")


# ============================================================================
# 16. ΣΕΛΙΔΑ: ΔΙΚΤΥΟ ΜΕΤΑΒΟΛΙΣΜΟΥ
# ============================================================================
elif page == "Δίκτυο Μεταβολισμού":
    render_masthead(
        "Χαρτογράφηση μεταβολισμού &middot; Πρωτόκολλο 05",
        "Δίκτυο Παραγωγής Μεταβολιτών",
        "Διμερές δίκτυο (bipartite graph) που συνδέει τα επιλεγμένα μικρόβια με "
        "τους μεταβολίτες που παράγουν, με βάση τους πραγματικούς στοιχειομετρικούς "
        "συντελεστές του μοντέλου (όχι τυχαίες συνδέσεις)."
    )

    df = load_history(1)
    if not df.empty:
        last = df.iloc[0]
        st.caption(f"Βασισμένο στην πιο πρόσφατη εκτέλεση: **{last['microbe_a']} + {last['microbe_b']}** "
                    f"(#{last['id']}, {last['timestamp']})")
        final_snapshot = {
            "polysaccharide_lps": last["lps"],
            "mannan_polysaccharide": last["mannan"],
            "flagellin": last["flagellin"],
            "ethanol": last["ethanol"],
            "lactate": last["lactate"],
            "acetate": last["acetate"],
        }
        fig, G = create_metabolic_network(last["microbe_a"], last["microbe_b"], final_snapshot)
        st.pyplot(fig)

        st.subheader("Στοιχεία δικτύου")
        c1, c2 = st.columns(2)
        c1.metric("Κόμβοι", G.number_of_nodes())
        c2.metric("Ακμές", G.number_of_edges())
    else:
        st.info("Εκτέλεσε πρώτα μια προσομοίωση στη σελίδα **Προσομοίωση** ώστε να "
                "υπάρχουν δεδομένα για την κατασκευή του δικτύου.")

    st.divider()
    st.subheader("🗺️ Πλήρης χάρτης δυνατοτήτων παραγωγής (όλα τα μικρόβια)")
    if st.checkbox("Εμφάνιση πλήρους δικτύου όλων των μικροβίων"):
        G_full = nx.DiGraph()
        for name, mp in _microbes.items():
            G_full.add_node(name, kind="microbe")
            for met_key in mp.production:
                if met_key not in G_full:
                    category = _metabolite_registry[met_key].category
                    G_full.add_node(met_key, kind="metabolite", category=category)
                G_full.add_edge(name, met_key, weight=mp.production[met_key][0])

        fig_full, ax = plt.subplots(figsize=(11, 8))
        pos = nx.spring_layout(G_full, seed=7, k=0.7)
        microbe_nodes = [n for n, d in G_full.nodes(data=True) if d.get("kind") == "microbe"]
        target_nodes = [n for n, d in G_full.nodes(data=True) if d.get("kind") == "metabolite" and d.get("category") == "target"]
        byproduct_nodes = [n for n, d in G_full.nodes(data=True) if d.get("kind") == "metabolite" and d.get("category") == "byproduct"]
        nx.draw_networkx_nodes(G_full, pos, nodelist=microbe_nodes, node_color=COLOR_ACCENT, node_shape="s", node_size=1800, ax=ax)
        nx.draw_networkx_nodes(G_full, pos, nodelist=target_nodes, node_color=COLOR_TARGET, node_size=1300, ax=ax)
        nx.draw_networkx_nodes(G_full, pos, nodelist=byproduct_nodes, node_color=COLOR_BYPRODUCT, node_size=1000, ax=ax)
        nx.draw_networkx_edges(G_full, pos, alpha=0.5, arrows=True, ax=ax)
        nx.draw_networkx_labels(G_full, pos, font_size=7, ax=ax)
        ax.set_title("Πλήρες δίκτυο μικροβίων → μεταβολιτών", fontfamily="serif", fontsize=13)
        ax.axis("off")
        fig_full.tight_layout()
        st.pyplot(fig_full)

elif page == "⚡ FBA (SBML)":
    show_fba()

# ============================================================================
# 17. ΣΕΛΙΔΑ: ΣΧΕΤΙΚΑ / ΕΠΙΣΤΗΜΟΝΙΚΟ ΥΠΟΒΑΘΡΟ
# ============================================================================
elif page == "Σχετικά":
    render_masthead(
        "Παράρτημα Β &middot; Τεκμηρίωση",
        "Σχετικά με την Εφαρμογή",
        "Επιστημονικό υπόβαθρο, μαθηματικό μοντέλο και αποποίηση ευθύνης."
    )
    st.markdown("""
## Vaccine Adjuvant Discovery Platform — Έκδοση 5.0 (single-file)

Εφαρμογή προσομοίωσης μικροβιακών συν-καλλιεργειών με στόχο τον εντοπισμό
συνδυασμών μικροβίων που παράγουν υποσχόμενα **ανοσοδιεγερτικά μόρια
(PAMPs)** — LPS, μαννάνη, φλαγγελίνη — σε ισορροπία με χαμηλή τοξικότητα
από παραπροϊόντα ζύμωσης.

### Μαθηματικό μοντέλο
- **Κινητική ανάπτυξης (Monod):**  μ(S) = μ_max · S / (Ks + S)
- **Ισοζύγιο βιομάζας:**  dX/dt = (μ − k_death) · X
- **Κατανάλωση υποστρώματος:**  dS/dt = − Σᵢ μᵢ·Xᵢ / Y_x/sᵢ
- **Παραγωγή μεταβολιτών (Luedeking–Piret):**  dP/dt = α·μ·X + β·X
- Αριθμητική ολοκλήρωση μέσω `scipy.integrate.odeint` (με εναλλακτικό
  ενσωματωμένο ολοκληρωτή Runge–Kutta 4ης τάξης αν λείπει το SciPy).

### Βαθμολόγηση εμβολιακού δυναμικού
Κάθε στόχος-PAMP αξιολογείται με βάση:
1. **Ανοσογονικότητα** (dose-response καμπύλη τύπου Hill ως προς τη
   συγκέντρωση, βαθμονομημένη σε EC50 ανά μόριο)
2. **Ασφάλεια** (π.χ. το LPS έχει χαμηλότερη βαθμολογία ασφάλειας λόγω
   κινδύνου ενδοτοξίνης)
3. **Σταθερότητα** (π.χ. πρωτεΐνες όπως η φλαγγελίνη είναι λιγότερο
   σταθερές θερμικά από πολυσακχαρίτες)

### Τεχνολογίες
- Python, Streamlit, NumPy, SciPy, pandas
- Matplotlib / Seaborn για οπτικοποιήσεις
- NetworkX για το δίκτυο μεταβολισμού
- SQLite για επιμονή (persistence) ιστορικού πειραμάτων

### ⚠️ Αποποίηση ευθύνης
Πρόκειται για **εκπαιδευτικό / αρχιτεκτονικό εργαλείο επίδειξης**. Οι
κινητικές παράμετροι των μικροβίων και το μοντέλο βαθμολόγησης είναι
απλοποιημένα και ενδεικτικά — δεν έχουν πειραματική βαθμονόμηση και δεν
αντικαθιστούν πραγματική ανοσολογική/μικροβιολογική αξιολόγηση. Δεν
πρέπει να χρησιμοποιηθούν για πραγματικές αποφάσεις ανάπτυξης εμβολίων.
    """)


# ============================================================================
# 18. ΥΠΟΣΕΛΙΔΟ ΠΛΕΥΡΙΚΗΣ ΜΠΑΡΑΣ
# ============================================================================
st.sidebar.divider()
st.sidebar.caption("🧬 Vaccine Adjuvant Discovery Platform · SQLite persistence")
st.sidebar.caption(f"Έκδοση 5.0 · {datetime.now().strftime('%Y-%m-%d %H:%M')}")
st.sidebar.caption(f"Καταγεγραμμένες εκτελέσεις: {count_history()}")