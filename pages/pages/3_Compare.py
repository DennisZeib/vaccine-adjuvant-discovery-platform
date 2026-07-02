import streamlit as st
import pandas as pd
from database import load_history

st.title("📊 Σύγκριση Αποτελεσμάτων")

df = load_history()
if not df.empty:
    st.dataframe(df[["microbe_a", "microbe_b", "glucose", "score", "toxicity"]], use_container_width=True)

    # Γράφημα σύγκρισης
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(df["score"], df["toxicity"], alpha=0.7, s=100, color="#3498db")
    for i, row in df.iterrows():
        ax.annotate(f"{row['microbe_a']}+{row['microbe_b']}", (row["score"], row["toxicity"]), fontsize=8)
    ax.set_xlabel("Αποτελεσματικότητα")
    ax.set_ylabel("Τοξικότητα")
    ax.set_title("Σύγκριση εκτελέσεων")
    ax.grid(alpha=0.3)
    st.pyplot(fig)
else:
    st.info("Δεν υπάρχουν δεδομένα.")