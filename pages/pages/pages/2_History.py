import streamlit as st
import pandas as pd
from database import load_history, delete_all

st.title("📜 Ιστορικό Εκτελέσεων")

df = load_history()
if not df.empty:
    st.dataframe(df, use_container_width=True)

    if st.button("🗑️ Διαγραφή όλων"):
        delete_all()
        st.rerun()

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("📥 Λήψη CSV", data=csv, file_name="history.csv", mime="text/csv")
else:
    st.info("Δεν υπάρχουν αποθηκευμένες εκτελέσεις.")