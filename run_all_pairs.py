"""
run_all_pairs.py – Τρέχει όλους τους συνδυασμούς των μικροβίων
και βγάζει τους TOP 10 με βάση βαθμολογία.
"""

import pandas as pd
from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator
from micromol_scorer import create_default_scorer
from itertools import combinations
import time

# Φόρτωση μικροβίων και μεταβολιτών
microbes = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes, metabolites)
scorer = create_default_scorer()

# Λίστα από όλα τα ονόματα μικροβίων
microbe_names = list(microbes.keys())
print(f"🧫 Βρέθηκαν {len(microbe_names)} μικρόβια. Δημιουργία όλων των ζευγαριών...")

# Όλοι οι συνδυασμοί (2)
pairs = list(combinations(microbe_names, 2))
print(f"🔬 Θα τρέξουν {len(pairs)} προσομοιώσεις... (θα πάρει λίγη ώρα)")

results = []
start_total = time.time()

for i, (a, b) in enumerate(pairs, 1):
    print(f"  [{i}/{len(pairs)}] {a} + {b} ...", end=" ", flush=True)
    try:
        # Προσομοίωση με 100 mM γλυκόζη, 48 ώρες
        initial = {"glucose": 100.0}
        res = sim.simulate_coculture([a, b], initial, duration=48.0)
        final = res.get("final_concentrations", {})

        # Βαθμολόγηση
        score = 0.0
        for met in ["polysaccharide_lps", "mannan_polysaccharide", "flagellin"]:
            if met in final:
                sc = scorer.score_molecule(met, met, "", None)
                score += sc.overall_vaccine_score

        # Τοξικότητα
        toxicity = sum(final.get(t, 0) for t in ["ethanol", "lactate", "acetate"])
        combined_score = score - toxicity * 0.5  # απλό μέτρο

        results.append({
            "microbe_a": a,
            "microbe_b": b,
            "LPS": final.get("polysaccharide_lps", 0),
            "mannan": final.get("mannan_polysaccharide", 0),
            "flagellin": final.get("flagellin", 0),
            "score": combined_score,
            "toxicity": toxicity
        })
        print("✅")
    except Exception as e:
        print(f"❌ Λάθος: {e}")

elapsed = time.time() - start_total
print(f"\n✅ Ολοκληρώθηκαν {len(results)} προσομοιώσεις σε {elapsed:.2f} δευτερόλεπτα.")

# Ταξινόμηση και αποθήκευση
df = pd.DataFrame(results)
df_sorted = df.sort_values("score", ascending=False)

# Τοπ 10
top10 = df_sorted.head(10)
print("\n🏆 TOP 10 ΣΥΝΔΥΑΣΜΟΙ:")
print(top10[["microbe_a", "microbe_b", "score", "LPS", "mannan", "flagellin", "toxicity"]].to_string(index=False))

# Αποθήκευση σε CSV
top10.to_csv("top10_pairs.csv", index=False)
df_sorted.to_csv("all_pairs_results.csv", index=False)
print("\n📁 Αποτελέσματα αποθηκεύτηκαν σε:")
print("  - top10_pairs.csv")
print("  - all_pairs_results.csv")
