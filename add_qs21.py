"""
Add QS-21 to _KNOWLEDGE_BASE
==============================
QS-21 is a standalone candidate (not produced by any microbe in the
co-culture simulator), so it's scored independently via
scorer.score_molecule() at a chosen dose, not through run_monte_carlo().
"""

with open("app.py", encoding="utf-8") as f:
    content = f.read()

marker = '"mannan_polysaccharide": {'
if marker not in content:
    print("ERROR: could not find insertion point (mannan_polysaccharide entry). Aborting - no changes made.")
else:
    qs21_entry = '''"saponin_qs21": {
        "immunogenicity": 88.0,
        "safety": 45.0,
        "stability": 40.0,
        "receptor": "NLRP3 inflammasome / lysosomal destabilization",
        "ec50": 10.0,
        "hill_n": 1.5,
        "note": "Standalone saponin adjuvant from Quillaja saponaria bark. "
                "NOT produced by any microbe in this simulator - score "
                "independently via scorer.score_molecule() at a chosen dose, "
                "not through run_monte_carlo(). Component of licensed "
                "adjuvant systems AS01 (Shingrix) and Mosquirix.",
        "citation": "Mouse dose-limiting toxicity: 20mcg/dose caused ~10% body weight loss (2020 conference abstract, no DOI). Human clinical dose: 25-50mcg per dose, well tolerated in Phase 1 (DOI 10.1016/s2666-5247(23)00410-x). Known liabilities: hemolytic toxicity, hydrolytic instability (DOI 10.1016/j.ejmech.2025.118223). ec50 above is simulator-internal placeholder - no real EC50/potency-curve data was found for QS-21 in toy-unit-compatible form.",
    },
    ''' + marker

    content = content.replace(marker, qs21_entry, 1)

    with open("app.py", "w", encoding="utf-8") as f:
        f.write(content)

    print("Done. QS-21 added to _KNOWLEDGE_BASE as a standalone candidate.")