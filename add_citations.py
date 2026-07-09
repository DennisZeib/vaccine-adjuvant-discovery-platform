import re

with open("app.py", encoding="utf-8") as f:
    content = f.read()

replacements = [
    (
        '"note": "Ισχυρό PAMP· υψηλή ανοσογονικότητα αντισταθμίζεται από "\n                "κίνδυνο τοξικότητας ενδοτοξίνης σε υψηλές δόσεις.",',
        '"note": "Ισχυρό PAMP· υψηλή ανοσογονικότητα αντισταθμίζεται από "\n                "κίνδυνο τοξικότητας ενδοτοξίνης σε υψηλές δόσεις.",\n        "citation": "ec50 is simulator-internal (toy units), not a real-world value. Real TLR4/LPS potency data was not found in literature search as of 2026-07.",'
    ),
    (
        '"note": "Καλά τεκμηριωμένο ανοσοδιεγερτικό μόριο (χρησιμοποιείται σε "\n                "ερευνητικές πλατφόρμες εμβολίων). Μειωμένη σταθερότητα λόγω "\n                "πρωτεϊνικής φύσης.",',
        '"note": "Καλά τεκμηριωμένο ανοσοδιεγερτικό μόριο (χρησιμοποιείται σε "\n                "ερευνητικές πλατφόρμες εμβολίων). Μειωμένη σταθερότητα λόγω "\n                "πρωτεϊνικής φύσης.",\n        "citation": "Real TLR5-binding EC50 = 2.4 pM (Smith et al. 2013, DOI 10.1002/bit.24903). Human vaccine dose-response threshold: 6-10mcg effective (DOI 10.1016/j.vaccine.2017.09.070). ec50 field above is simulator-internal toy units, NOT this real value \\u2014 unit systems are incompatible without full model recalibration.",'
    ),
    (
        '"note": "Πολυσακχαρίτης κυτταρικού τοιχώματος ζυμών. Καλά "\n                         "ανεκτός, χρησιμοποιείται σε προσεγγίσεις adjuvant.",',
        '"note": "Πολυσακχαρίτης κυτταρικού τοιχώματος ζυμών. Καλά "\n                         "ανεκτός, χρησιμοποιείται σε προσεγγίσεις adjuvant.",\n        "citation": "Mouse dose-response: 1-7mcg -> cellular immunity, >7mcg -> humoral dominates (DOI 10.1007/s002620050449). ec50 above is simulator-internal toy units.",'
    ),
]

count = 0
for old, new in replacements:
    if old in content:
        content = content.replace(old, new)
        count += 1
    else:
        print(f"WARNING: pattern not found (skipped): {old[:60]}...")

with open("app.py", "w", encoding="utf-8") as f:
    f.write(content)

print(f"Done. {count}/3 citation blocks added to app.py")