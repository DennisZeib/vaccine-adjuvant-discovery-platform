"""
Recalibrate Salmonella enterica mu_max using real literature data
=====================================================================
Real data found (July 2026):
- S. typhimurium in glucose-minimal medium: doubling time ~48 min
  -> mu_max = ln(2)/0.8h = ~0.87/h
  (DOI: 10.1128/jb.114.3.966-973.1973)
- S. typhimurium on beef/chicken, 30C: doubling time 0.74h
  -> mu_max = ln(2)/0.74h = ~0.94/h
  (DOI: 10.1128/aem.58.11.3482-3487.1992, different conditions but
  consistent order of magnitude)

Current simulator value: mu_max = 0.70/h - about 20-25% lower than
the glucose-minimal-medium real value, which is the most directly
comparable condition to the simulator's own glucose-based model.

This script updates mu_max to 0.85 (close to the glucose-minimal-medium
real value, slightly conservative given only one directly-comparable
data point) and adds a citation.
"""

with open("app.py", encoding="utf-8") as f:
    content = f.read()

old_block = '''            mu_max=0.70,
            Ks=8.0,
            Y_xs=0.45,
            k_death=0.018,
            production={
                "polysaccharide_lps": (1.10, 0.08),
                "flagellin": (0.55, 0.04),
                "acetate": (0.45, 0.02),
            },
            description="Gram-αρνητικό, πολύ κινητό. Υψηλή παραγωγή LPS & "
                         "φλαγγελίνης· κλασικό μοντέλο μελέτης PAMPs.",'''

new_block = '''            mu_max=0.85,  # RECALIBRATED (was 0.70) - see citation below
            Ks=8.0,
            Y_xs=0.45,
            k_death=0.018,
            production={
                "polysaccharide_lps": (1.10, 0.08),
                "flagellin": (0.55, 0.04),
                "acetate": (0.45, 0.02),
            },
            description="Gram-αρνητικό, πολύ κινητό. Υψηλή παραγωγή LPS & "
                         "φλαγγελίνης· κλασικό μοντέλο μελέτης PAMPs. "
                         "[mu_max recalibrated July 2026: real S. typhimurium "
                         "doubling time in glucose-minimal medium is ~48min "
                         "(mu_max ~0.87/h), DOI 10.1128/jb.114.3.966-973.1973. "
                         "Original toy value of 0.70 was ~20-25% low. Set to "
                         "0.85 as a close, slightly conservative estimate.]",'''

if old_block not in content:
    print("ERROR: could not find exact Salmonella block to replace. No changes made.")
else:
    content = content.replace(old_block, new_block, 1)
    with open("app.py", "w", encoding="utf-8") as f:
        f.write(content)
    print("Done. Salmonella mu_max recalibrated from 0.70 to 0.85 with citation.")