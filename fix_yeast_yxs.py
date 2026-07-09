"""
Recalibrate Saccharomyces cerevisiae Y_xs using real literature data
========================================================================
Real data found (July 2026):
- Wild-type CEN.PK 113-5D biomass yield: 0.11 g/g glucose (fermentative,
  Crabtree-positive conditions) - DOI: 10.1016/j.synbio.2022.06.004
- Engineered Crabtree-relieved strain: 0.21 g/g glucose (same source)

Current simulator value: Y_xs = 0.40 (roughly 2-4x higher than real
wild-type yeast under normal Crabtree-positive conditions).

This script updates Y_xs to 0.20 (upper end of the real range, closer
to the engineered/best-case strain) as a defensible middle estimate,
and adds a citation documenting the real range and the approximation
made. This is NOT a precise fit - it's a documented, literature-
grounded correction replacing an unvalidated placeholder.
"""

with open("app.py", encoding="utf-8") as f:
    content = f.read()

old_block = '''            mu_max=0.45,
            Ks=10.0,
            Y_xs=0.40,
            k_death=0.009,
            production={
                "mannan_polysaccharide": (1.30, 0.08),
                "ethanol": (1.50, 0.10),
            },
            description="Ευκαρυωτική ζύμη. Το κυτταρικό της τοίχωμα είναι "
                         "πλούσιο σε μαννάνη (Dectin-1/TLR2)· παράγει αιθανόλη "
                         "μέσω ζύμωσης (υψηλή τοξικότητα σε υψηλή γλυκόζη).",'''

new_block = '''            mu_max=0.45,
            Ks=10.0,
            Y_xs=0.20,  # RECALIBRATED (was 0.40) - see citation below
            k_death=0.009,
            production={
                "mannan_polysaccharide": (1.30, 0.08),
                "ethanol": (1.50, 0.10),
            },
            description="Ευκαρυωτική ζύμη. Το κυτταρικό της τοίχωμα είναι "
                         "πλούσιο σε μαννάνη (Dectin-1/TLR2)· παράγει αιθανόλη "
                         "μέσω ζύμωσης (υψηλή τοξικότητα σε υψηλή γλυκόζη). "
                         "[Y_xs recalibrated July 2026: real wild-type "
                         "S. cerevisiae biomass yield is ~0.11 g/g glucose "
                         "under Crabtree-positive conditions, up to 0.21 g/g "
                         "in engineered Crabtree-relieved strains (DOI "
                         "10.1016/j.synbio.2022.06.004). Original toy value "
                         "of 0.40 was ~2-4x too high. Set to 0.20 as a "
                         "documented middle estimate, not a precise fit.]",'''

if old_block not in content:
    print("ERROR: could not find exact yeast block to replace. No changes made.")
    print("This likely means whitespace differs from what was extracted earlier.")
else:
    content = content.replace(old_block, new_block, 1)
    with open("app.py", "w", encoding="utf-8") as f:
        f.write(content)
    print("Done. Yeast Y_xs recalibrated from 0.40 to 0.20 with citation.")