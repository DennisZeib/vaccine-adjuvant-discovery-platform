"""
download_bigg_models.py

Download genome-scale metabolic models from BiGG database.
"""

from cobra.io import load_model

# List of BiGG model IDs to download
model_ids = [
    "iML1515",   # E. coli
    "iND750",    # Yeast
    "iB21_1397", # E. coli BL21
    "iECD_1391", # E. coli BL21(DE3)
    "iBWG_1329", # E. coli BL21(DE3)
    "iSB619",    # Bacillus subtilis
    "iYO844",    # Lactobacillus plantarum
]

for model_id in model_ids:
    try:
        print(f"Downloading {model_id}...")
        model = load_model(model_id)
        model.save(f"{model_id}.xml", overwrite=True)
        print(f"  Saved: {model_id}.xml")
    except Exception as e:
        print(f"  Failed: {model_id} - {e}")