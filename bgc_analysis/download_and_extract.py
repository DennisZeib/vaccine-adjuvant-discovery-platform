import os
import urllib.request
import gzip
import shutil

url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/686/825/GCA_000686825.1_ASM68682v1/GCA_000686825.1_ASM68682v1_genomic.gbff.gz'
base = os.path.join(os.getcwd(), 'bgc_analysis')
genomes = os.path.join(base, 'genomes')
os.makedirs(genomes, exist_ok=True)

gz_path = os.path.join(genomes, 'GCA_000686825.1_ASM68682v1_genomic.gbff.gz')
out_path = os.path.join(genomes, 'GCA_000686825.1_ASM68682v1_genomic.gbff')

print('Downloading', url)
urllib.request.urlretrieve(url, gz_path)
print('Downloaded to', gz_path)

print('Decompressing...')
with gzip.open(gz_path, 'rb') as f_in:
    with open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
print('Decompressed to', out_path)

print('Files in genomes:')
for name in sorted(os.listdir(genomes)):
    p = os.path.join(genomes, name)
    print(name, os.path.getsize(p))
