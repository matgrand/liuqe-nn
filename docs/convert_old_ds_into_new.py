# convert old ds into new ds format
import os
from scipy.io import loadmat, savemat
from tqdm import tqdm

old_full = 'og_raw_ds.tar.gz'
old_smol = 'og_raw_ds_smol.tar.gz'
new_full = 'raw_ds.tar.gz'
new_smol = 'raw_ds_smol.tar.gz'

# datasets = [old_smol]
datasets = [old_smol, old_full]
to_change = [
    ('Bm0','Bm'),
    ('Ff0','Ff'),
    ('Ia0','Ia'),
    ('Ip0','Ip'),
]

for ds in datasets:
    print(f"Processing {ds} ...")
    os.system(f"rm -rf raw_ds")
    print(f"Extracting {ds} ...")
    os.system(f"tar -xf {ds}")
    print(f'Extracted {ds}. Now processing .mat files ...')
    files = os.listdir('raw_ds')
    for f in tqdm(files):
        assert f.endswith('.mat')
        data = loadmat(f'raw_ds/{f}')
        # change keys, but keep also the old ones
        for old, new in to_change:
            data[new] = data[old]
        savemat(f'raw_ds/{f}', data)
    print(f'Processed {len(files)} .mat files. Now compressing ...')
    os.system(f"tar -czf {ds.replace('og_','ogc_')} raw_ds")

