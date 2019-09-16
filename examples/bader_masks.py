from mmmg_tools.objects import Charge
from pickle import dump
import gzip

atom_list = [0, 8]
print('Reading charge file.')
chg = Charge.from_file('notebooks/CHGCAR')
print('Charge file read.')
chg.bader_calc(args='-p')
chg.bader.write(prefix='CHGCAR-')

for site in atom_list:
    mask = chg.bader.masks.ion_mask(site, chg.bader.bader_ion)
    mask.apply(chg).write_file(f"{site}-CHGCAR")

print(f"\nwriting gzip'd pickle: ./bdr.pgz")
with gzip.open('bdr.pgz','wb') as f:
    dump(chg,f)


