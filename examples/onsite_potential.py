from mmmg_tools.objects import Potential

pot = Potential.from_file('LOCPOT')
r = {'Mg': 1.4, 'O': 1.6}
Z = 1

print(pot.onsite_electrostatic(Z=Z, r=r, sites=None, upsample=2))
