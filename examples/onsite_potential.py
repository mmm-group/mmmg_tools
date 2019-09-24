from mmmg_tools.objects import Potential

pot = Potential.from_file('LOCPOT')
Z = {'Mg':4.825, 'O':4.775}
r = 1.6

print(pot.onsite_electrostatic(Z=Z,r=r,sites=None,upsample=2))
