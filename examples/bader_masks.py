from mmmg_tools.objects import Charge

atom_list = [0, 8]
chg = Charge.from_file('CHGCAR')
chg.bader_calc(spin_flag=True, export_mode=('atoms', atom_list))
