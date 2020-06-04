from mmmg_tools.objects import Charge

chg = Charge.from_file('CHGCAR')
chg.bader_calc(threads=16)
