from mmmg_tools.objects import Structure

structure = Structure.from_file('CONTCAR')
for i in range(8,16):
    structure.set_vacuum(i)
    structure.write_file(f'{i}-vac-POSCAR')
