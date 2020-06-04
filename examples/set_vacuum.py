from mmmg_tools.objects import Structure

structure = Structure.from_file('CONTCAR')
for i in range(3):
    print(f'{i} No PB')
    s = Structure.from_structure(structure.structure)
    s.set_vacuum(15, axis=i)
    s.write_file(f'{i}-nopb-POSCAR')
    print(f'{i} PB')
    s = Structure.from_structure(structure.structure)
    s.set_vacuum(15, axis=i, add_pb_atoms=True)
    s.write_file(f'{i}-pb-POSCAR')
