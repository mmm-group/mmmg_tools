from mmmg_tools.objects import Charge
from pickle import dump
import argparse
import gzip

if __name__ == '__main__':
    psr = argparse.ArgumentParser(
        description='Python wrapped Bader charge analysis.')
    psr.add_argument('filename', metavar='CHGCAR', type=str,
        help='Path to charge density file.')
    psr.add_argument('-ref', required=False, type=str,
        help='Path to reference charge density file.')
    psr.add_argument('-b', required=False, type=str.lower, 
        choices=['neargrid','ongrid','weight'],
        help=('Use the default near-grid bader partitioning, ' 
            'the original on-grid based algorithm, ' 
            'or the weight method of Yu and Trinkle.'))
    psr.add_argument('-r', required=False, type=int,
        help=('Number of refine edge iterations. To access the ' 
            'old iteration method, checking all edge points, '
            '-2 must be passed.'))
    psr.add_argument('-vac', nargs='?', required=False, const=1e-3, type=float,
        help=('Assign low density points to vacuum.'))
    psr.add_argument('-m', required=False, type=str.lower,
        choices=['known','max'], help='Determines how trajectories terminate.')
    psr.add_argument('-s', required=False, type=float,
        help='Stepsize.')
    psr.add_argument('-t', required=False, type=float,
        help='Significant Bader volume tolerance.')
    psr.add_argument('-p', required=False, const='', action='store_const',
        help='Output masks of Bader volumes.')
    args = psr.parse_args()

    calc_args = ''
    for k, v in args.__dict__.items():
        if k != 'filename' and k != 'ref' and v is not None:
            calc_args += f'-{k} {v} '
    print(f'Command line args:\n  {calc_args}')

    print('Reading charge file.')
    chg = Charge.from_file(args.filename)
    print('Charge file read.')
    ref = None
    if args.ref is not None:
        print('Reading reference charge file.')
        refchg = Charge.from_file(args.ref)
        print('Reference charge file read.')
        ref = refchg.data['total']
    chg.bader_calc(ref=ref, args=calc_args)
    chg.bader.write(prefix=args.filename+'-')
    print(f"\nwriting gzip'd pickle: ./bdr.pgz")
    with gzip.open('bdr.pgz','wb') as f:
        dump(chg,f)


