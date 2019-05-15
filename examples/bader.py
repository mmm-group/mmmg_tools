from mmmg_tools.objects import Charge
from pickle import dump
from sys import argv

if __name__ == '__main__':
    argv.append('')
    chg = Charge.from_file(argv[1])
    chg.bader_calc(args=' '.join(argv[2:]))
    chg.bader.write(prefix=argv[1]+'-')
    print(f'\n  WRITING PICKLE: ./bdr.p')
    with open('bdr.p','+wb') as f:
        dump(chg,f)

