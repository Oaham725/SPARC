# @Author:MH, LXY
# -*- coding = 'utf-8' -*-
# @Time:2020/11/20 08:39
# @File: main.py
# @Software: PyCharm


'''
In this programm, pdb files are read in to generate an instance of the protein data structure.
The direction of a residue can be further calculated.
'''

import re
import glob
import warnings
from utils.positioner import Protein

warnings.filterwarnings('ignore')

if __name__ == '__main__':
    pdbFiles = glob.glob('inpFiles/*.pdb')
    xlsxDir = r'utils/XlsxDir'
    txtDir = r'outFiles'
    SavePttern = re.compile(r'(.*?).pdb')

    for file in pdbFiles:
        filename = file.split('.pdb')[0]
        saveNames = f"{filename.split('/')[-1]}.out"
        myStructure = Protein(file)
        myStructure.toStructure()
        myStructure.getResidue(xlsxDir, txtDir, saveNames, mode='sparc', z_init=0.0, z_end=0.0)

    print('Finished!')
