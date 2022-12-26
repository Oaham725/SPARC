"""
This programm reads in .pdb files, create a structure object to extract information about the position of certain atoms of proteins.
Then calculate the SERS spectra of protein of interest according to the information extracted as well as invoked files.
More about the structure class, please refer to:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc159
"""

import numpy as np
from Bio.PDB.PDBParser import PDBParser
import utils.rotateMatrix as RM
import utils.matrix as Matrix
from .gpostGetZ import get_z


def Centre_Two_vector(x, y):
    m = x + y


def SaveFiles(DataSet, SavePath):
    np.savetxt(SavePath, DataSet, fmt='%.20f')
    with open(SavePath, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("{} 2\n".format(DataSet.shape[0]) + content)
        f.close()


def MatrixPerforming(MatrixDir, res_name, Z, T_set, avg_atoms, mode, z_init=0.0, z_end=0.0):

    MatrixSet, Label, WaveNumber, Width, v_mode = Matrix.LoadMartix(FilePath=f'{MatrixDir}/{res_name}.xlsx')
    
    S_Set, WaveSet = Matrix.ProcessMartix(
        MatrixSet=MatrixSet,
        WaveNumber=WaveNumber,
        Label=Label,
        v_mode=v_mode,
        Z=Z,
        z_init=z_init,
        z_end=z_end,
        T_Set=T_set,
        avg_atoms=avg_atoms,
        res_name=res_name,
        mode=mode,
    )
    temp = Matrix.SaveData(S_Set, WaveSet, Width)
    return temp


class Protein:

    def __init__(self, filename, structure_id="myProtein"):
        self.filename = filename
        self.structure_id = structure_id

    def toStructure(self):
        parser = PDBParser(PERMISSIVE=1)
        self.structure = parser.get_structure(self.structure_id, self.filename)

    def getResidue(self, MatrixDir, TxtSaveDir, SaveName, mode='sparc', z_init=0.0, z_end=0.0):
        # f = open("POI.txt", "+a")

        # structure contains model?
        if len(self.structure) != 1:
            return
        else:
            self.model = self.structure[0]
        # model contains chain?
        if len(self.model) != 1:
            return
        else:
            self.chain = self.model["X"]
        # print residues (sequence)
        self.aaList = self.chain.get_list()
        print('组成该蛋白的氨基酸个数：', len(self.aaList), '\nfile = %s' % self.filename)

        T_zeros = np.zeros((3, 3))
        DataSet = np.zeros((1, 3))
        for aa in self.aaList:
            atoms = list(aa.get_atoms())  # import all message of atoms in this protein(peptide)
            Data_TRP = np.zeros((1, 3))
            Data_HIS = np.zeros((1, 3))
            Data_ASP = np.zeros((1, 3))
            Data_ALA = np.zeros((1, 3))
            Data_GLU = np.zeros((1, 3))
            Data_PHE = np.zeros((1, 3))
            Data_ARG = np.zeros((1, 3))
            Data_SER = np.zeros((1, 3))
            Data_GLY = np.zeros((1, 3))
            Data_TYR = np.zeros((1, 3))
            Data_VAL = np.zeros((1, 3))
            Data_THR = np.zeros((1, 3))
            Data_ILE = np.zeros((1, 3))
            Data_PRO = np.zeros((1, 3))
            Data_MET = np.zeros((1, 3))
            Data_ASN = np.zeros((1, 3))
            Data_CYS = np.zeros((1, 3))
            Data_GLN = np.zeros((1, 3))
            Data_LEU = np.zeros((1, 3))
            Data_LYS = np.zeros((1, 3))
            Data_HEM = np.zeros((1, 3))
            Data_DA = np.zeros((1, 3))
            Data_DG = np.zeros((1, 3))
            Data_DC = np.zeros((1, 3))
            Data_DT = np.zeros((1, 3))
            if aa.get_resname(
            ) == "HIS":  # if residue name == HIS, then catch several atoms as follows:
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[6].get_vector()  # atom：CD2
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CB
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[8].get_vector()  # atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                '''Standard orientation'''
                S_Atom1 = np.array([-1.291888,   -0.676124,    0.269808])  # atom：CA
                S_Atom2 = np.array([2.384911,   -1.094687,   -0.002806])  # atom：CD2
                S_Atom3 = np.array([1.226849,   -0.410479,   -0.289245])  # atom: CG
                S_Atom4 = np.array([-0.117089,   -0.898493,   -0.725900])  # atom: CB
                S_Atom5 = np.array([-2.262978,   -1.758756,    0.136786])  # atom: N
                S_Atom6 = np.array([-2.042862,    0.634826,    0.082444])  # atom: C

                # vector_1 = (Atom6 - Atom1)[0:3]  # import vector for c-c(o)
                vector_2 = (Atom5 - Atom1)[0:3]  # import vector for c-n, [-0.97, -1.083, -0,131]
                vector_3 = (Atom4 - Atom1)[0:3]  # import vector for ch2, [1.175, -0.224, -0.996]
                vector_4 = (Atom4 - Atom3)[0:3]  # import vector 1 for imidazole ring, [-1.344, -0.488, -0.436]
                vector_5 = (Atom3 - Atom2)[0:3]  # import vector 2 for imidazole ring, [-1.159, 0.684, -0.284]
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom4 - S_Atom1)
                S_vector_4 = (S_Atom4 - S_Atom3)
                S_vector_5 = (S_Atom3 - S_Atom2)

                T1 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)
                T2 = RM.Calculate(S_vector_4, vector_4, S_vector_3, vector_3)
                T3 = RM.Calculate(S_vector_2, vector_2, S_vector_3, vector_3)

                T_Set = [T1, T2, T3]
                #计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                #计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2] - 1
                # # For amino acid which do not need 5 transforming matrix
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''Matrix performing'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_HIS = np.append(Data_HIS, temp, axis=0).reshape(-1, 3)

            if aa.get_resname() == "TYR":  # 如果残基为酪氨酸，则采一部分原子
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[10].get_vector()  # tyr_atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([2.113680,    0.655217,   0.256669])  # atom：CA
                S_Atom2 = np.array([1.073072,    1.135336,   -0.775158])  # atom：CB
                S_Atom3 = np.array([-0.343785,    0.692537,   -0.485414])  # atom: CG
                S_Atom4 = np.array([-0.901573,   -0.423641,   -1.113691])  # atom: CD1
                S_Atom5 = np.array([3.420072,    1.234632,   -0.066485])  # atom: N
                S_Atom6 = np.array([2.169120,   -0.873387,    0.283609])  # tyr_atom: C
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.066, -1.528, 0.013]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n, [1.306, 0.59, -0.309]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future, [-1.037, 0.484, -1.033]
                vector_4 = (Atom3 - Atom2)[0:3]  # 苯环向量1, [-1.418, -0.44, 0.288]
                vector_5 = (Atom3 - Atom4)[0:3]  # 苯环向量2, [0.556, 1.123, 0.619]
                S_vector_1 = (S_Atom6 - S_Atom1)  
                S_vector_2 = (S_Atom5 - S_Atom1)  
                S_vector_3 = (S_Atom2 - S_Atom1)  
                S_vector_4 = (S_Atom3 - S_Atom2)  
                S_vector_5 = (S_Atom3 - S_Atom4)  
                T1 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)
                T2 = RM.Calculate(S_vector_4, vector_4, S_vector_3, vector_3)
                T3 = RM.Calculate(S_vector_1, vector_1, S_vector_2, vector_2)
                T_Set = [T1, T2, T3]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2] - 1
                # # For amino acid which do not need 5 transforming matrix
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_TYR = np.append(Data_TYR, temp).reshape(-1, 3)

            if aa.get_resname() == "PHE":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[9].get_vector()  # atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([-1.661550,   -0.671746,    0.106417])  # atom：CA
                S_Atom2 = np.array([-0.512005,   -0.792883,   -0.930952])  # atom：CB
                S_Atom3 = np.array([0.873190,   -0.411015,   -0.456836])  # atom: CG
                S_Atom4 = np.array([1.556041,    0.675593,   -1.009781])  # atom: CD1
                S_Atom5 = np.array([-2.794987,   -1.482730,   -0.329924])  # atom: N
                S_Atom6 = np.array([-2.159217,    0.750972,    0.327407])  # atom: C
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.487, 1.424, 0.233]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n, [-1.137, -0.799, -0.45]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.151, -0.118, -1.036]
                vector_4 = (Atom3 - Atom2)[0:3]  # 苯环向量1, [1.384, 0.386, 0.474]
                vector_5 = (Atom3 - Atom4)[0:3]  # 苯环向量2, [0.686, 1.078, -0.567]
                # 苯环标准向量1 苯环标准向量2
                '''old version
                phStd1 = np.array([1.384, 0.386, 0.474])
                phStd2 = np.array([0.686, 1.078, -0.567])
                CH2td1 = np.array([1.151, -0.118, -1.036])
                # 苯环T记为1
                T1 = RM.Calculate(phStd1, vector_4, phStd2, vector_5)
                # print(T1 @ T1.T)
                T2 = RM.Calculate(phStd1, vector_4, CH2td1, vector_3)
                # c - c(o)
                CCO_Std = np.array([-0.487, 1.424, 0.233])
                # c-n
                CN_Std = np.array([-1.137, -0.799, -0.45])
                T3 = RM.Calculate(CCO_Std, vector_1, CN_Std, vector_2)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o), [-0.487, 1.424, 0.233]
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n, [-1.137, -0.799, -0.45]
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用, [1.151, -0.118, -1.036]
                S_vector_4 = (S_Atom3 - S_Atom2)  # 苯环向量1, [1.384, 0.386, 0.474]
                S_vector_5 = (S_Atom3 - S_Atom4)  # 苯环向量2, [0.686, 1.078, -0.567]
                T1 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)
                T2 = RM.Calculate(S_vector_4, vector_4, S_vector_3, vector_3)
                T3 = RM.Calculate(S_vector_1, vector_1, S_vector_2, vector_2)
                T_Set = [T1, T2, T3]
                # 计算均匀场模型所用坐标：
                # X = 1.5
                # Y = 1.5
                # Z = 1.5
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]   #在计算界面SERS时需通过减小距离，模拟其散射界面大于其他氨基酸
                # T最多是5个，不够可以用零矩阵占位
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_PHE = np.append(Data_PHE, temp).reshape(-1, 3)

            if aa.get_resname() == "TRP":
                '''计算Atoms 和 vector'''
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[11].get_vector()  # atom: CD2
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[12].get_vector()  # atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                S_Atom1 = np.array([-2.113650,   -0.895615,   -0.068390]) # atom：CA
                S_Atom2 = np.array([-1.238209,   -0.414289,   -1.257396]) # atom：CB
                S_Atom3 = np.array([0.004199,    0.327743,   -0.863152]) # atom: CG
                S_Atom4 = np.array([1.232942,   -0.236077,   -0.3574091])  # atom: CD2
                S_Atom5 = np.array([-3.296186,   -1.671994,   -0.450450]) # atom: N
                S_Atom6 = np.array([-2.625376,    0.306965,    0.713989])  # atom: C
                vector_1 = (Atom6 - Atom1)[0: 3]  # import vector for c-c(o), [-0.515, 1.208, 0.771]
                vector_2 = (Atom5 - Atom1)[0: 3]  # import vector for c-n, [-1.18, -0.781, -0.38]
                vector_3 = (Atom2 - Atom1)[0: 3]  # import vector for ch2 [0.88, 0.472, -1.189]
                vector_4 = (Atom4 - Atom3)[0: 3]  # import vector 1 for aromatic ring , [1.23, -0.564, 0.501]
                vector_5 = (Atom3 - Atom2)[0: 3]  # import vector 2 for aromatic ring, [1.241, 0.745, 0.395]
                '''计算不同的T,old version 
                # 苯环标准向量1 苯环标准向量2
                phStd1 = np.array([1.23, -0.564, 0.501])  # standard vector 1 for aromatic ring
                phStd2 = np.array([1.241, 0.745, 0.395])  # standard vector 2 for aromatic ring
                CH2td1 = np.array([0.88, 0.472, -1.189])  # standard vector for ch2
                CCO_Std = np.array([-0.515, 1.208, 0.771])  # standard vector for c-c(o)
                CN_Std = np.array([-1.18, -0.781, -0.38])  # standard vector for c-n

                T1 = RM.Calculate(phStd1, vector_4, phStd2, vector_5)  # transformation of aromatic ring
                T2 = RM.Calculate(phStd1, vector_4, CH2td1, vector_3)  # orientation of ch2
                T3 = RM.Calculate(CCO_Std, vector_1, CN_Std, vector_2)  # orientation of alpha carbon
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # import vector for c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # import vector for c-n
                S_vector_3 = (S_Atom2 - S_Atom1)  # import vector for ch2
                S_vector_4 = (S_Atom4 - S_Atom3)  # import vector 1 for aromatic ring 
                S_vector_5 = (S_Atom3 - S_Atom2)  # import vector 2 for aromatic ring
                T1 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)  # transformation of aromatic ring
                T2 = RM.Calculate(S_vector_4, vector_4, S_vector_3, vector_3)  # orientation of ch2
                T3 = RM.Calculate(S_vector_1, vector_1, S_vector_2, vector_2)  # orientation of alpha carbon
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]-1
                T_Set = [T1, T2, T3]
                # For amino acid which do not need 5 transforming matrix
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''Matrix performing'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_TRP = np.append(Data_TRP, temp).reshape(-1, 3)

            if aa.get_resname() == "THR":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG2
                Atom4 = atoms[4].get_vector()  # atom: OG1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                # Atom7 = atoms[5].get_vector()  # atom: O
                '''standard orientation'''
                S_Atom1 = np.array([-0.014306,    0.400817,   -0.577393])  # atom：CA
                S_Atom2 = np.array([1.053056,   -0.574984,   -0.010551])  # atom：CB
                S_Atom3 = np.array([2.460354,   -0.276346,   -0.518983])  # atom: CG2
                S_Atom4 = np.array([0.969982,   -0.468351,    1.412936])  # atom: OG1
                S_Atom5 = np.array([0.143683,    1.808834,   -0.262623])  # atom: N
                S_Atom6 = np.array([-1.393039,   -0.046916,   -0.105532])  # atom: C
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-1.38, -0.451, 0.466]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n, [0.156, 1.407, 0.32]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.067, -0.976, 0.568]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c, [1.406, 0.292, -0.517]
                vector_5 = (Atom4 - Atom2)[0:3]  # c-o, [-0.076, 0.116, 1.423]
                '''old version
                T1td1 = np.array([-1.38, -0.451, 0.466])
                T1td2 = np.array([0.156, 1.407, 0.32])
                T2td1 = np.array([1.067, -0.976, 0.568])
                CO_Std = np.array([-0.076, 0.116, 1.423])
                # Cc_Std = np.array([1.406, 0.292, -0.517])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CO_Std, vector_5)
                # T3 = RM.Calculate(CCO_Std, vector_1, CN_Std, vector_2)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c
                S_vector_5 = (S_Atom4 - S_Atom2)  # c-o
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_5, vector_5) 
                T_Set = [T1, T2]
                # T最多是5个，不够可以用零矩阵占位
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_THR = np.append(Data_THR, temp).reshape(-1, 3)

            if aa.get_resname() == "ALA":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                # Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: C
                Atom5 = atoms[0].get_vector()  # atom: N
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom4[:], Atom5[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([0.676950,    0.032249,    0.362285])  # atom：CA
                S_Atom2 = np.array([1.288317,    1.286925,   -0.271513])  # atom：CB
                # S_Atom3 = atoms[3].get_vector()  # atom: CG
                S_Atom4 = np.array([-0.797366,   -0.144218,    0.031983])  # atom: C
                S_Atom5 = np.array([1.396179,   -1.161570,   -0.072296])  # atom: N
                # S_Atom6 = atoms[10].get_vector()  # atom: C

                vector_1 = (Atom4 - Atom1)[0:3]  # 键轴方向c-c(o), [1.474, -0.176, 0.33]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.719, -1.194, 0.434]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-0.611, 1.255, 0.634]
                '''old version'''
                # T1td1 = np.array([-0.719, -1.194, 0.434])
                # T1td2 = np.array([-0.611, 1.255, 0.634])
                # T1 = RM.Calculate(T1td2, vector_3, T1td1, vector_2)
                S_vector_1 = (S_Atom4 - S_Atom1)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                T1 = RM.Calculate(S_vector_3, vector_3, S_vector_2, vector_2)

                T_Set = [T1]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom1[0]
                # Y = Atom1[1]
                Z = Atom1[2]
                # T最多是5个，不够可以用零矩阵占位
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_ALA = np.append(Data_ALA, temp).reshape(-1, 3)

            if aa.get_resname() == "ILE":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG2
                Atom4 = atoms[4].get_vector()  # atom: CG1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # C
                Atom7 = atoms[5].get_vector()  # atom: CD
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([0.625228,   -0.508605,   -0.632754])  # atom：CA
                S_Atom2 = np.array([-0.797304,   -0.573033,    0.006618])  # atom：CB
                S_Atom3 = np.array([-0.720466,   -0.619522,    1.536719])  # atom: CG2
                S_Atom4 = np.array([-1.760482,    0.501941,   -0.528087])  # atom: CG1
                S_Atom5 = np.array([1.297313,   -1.796603,   -0.468075])  # atom: N
                S_Atom6 = np.array([1.527308,    0.576302,   -0.060460])  # C
                S_Atom7 = np.array([-3.228612,    0.243851,   -0.174319])  # atom: CD
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.898, 1.084, 0.581]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.674, -1.288, 0.162]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.424, -0.07, 0.635]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [0.072, -0.046, 1.53]
                vector_5 = (Atom4 - Atom2)[0:3]  # c-c 2, [-0.964, 1.071, -0.539]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-c 3, [-1.466, -0.242, 0.374]
                '''old version#还有提升的空间，可用4个甚至更多T来描述。
                T1td1 = np.array([0.898, 1.084, 0.581])
                T1td2 = np.array([0.674, -1.288, 0.162])
                T2td1 = np.array([-1.424, -0.07, 0.635])
                CC1_Std = np.array([0.072, -0.046, 1.53])
                CC2_Std = np.array([-0.964, 1.071, -0.539])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CC2_Std, vector_5, CC1_Std, vector_4)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1
                S_vector_5 = (S_Atom4 - S_Atom2)  # c-c 2
                S_vector_6 = (S_Atom7 - S_Atom4)  # c-c 3
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_5, vector_5, S_vector_4, vector_4)
                T_Set = [T1, T2, T3]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_ILE = np.append(Data_ILE, temp).reshape(-1, 3)

            if aa.get_resname() == "LEU":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # C
                Atom7 = atoms[5].get_vector()  # atom: CD2
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([0.770509,    0.661364,    0.141644])  # atom：CA
                S_Atom2 = np.array([-0.548286,    0.456663,   -0.639337])  # atom：CB
                S_Atom3 = np.array([-1.693648,   -0.266233,    0.094971])  # atom: CG
                S_Atom4 = np.array([-2.270513,    0.580806,    1.238406])  # atom: CD1
                S_Atom5 = np.array([1.496559,    1.809670,   -0.400514])  # atom: N
                S_Atom6 = np.array([1.710057,   -0.536330,    0.107059])  # C
                S_Atom7 = np.array([-2.798827,   -0.637075,   -0.902916])  # atom: CD2
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.937, -1.201, -0.028]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.728, 1.143, -0.551]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.319, -0.207, -0.78]
                vector_4 = (Atom4 - Atom3)[0:3]  # c-c 1, [-0.575, 0.838, 1.152]
                vector_5 = (Atom7 - Atom3)[0:3]  # c-c 2, [-1.107, -0.362, -0.999]
                '''old version
                T1td1 = np.array([0.937, -1.201, -0.028])
                T1td2 = np.array([0.728, 1.143, -0.551])
                T2td1 = np.array([-1.319, -0.207, -0.78])
                CC1_Std = np.array([-0.575, 0.838, 1.152])
                CC2_Std = np.array([-1.107, -0.362, -0.999])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CC2_Std, vector_5, CC1_Std, vector_4)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用
                S_vector_4 = (S_Atom4 - S_Atom3)  # c-c 1
                S_vector_5 = (S_Atom7 - S_Atom3)  # c-c 2
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_5, vector_5, S_vector_4, vector_4)
                T_Set = [T1, T2, T3]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_LEU = np.append(Data_LEU, temp).reshape(-1, 3)

            if aa.get_resname() == "VAL":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG1
                Atom4 = atoms[4].get_vector()  # atom: CG2
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:],
                                  Atom6[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([-0.169937,    0.515682,   -0.636527])  # atom：CA
                S_Atom2 = np.array([1.280367,    0.380676,   -0.077326])  # atom：CB
                S_Atom3 = np.array([1.309693,    0.414236,    1.453780])  # atom: CG1
                S_Atom4 = np.array([2.058252,   -0.810514,   -0.647361])  # atom: CG2
                S_Atom5 = np.array([-0.653621,    1.878872,   -0.422876])  # atom: N
                S_Atom6 = np.array([-1.174014,   -0.449758,   -0.022226])  # C
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-1.004, -0.966, 0.614]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.484, 1.363, 0.212]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.45, -0.135, 0.559]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [0.03, 0.032, 1.532]
                vector_5 = (Atom4 - Atom2)[0:3]  # c-c 2, [-.778, -1.19, -0.571]
                '''old version
                T1td1 = np.array([-1.004, -0.966, 0.614])
                T1td2 = np.array([-0.484, 1.363, 0.212])
                T2td1 = np.array([1.45, -0.135, 0.559])
                CC1_Std = np.array([0.03, 0.032, 1.532])
                # CC2_Std = np.array([-1.107, -0.362, -0.999])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                # T3 = RM.Calculate(CC2_Std, vector_5, CC1_Std, vector_4)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1
                S_vector_5 = (S_Atom4 - S_Atom2)  # c-c 2
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)

                T_Set = [T1, T2]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_VAL = np.append(Data_VAL, temp).reshape(-1, 3)

            if aa.get_resname() == "LYS":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: CE
                Atom8 = atoms[6].get_vector()  # atom: NZ
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([1.785591,    0.152509,    0.3992635])  # atom：CA
                S_Atom2 = np.array([1.075474,    1.488812,    0.037818])  # atom：CB
                S_Atom3 = np.array([-0.410085,    1.612429,    0.444633])  # atom: CG
                S_Atom4 = np.array([-1.418226,    1.125598,   -0.632495])  # atom: CD
                S_Atom5 = np.array([3.068216,    0.092346,   -0.292257])  # atom: N
                S_Atom6 = np.array([0.828002,   -0.993279,    0.085400])  # atom: C
                S_Atom7 = np.array([-2.670424,    0.397817,   -0.1111775])  # atom: CE
                S_Atom8 = np.array([-2.454859,   -1.066978,    0.021610])  # atom: NZ
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.709, 1.304, -0.348]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.877, -1.1, -0.395]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.386, -0.174, -0.659]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [1.161, 0.719, 0.693]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-c 2, [-1.345, 0.283, 0.685]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-c 3, [1.135, 0.715, 0.728]
                vector_7 = (Atom7 - Atom8)[0:3]  # c-n 2, [-1.348, 0.236, 0.637]
                '''old version
                T1td1 = np.array([-0.709, 1.304, -0.348])
                T1td2 = np.array([-0.877, -1.1, -0.395])
                T2td1 = np.array([1.386, -0.174, -0.659])
                CC1_Std = np.array([1.161, 0.719, 0.693])
                CC2_Std = np.array([-1.345, 0.283, 0.685])
                CC3_Std = np.array([1.135, 0.715, 0.728])
                CN_Std = np.array([-1.348, 0.236, 0.637])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CC2_Std, vector_5, CC1_Std, vector_4)
                T4 = RM.Calculate(CC3_Std, vector_6, CC2_Std, vector_5)
                T5 = RM.Calculate(CC3_Std, vector_6, CN_Std, vector_7)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用,
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1
                S_vector_5 = (S_Atom3 - S_Atom4)  # c-c 2
                S_vector_6 = (S_Atom7 - S_Atom4)  # c-c 3
                S_vector_7 = (S_Atom7 - S_Atom8)  # c-n 2
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_5, vector_5, S_vector_4, vector_4)
                T4 = RM.Calculate(S_vector_6, vector_6, S_vector_5, vector_5)
                T5 = RM.Calculate(S_vector_6, vector_6, S_vector_7, vector_7)

                T_Set = [T1, T2, T3, T4, T5]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_LYS = np.append(Data_LYS, temp).reshape(-1, 3)

            ###########################

            if aa.get_resname() == "ARG":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[9].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: NE
                Atom8 = atoms[6].get_vector()  # atom: CZ
                Atom9 = atoms[7].get_vector()  # atom: NH1
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:], Atom9[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([2.562007,   -0.515202,   -0.477131])  # atom：CA
                S_Atom2 = np.array([1.161198,   -0.917117,    0.054729])  # atom：CB
                S_Atom3 = np.array([-0.019817,   -0.075583,   -0.442142])  # atom: CG
                S_Atom4 = np.array([-1.344587,   -0.639998,    0.068853])  # atom: CD
                S_Atom5 = np.array([3.488172,   -1.613497,   -0.238367])  # atom: N
                S_Atom6 = np.array([3.115644,    0.746547,    0.180016])  # atom: C
                S_Atom7 = np.array([-2.470912,    0.188447,   -0.409034])  # atom: NE
                S_Atom8 = np.array([-3.725467,    0.110992,    0.024296])  # atom: CZ
                S_Atom9 = np.array([-4.100196,   -0.878145,    0.844395])  # atom: NH1
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.75, 1.296, 0.281]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.874, -1.107, 0.281]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.36, -0.131, 0.72]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [-1.178, 0,708, 0.043]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-c 2, [1.314, 0.239, -0.742]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-n 1, [-1.115, 0.684, -0.742]
                vector_7 = (Atom7 - Atom8)[0:3]  # c-n 2, [1.318, 0.31, -0.404]
                vector_8 = (Atom9 - Atom8)[0:3]  # c-n 3, [-0.358, -1.334, -0.269]
                '''old version
                T1td1 = np.array([0.75, 1.296, 0.281])
                T1td2 = np.array([0.874, -1.107, 0.281])
                T2td1 = np.array([-1.36, -0.131, 0.72])
                CC1_Std = np.array([-1.178, 0.708, 0.043])
                CC2_Std = np.array([1.314, 0.239, -0.742])
                CN1_Std = np.array([-1.115, 0.684, -0.742])
                CN2_Std = np.array([1.318, 0.31, -0.404])
                CN3_Std = np.array([-0.358, -1.334, -0.269])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(CC1_Std, vector_4, T2td1, vector_3)
                T3 = RM.Calculate(CC2_Std, vector_5, CC1_Std, vector_4)
                T4 = RM.Calculate(CN1_Std, vector_6, CC2_Std, vector_5)
                T5 = RM.Calculate(CN3_Std, vector_8, CN2_Std, vector_7)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o), [0.75, 1.296, 0.281]
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用, [0.874, -1.107, 0.281]
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用, [-1.36, -0.131, 0.72]
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1, [-1.178, 0,708, 0.043]
                S_vector_5 = (S_Atom3 - S_Atom4)  # c-c 2, [1.314, 0.239, -0.742]
                S_vector_6 = (S_Atom7 - S_Atom4)  # c-n 1, [-1.115, 0.684, -0.742]
                S_vector_7 = (S_Atom7 - S_Atom8)  # c-n 2, [1.318, 0.31, -0.404]
                S_vector_8 = (S_Atom9 - S_Atom8)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_4, vector_4, S_vector_3, vector_3)
                T3 = RM.Calculate(S_vector_5, vector_5, S_vector_4, vector_4)
                T4 = RM.Calculate(S_vector_6, vector_6, S_vector_5, vector_5)
                T5 = RM.Calculate(S_vector_8, vector_8, S_vector_7, vector_7)
                T_Set = [T1, T2, T3, T4, T5]   #, T2, T3, T4, T5]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_ARG = np.append(Data_ARG, temp).reshape(-1, 3)

            if aa.get_resname() == "ASN":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: OD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: ND2
                # Atom8 = atoms[6].get_vector()  # atom: NZ
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([0.742603,    0.690993,    0.164065])  # atom：CA
                S_Atom2 = np.array([-0.573188,    0.784537,   -0.627202])  # atom：CB
                S_Atom3 = np.array([-1.766833,    0.191004,    0.115614])  # atom: CG
                S_Atom4 = np.array([-1.933497,    0.365950,    1.310752])  # atom: OD1
                S_Atom5 = np.array([1.693484,    1.698166,   -0.307639])  # atom: N
                S_Atom6 = np.array([1.410997,   -0.672376,    0.077368])  # atom: C
                S_Atom7 = np.array([-2.674230,   -0.470093,   -0.659954])  # atom: ND2
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.668, -1.363, -0.087]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.951, 1.007, -0.472]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.316, 0.094, -0.791]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [-1.193, -0.595, 0.743]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-o 2, [0.167, -0.173, -1.195]
                vector_6 = (Atom3 - Atom7)[0:3]  # c-n 3, [0.908, 0.66, 0.776]
                '''old version
                T1td1 = np.array([0.668, -1.363, -0.087])
                T1td2 = np.array([0.951, 1.007, -0.472])
                T2td1 = np.array([-1.316, 0.094, -0.791])
                CC1_Std = np.array([-1.193, -0.595, 0.743])
                CO_Std = np.array([0.167, -0.173, -1.195])
                CN_Std = np.array([0.908, 0.66, 0.776])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CO_Std, vector_5, CN_Std, vector_6)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o)
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1
                S_vector_5 = (S_Atom3 - S_Atom4)  # c-o 2
                S_vector_6 = (S_Atom3 - S_Atom7)  # c-n 3
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_5, vector_5, S_vector_6, vector_6)
                T_Set = [T1, T2, T3]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_ASN = np.append(Data_ASN, temp).reshape(-1, 3)

            if aa.get_resname() == "ASP":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: OD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OD2
                # Atom8 = atoms[6].get_vector()  # atom: NZ
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([0.593896,   -0.667680,   -0.370859])  # atom：CA
                S_Atom2 = np.array([-0.643946,   -0.715178,    0.577965])  # atom：CB
                S_Atom3 = np.array([-1.918663,   -0.017479,    0.097179])  # atom: CG
                S_Atom4 = np.array([-1.807928,    1.265738,   -0.309363])  # atom: OD1
                S_Atom5 = np.array([1.245054,   -1.965158,   -0.387587])  # atom: N
                S_Atom6 = np.array([1.541560,    0.472290,    0.019045])  # atom: C
                S_Atom7 = np.array([-2.985082,   -0.568616,   0.091646])  # atom: OD2
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [0.947, 1.141, 0.389]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.651, -1.297, -0.017]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.238,-0.047 ,0.948]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [-1.275, 0.698, -0.482]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-o 2, [-0.112, -1.284, 0.403]
                vector_6 = (Atom3 - Atom7)[0:3]  # c-o 3, [1.066, 0.552, 0.008]
                '''old version
                T1td1 = np.array([0.947, 1.141, 0.389])
                T1td2 = np.array([0.651, -1.297, -0.017])
                T2td1 = np.array([-1.238, -0.047, 0.948])
                CC1_Std = np.array([-1.275, 0.698, -0.482])
                #还有提升的空间
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o), [0.947, 1.141, 0.389]
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用, [0.651, -1.297, -0.017]
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用, [-1.238,-0.047 ,0.948]
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1, [-1.275, 0.698, -0.482]
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T_Set = [T1, T2]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_ASP = np.append(Data_ASP, temp).reshape(-1, 3)

            if aa.get_resname() == "CYS":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: s
                Atom4 = atoms[4].get_vector()  # atom: C
                Atom5 = atoms[0].get_vector()  # atom: N
                avg_atoms = np.c_[
                    Atom1[:], Atom2[:], Atom3[:],
                    Atom4[:], Atom5[:]
                    ].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([-0.405904,    0.679787,   -0.447486])  # atom：CA
                S_Atom2 = np.array([0.923658,    0.924040,    0.289004])  # atom：CB
                S_Atom3 = np.array([2.147072,   -0.412974,   -0.086756])  # atom: s
                S_Atom4 = np.array([-1.096538,   -0.562798,    0.108867])  # atom: C
                S_Atom5 = np.array([-1.365151,    1.776402,   -0.312738])  # atom: N
                vector_1 = (Atom4 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.691, -1.243, 0.556]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.96, 1.096, 0.133]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.329, 0.245, 0.737]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-s 1, [1.224, -1.337, -0.375]
                '''old version
                T1td1 = np.array([-0.691, -1.243, 0.556])
                T1td2 = np.array([-0.96, 1.096, 0.133])
                T2td1 = np.array([1.329, 0.245, 0.737])
                CS_Std = np.array([1.224, -1.337, -0.375])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CS_Std, vector_4)'''
                S_vector_1 = (S_Atom4 - S_Atom1)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                S_vector_4 = (S_Atom3 - S_Atom2)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T_Set = [T1, T2]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode, z_init, z_end)
                Data_CYS = np.append(Data_CYS, temp).reshape(-1, 3)

            if aa.get_resname() == "GLN":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OE1
                Atom8 = atoms[6].get_vector()  # atom: NE2
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:]].mean(-1)
                '''standart orientation'''
                S_Atom1 = np.array([-1.289377,   -0.637556,    0.338258])  # atom：CA
                S_Atom2 = np.array([ 0.131356,   -0.634882,   -0.262171])  # atom：CB
                S_Atom3 = np.array([1.143626,    0.286403,    0.424951])  # atom: CG
                S_Atom4 = np.array([2.575856,   -0.054985,    0.020918])  # atom: CD
                S_Atom5 = np.array([-2.003582,   -1.823855,   -0.12208])  # atom: N
                S_Atom6 = np.array([-2.140482,    0.577419,   -0.007575])  # atom: C
                S_Atom7 = np.array([2.972614,   -1.202034,   -0.092498])  # atom: OE1
                S_Atom8 = np.array([3.398041,    1.017907,   -0.177176])  # atom: NE2

                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.842, 1.215, -0.348]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.706, -1.187, -0.459]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.43, 0.002, -0.599]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [1.013, 0.914, 0.695]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-c 2, [-1.432, 0.335, 0.412]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-o 1, [0.4, -1.146, -0.117]
                vector_7 = (Atom8 - Atom4)[0:3]  # c-n 2, [0.727, 1.076, -0.2]
                '''old version
                T1td1 = np.array([-0.842, 1.215, -0.348])
                T1td2 = np.array([-0.706, -1.187, -0.459])
                T2td1 = np.array([1.43, 0.002, -0.599])
                CC1_Std = np.array([1.013, 0.914, 0.695])
                CC2_Std = np.array([-1.432, 0.335, 0.412])
                CO_Std = np.array([0.4, -1.146, -0.117])
                CN_Std = np.array([0.727, 1.076, -0.2])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CC1_Std, vector_4, CC2_Std, vector_5)
                T4 = RM.Calculate(CO_Std, vector_6, CN_Std, vector_7)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                S_vector_4 = (S_Atom3 - S_Atom2)
                S_vector_5 = (S_Atom3 - S_Atom4)
                S_vector_6 = (S_Atom7 - S_Atom4)
                S_vector_7 = (S_Atom8 - S_Atom4)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)
                T4 = RM.Calculate(S_vector_6, vector_6, S_vector_7, vector_7)
                T_Set = [T1, T2, T3, T4]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_GLN = np.append(Data_GLN, temp).reshape(-1, 3)

            if aa.get_resname() == "GLU":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OE1
                Atom8 = atoms[6].get_vector()  # atom: OE2
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:]].mean(-1)
                '''standard orientation'''
                S_Atom1 = np.array([1.177936,    0.602759,    0.330390])  # atom：CA
                S_Atom2 = np.array([-0.207478,    0.518756,   -0.344511])  # atom：CB
                S_Atom3 = np.array([-1.211321,   -0.472563,    0.246027])  # atom: CG
                S_Atom4 = np.array([-2.713205,   -0.071682,   -0.005201])  # atom: CD
                S_Atom5 = np.array([1.831500,    1.869342,   -0.018289])  # atom: N
                S_Atom6 = np.array([2.146297,   -0.517642,   -0.002824])  # atom: C
                S_Atom7 = np.array([-3.518526,   -1.031001,   -0.024000])  # atom: OE1
                S_Atom8 = np.array([ -2.945161,    1.157874,   -0.117145])  # atom: OE2
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.84, -1.185, 0.466]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.71, 1.221, 0.357]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.425, 0.049, 0.599]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [0.999, -0.988, -0.613]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-c 2, [-1.413, -0.253, -0.613]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-o 1, [0.923, -0.854, -0.105]
                vector_7 = (Atom8 - Atom4)[0:3]  # c-o 2, [0.425, 1.078, 0.484]
                '''old version
                T1td1 = np.array([-0.84, -1.185, 0.466])
                T1td2 = np.array([-0.71, 1.221, 0.357])
                T2td1 = np.array([1.425, 0.049, 0.599])
                CC1_Std = np.array([0.999, -0.988, -0.613])
                CC2_Std = np.array([-1.413, -0.253, -0.613])
                CO_Std = np.array([0.923, -0.854, -0.105])
                CO2_Std = np.array([0.425, 1.078, 0.484])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T2 = RM.Calculate(T2td1, vector_3, CC1_Std, vector_4)
                T3 = RM.Calculate(CC1_Std, vector_4, CC2_Std, vector_5)
                T4 = RM.Calculate(CO_Std, vector_6, CO2_Std, vector_7)
                '''
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o), [-0.84, -1.185, 0.466]
                S_vector_2 = (S_Atom5 - S_Atom1)  # c-n,暂时没用, [-0.71, 1.221, 0.357]
                S_vector_3 = (S_Atom2 - S_Atom1)  # for the future,暂时没用, [1.425, 0.049, 0.599]
                S_vector_4 = (S_Atom3 - S_Atom2)  # c-c 1, [0.999, -0.988, -0.613]
                S_vector_5 = (S_Atom3 - S_Atom4)  # c-c 2, [-1.413, -0.253, -0.613]
                S_vector_6 = (S_Atom7 - S_Atom4)  # c-o 1, [0.923, -0.854, -0.105]
                S_vector_7 = (S_Atom8 - S_Atom4)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)
                T4 = RM.Calculate(S_vector_6, vector_6, S_vector_7, vector_7)
                T_Set = [T1, T2, T3, T4]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_GLU = np.append(Data_GLU, temp).reshape(-1, 3)

            if aa.get_resname() == "GLY":
                Atom1 = atoms[1].get_vector()  # atom: CA
                # Atom2 = atoms[2].get_vector()  # atom：CB
                # Atom3 = atoms[3].get_vector()  # atom: CG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[2].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: NE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1

                avg_atoms = np.c_[Atom1[:], Atom5[:], Atom6[:], ].mean(-1)

                '''standard orientation'''
                S_Atom1 = np.array([-0.728530,   -0.709418,    0.131228])  # atom: CA
                S_Atom5 = np.array([-1.886203,    0.102096,   -0.210543])  # atom: N
                S_Atom6 = np.array([0.539175,    0.105452,    0.017471])  # atom: C
                '''old version
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [1.268, 0.814, -0.114]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-1.157, 0.811, -0.342]
                T1td1 = np.array([1.268, 0.814, -0.114])
                T1td2 = np.array([-1.157, 0.811, -0.342])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)'''
                vector_1 = (Atom6 - Atom1)[0:3]
                vector_2 = (Atom5 - Atom1)[0:3]
                S_vector_1 = (S_Atom6 - S_Atom1)  # 键轴方向c-c(o), [1.268, 0.814, -0.114]
                S_vector_2 = (S_Atom5 - S_Atom1)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T_Set = [T1]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom1[0]
                # Y = Atom1[1]
                Z = Atom1[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_GLY = np.append(Data_GLY, temp).reshape(-1, 3)

            if aa.get_resname() == "MET":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: SD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: CE
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:]].mean(-1)

                '''standard orientation'''
                S_Atom1 = np.array([-1.458483,   -0.595563,    0.406257])  # atom：CA
                S_Atom2 = np.array([-0.078961,   -0.790572,   -0.266162])  # atom：CB
                S_Atom3 = np.array([1.062591,    0.021109,    0.344617])  # atom: CG
                S_Atom4 = np.array([2.647310,   -0.441116,   -0.461079])  # atom: SD
                S_Atom5 = np.array([-2.315101,   -1.736402,    0.096081])  # atom: N
                S_Atom6 = np.array([-2.195780,    0.668459,   -0.017972])  # atom: C
                S_Atom7 = np.array([ 3.768422,    0.697215,    0.418044])  # atom: CE
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [-0.738, 1.263, -0.426]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [-0.857, -1.141, -0.308]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [1.379, -0.195, -0.672]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-c 1, [1.142, 0.81, 0.613]
                vector_5 = (Atom3 - Atom4)[0:3]  # c-s 1, [-1.583, 0.456, 0.812]
                vector_6 = (Atom7 - Atom4)[0:3]  # c-s 2, [1.123, 1.133, 0.884]
                # T1td1 = np.array([-0.738, 1.263, -0.426])
                # T1td2 = np.array([-0.857, -1.141, -0.308])
                # T2td1 = np.array([1.379, -0.195, -0.672])
                # CC1_Std = np.array([1.142, 0.81, 0.613])
                # CS_Std = np.array([-1.583, 0.456, 0.812])
                # # CS2_Std = np.array([1.123, 1.133, 0.884])
                S_vector_1 = (S_Atom6 - S_Atom1)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                S_vector_4 = (S_Atom3 - S_Atom2)
                S_vector_5 = (S_Atom3 - S_Atom4)
                S_vector_6 = (S_Atom7 - S_Atom4)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)
                T3 = RM.Calculate(S_vector_4, vector_4, S_vector_5, vector_5)

                T_Set = [T1, T2, T3]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_MET = np.append(Data_MET, temp).reshape(-1, 3)

            if aa.get_resname() == "SER":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: OG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[4].get_vector()  # atom: C
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom5[:], Atom6[:], ].mean(-1)

                '''standard orientation'''
                S_Atom1 = np.array([-0.256214,    0.574018,    0.444502])  # atom：CA
                S_Atom2 = np.array([-1.549050,    0.061470,   -0.186637])  # atom：CB
                S_Atom3 = np.array([-1.653139,   -1.329190,    0.104625])  # atom: OG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                S_Atom5 = np.array([-0.139612,    2.008145,    0.180827])  # atom: N
                S_Atom6 = np.array([0.930369,   -0.215298,   -0.111776])  # atom: C
                vector_1 = (Atom6 - Atom1)[0:3]  # 键轴方向c-c(o), [1.187, -0.789, -0.556]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n,暂时没用, [0.116, 1.434, -0.263]
                vector_3 = (Atom2 - Atom1)[0:3]  # for the future,暂时没用, [-1.293, -0.513, -0.631]
                vector_4 = (Atom3 - Atom2)[0:3]  # c-o 1, [-0.105, -1.39, 0.294]
                # T1td1 = np.array([1.187, -0.789, -0.556])
                # T1td2 = np.array([0.116, 1.434, -0.263])
                # T2td1 = np.array([-1.293, -0.513, -0.631])
                # CO_Std = np.array([-0.105, -1.39, 0.294])
                S_vector_1 = (S_Atom6 - S_Atom1)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                S_vector_4 = (S_Atom3 - S_Atom2)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_1, vector_1)
                T2 = RM.Calculate(S_vector_3, vector_3, S_vector_4, vector_4)

                T_Set = [T1, T2]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom3[0]
                # Y = Atom3[1]
                Z = Atom3[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_SER = np.append(Data_SER, temp).reshape(-1, 3)

            if aa.get_resname() == "PRO":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CD
                Atom2 = atoms[4].get_vector()  # atom：CA
                # Atom3 = atoms[3].get_vector()  # atom: CG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: NE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom5[:], Atom6[:], ].mean(-1)

                '''standard orientation'''
                S_Atom1 = np.array([-1.932552,    0.794676,   -0.248481])  # atom：CD
                S_Atom2 = np.array([ 0.089722,   -0.161132,    0.605864])  # atom：CA
                S_Atom5 = np.array([-0.798135,    1.003972,    0.665051])  # atom: N
                S_Atom6 = np.array([1.441581,    0.155540,   -0.007587])  # atom: C

                vector_1 = (Atom6 - Atom2)[0:3]  # 键轴方向c-c(o), [1.352, 0.315, -0.613]
                vector_2 = (Atom5 - Atom1)[0:3]  # c-n, 环1, [1.136, 0.21, 0.912]
                vector_3 = (Atom2 - Atom1)[0:3]  # 环2, [2.024, -0.955, 0.853]
                # T1td1 = np.array([1.136, 0.21, 0.912])
                # T1td2 = np.array([2.024, -0.955, 0.853])
                S_vector_1 = (S_Atom6 - S_Atom2)
                S_vector_2 = (S_Atom5 - S_Atom1)
                S_vector_3 = (S_Atom2 - S_Atom1)
                T1 = RM.Calculate(S_vector_2, vector_2, S_vector_3, vector_3)
                T_Set = [T1]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom1[0]
                # Y = Atom1[1]
                Z = Atom1[2]
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_PRO = np.append(Data_PRO, temp).reshape(-1, 3)

            if aa.get_resname() == "HEM":
                Atom1 = atoms[0].get_vector()  # atom：FE
                Atom2 = atoms[1].get_vector()  # atom：NA
                Atom3 = atoms[2].get_vector()  # atom: NB
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], ].mean(-1)

                vector_1 = (Atom2 - Atom1)[0:3]  # Fe-Na, [?]
                vector_2 = (Atom3 - Atom1)[0:3]  # Fe-Nb, [?]

                T1td1 = np.array([1.206, 1.585, 0.022])
                T1td2 = np.array([-1.587, 1.207, -0.047])
                T1 = RM.Calculate(T1td2, vector_2, T1td1, vector_1)
                T_Set = [T1]
                # 计算均匀场模型所用坐标：
                # X = 1
                # Y = 1
                # Z = 1
                # 计算SERS和扩散场所用坐标：
                # X = Atom1[0]
                # Y = Atom1[1]
                Z = Atom1[2] - 10  #此处因激发线不同需稍作修改（8）
                # print(T_Set)
                while len(T_Set) < 5:
                    T_Set.append(T_zeros)
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_HEM = np.append(Data_HEM, temp).reshape(-1, 3)

            if aa.get_resname() == "DA ":
                Atom1 = atoms[0].get_vector()  # atom：P: 3.792	-1.362	-0.085
                Atom2 = atoms[1].get_vector()  # atom：O1P: 5.203	-1.004	-0.402
                Atom3 = atoms[3].get_vector()  # atom: O5': 2.769	-0.141	-0.007
                Atom4 = atoms[4].get_vector()  # atom: C5': 2.763	0.908	-1.019
                Atom5 = atoms[5].get_vector()  # atom: C4': 1.608	1.845	-0.747
                Atom6 = atoms[6].get_vector()  # atom: O4': 0.361	1.15	-0.966
                Atom7 = atoms[7].get_vector()  # atom: C1': -0.546	1.436	0.109
                Atom8 = atoms[8].get_vector()  # tyr_atom: N9: -1.46	0.322	0.261
                Atom9 = atoms[9].get_vector()  # atom: C8: -1.168	-0.951	0.705
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:], Atom9[:]].mean(-1)
                vector_1 = (Atom2 - Atom1)[0:3]
                vector_2 = (Atom3 - Atom1)[0:3]
                vector_3 = (Atom4 - Atom3)[0:3]
                vector_4 = (Atom5 - Atom4)[0:3]
                vector_5 = (Atom6 - Atom5)[0:3]
                vector_6 = (Atom7 - Atom6)[0:3]
                vector_7 = (Atom8 - Atom7)[0:3]
                vector_8 = (Atom9 - Atom8)[0:3]

                v_1 = np.array([1.411, 0.358, -0.317])  #new
                v_2 = np.array([-1.023, 1.221, 0.078])  #new
                v_3 = np.array([-0.006, 1.049, -1.012])  #new
                v_4 = np.array([-1.155, 0.937, 0.272])  #new
                v_5 = np.array([-1.247, -0.695, -0.219])  #new
                v_6 = np.array([-0.907, 0.286, 1.075])  #new
                v_7 = np.array([-0.914, -1.114, 0.152])  #new
                v_8 = np.array([0.292, -1.273, 0.444])  #new

                T1 = RM.Calculate(v_1, vector_1, v_2, vector_2)
                T2 = RM.Calculate(v_3, vector_3, v_4, vector_4)
                T3 = RM.Calculate(v_4, vector_4, v_5, vector_5)
                T4 = RM.Calculate(v_5, vector_5, v_6, vector_6)
                T5 = RM.Calculate(v_7, vector_7, v_8, vector_8)

                T_Set = [T1, T2, T3, T4, T5]
                # print(T_Set)
                Z = Atom7[2]
                '''Matrix performing'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_DA = np.append(Data_DA, temp, axis=0).reshape(-1, 3)

            if aa.get_resname() == "DG ":
                Atom1 = atoms[0].get_vector()  # atom: P: -3.309	-1.591	0.144
                Atom2 = atoms[1].get_vector()  # atom: O1P: -3.229	-2.347	-1.131
                Atom3 = atoms[3].get_vector()  # atom: O5': -3.131	-0.03	-0.133
                Atom4 = atoms[4].get_vector()  # atom: C5': -3.114	0.993	0.904
                Atom5 = atoms[5].get_vector()  # atom: C4': -1.968	1.951	0.656
                Atom6 = atoms[6].get_vector()  # atom: O4': -0.715	1.256	0.833
                Atom7 = atoms[7].get_vector()  # atom: C1: 0.149	1.487	-0.291
                Atom8 = atoms[8].get_vector()  # atom: N9: 0.96	0.302	-0.513
                Atom9 = atoms[9].get_vector()  # atom: C8: 0.549	-0.917	-1.029
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:], Atom9[:]].mean(-1)

                vector_1 = (Atom2 - Atom1)[0:3]
                vector_2 = (Atom3 - Atom1)[0:3]
                vector_3 = (Atom4 - Atom3)[0:3]
                vector_4 = (Atom5 - Atom4)[0:3]
                vector_5 = (Atom6 - Atom5)[0:3]
                vector_6 = (Atom7 - Atom6)[0:3]
                vector_7 = (Atom8 - Atom7)[0:3]
                vector_8 = (Atom9 - Atom8)[0:3]
                v_1 = np.array([0.08, -0.756, -1.275])  # new
                v_2 = np.array([0.178, 1.561, -0.277])  # new
                v_3 = np.array([0.017, 1.023, 1.037])  # new
                v_4 = np.array([1.146, 0.958, -0.248])  # new
                v_5 = np.array([1.253, -0.695, 0.177])  # new
                v_6 = np.array([0.864, 0.231, -1.124])  # new
                v_7 = np.array([0.811, -1.185, -0.222])  # new
                v_8 = np.array([-0.411, -1.219, -0.516])  # new

                T1 = RM.Calculate(v_1, vector_1, v_2, vector_2)
                T2 = RM.Calculate(v_3, vector_3, v_4, vector_4)
                T3 = RM.Calculate(v_4, vector_4, v_5, vector_5)
                T4 = RM.Calculate(v_5, vector_5, v_6, vector_6)
                T5 = RM.Calculate(v_7, vector_7, v_8, vector_8)

                T_Set = [T1, T2, T3, T4, T5]
                Z = Atom7[2]
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_DG = np.append(Data_DG, temp).reshape(-1, 3)

            if aa.get_resname() == "DC ":
                Atom1 = atoms[0].get_vector()  # atom：P : 3.114	-1.349	0.091
                Atom2 = atoms[1].get_vector()  # atom：O1P: 4.566	-1.019	0.12
                Atom3 = atoms[3].get_vector()  # atom: O5': 2.116	-0.109	0.011
                Atom4 = atoms[4].get_vector()  # atom: C5': 2.263	0.94	-0.99
                Atom5 = atoms[5].get_vector()  # atom: C4': 1.142	1.942	-0.815
                Atom6 = atoms[6].get_vector()  # atom: O4': -0.124	1.323	-1.137
                Atom7 = atoms[7].get_vector()  # atom: C1': -1.073	1.547	-0.085
                Atom8 = atoms[8].get_vector()  # atom: N1 : -1.954	0.378	-0.001
                Atom9 = atoms[9].get_vector()  # atom: C6 : -1.406	-0.867	0.142
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:], Atom9[:]].mean(-1)

                vector_1 = (Atom2 - Atom1)[0:3]
                vector_2 = (Atom3 - Atom1)[0:3]
                vector_3 = (Atom4 - Atom3)[0:3]
                vector_4 = (Atom5 - Atom4)[0:3]
                vector_5 = (Atom6 - Atom5)[0:3]
                vector_6 = (Atom7 - Atom6)[0:3]
                vector_7 = (Atom8 - Atom7)[0:3]
                vector_8 = (Atom9 - Atom8)[0:3]
                v_1 = np.array([1.452, 0.33, 0.029])  # new
                v_2 = np.array([-0.998, 1.24, -0.08])  # new
                v_3 = np.array([0.147, 1.049, -1.001])  # new
                v_4 = np.array([-1.121, 1.002, 0.175])  # new
                v_5 = np.array([-1.266, -0.619, -0.322])  # new
                v_6 = np.array([-0.949, 0.224, 1.052])  # new
                v_7 = np.array([-0.881, -1.169, 0.084])  # new
                v_8 = np.array([0.548, -1.245, 0.143])  # new

                T1 = RM.Calculate(v_1, vector_1, v_2, vector_2)
                T2 = RM.Calculate(v_3, vector_3, v_4, vector_4)
                T3 = RM.Calculate(v_4, vector_4, v_5, vector_5)
                T4 = RM.Calculate(v_5, vector_5, v_6, vector_6)
                T5 = RM.Calculate(v_7, vector_7, v_8, vector_8)

                T_Set = [T1, T2, T3, T4, T5]
                Z = Atom7[2]
                '''处理Matrix'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_DC = np.append(Data_DC, temp).reshape(-1, 3)

            if aa.get_resname() == "DT ":
                Atom1 = atoms[0].get_vector()  # atom：P: 2.972	-1.431	0.126
                Atom2 = atoms[1].get_vector()  # atom：O1P: 2.452	-2.336	1.181
                Atom3 = atoms[3].get_vector()  # atom: O5': 2.153	-0.069	0.078
                Atom4 = atoms[4].get_vector()  # atom: C5': 2.38	0.987	-0.899
                Atom5 = atoms[5].get_vector()  # atom: C4': 1.31	2.044	-0.745
                Atom6 = atoms[6].get_vector()  # atom: O4': 0.021	1.499	-1.113
                Atom7 = atoms[7].get_vector()  # atom: C1': -0.942	1.76	-0.087
                Atom8 = atoms[8].get_vector()  # atom: N1: -1.877	0.628	-0.035
                Atom9 = atoms[9].get_vector()  # atom: C6: -1.387	-0.654	0.128
                avg_atoms = np.c_[Atom1[:], Atom2[:], Atom3[:], Atom4[:], Atom5[:], Atom6[:],
                                  Atom7[:], Atom8[:], Atom9[:]].mean(-1)

                vector_1 = (Atom2 - Atom1)[0:3]
                vector_2 = (Atom3 - Atom1)[0:3]
                vector_3 = (Atom4 - Atom3)[0:3]
                vector_4 = (Atom5 - Atom4)[0:3]
                vector_5 = (Atom6 - Atom5)[0:3]
                vector_6 = (Atom7 - Atom6)[0:3]
                vector_7 = (Atom8 - Atom7)[0:3]
                vector_8 = (Atom9 - Atom8)[0:3]

                v_1 = np.array([-0.52, -0.905, 1.055])  # new
                v_2 = np.array([-0.819, 1.362, -0.048])  # new
                v_3 = np.array([0.227, 1.056, -0.977])  # new
                v_4 = np.array([-1.07, 1.057, 0.154])  # new
                v_5 = np.array([-1.289, -0.545, -0.368])  # new
                v_6 = np.array([-0.963, 0.261, 1.026])  # new
                v_7 = np.array([-0.935, -1.132, 0.052])  # new
                v_8 = np.array([0.49, -1.282, 0.163])  # new

                T1 = RM.Calculate(v_1, vector_1, v_2, vector_2)
                T2 = RM.Calculate(v_3, vector_3, v_4, vector_4)
                T3 = RM.Calculate(v_4, vector_4, v_5, vector_5)
                T4 = RM.Calculate(v_5, vector_5, v_6, vector_6)
                T5 = RM.Calculate(v_7, vector_7, v_8, vector_8)

                T_Set = [T1, T2, T3, T4, T5]
                Z = Atom7[2]
                '''Matrix performing'''
                temp = MatrixPerforming(MatrixDir, aa.get_resname(), Z, T_Set, avg_atoms, mode,
                                        z_init, z_end)
                Data_DT = np.append(Data_DT, temp).reshape(-1, 3)
            '''在这个地方加几行，把上面每个if返回的data合并，然后写入文件'''
            '''我的想法是，在每个if下加一个SaveData函数，保存计算出来的data，最后再拼接起来'''

            TempSet = np.concatenate(
                (Data_TRP, Data_HIS, Data_ASP, Data_ALA, Data_GLU, Data_PHE, Data_ARG, Data_SER,
                 Data_GLY, Data_TYR, Data_VAL, Data_THR, Data_ILE, Data_PRO, Data_MET, Data_ASN,
                 Data_CYS, Data_GLN, Data_LEU, Data_LYS, Data_HEM, Data_DA, Data_DG, Data_DC,
                 Data_DT),
                axis=0)
            DataSet = np.append(DataSet, TempSet).reshape(-1, 3)
        DataSet = DataSet[np.all(DataSet != 0, axis=1)]
        SaveFiles(DataSet=DataSet, SavePath=f'{TxtSaveDir}/{SaveName}')
