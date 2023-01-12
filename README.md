# SPARC
A separable PDs-based algorithm with rotational coarse-grained model (SPARC).

######  Background ######
 This algorithm is used for rapidly calculating SERS spectra of biomolecules and determining the 3D structure of protein, which has not been available yet.

###### Install ###### 
 This code is developed based on Python, which is quite simple to 'install'. Basically, the user need to change the path of PdbSet, TxtSet, and XlsxDir in 'main.py' and put the 'positioner.py' in the same folder. It can be performed on win and Linux without any requirment of non standard hardware.

e.g. In the current folder we provided, we put the folder named 'Code' in the Desktop. Thus, we need to modify the corresponding path to:
    pdbFiles = glob.glob(r"C:\Users\Hao\Desktop\Code\PdbSet\*.pdb")
    XlsxDir = r'C:\Users\Hao\Desktop\Code\XlsxDir\\'
    TxtDir = r'C:\Users\Hao\Desktop\Code\TxtSet' 

You can then install dependency package according to the instruction in the code. Generally, Numpy, Pandas, Math, xlrd(version < 2), and Biopython are needed, which would cost several minites to install. After installing these package, you can use the algorthm by running 'main.py'.

###### Usage ######
 To note that, this algoithm only produces SERS intensity of each vibrations, which means you need draw the spectra by another software (Multiwfn is used in this work).

i.e. We saved the empirical polarizability derivatives of 20 amino acids and AGCT in XlsxDir folder. Positioner.py will read the coordinates of PDB in PdbSet, calculate the rotation matrix, invoke the PDs needed from XlsxDir, and calculate the SERS intensity to TxtSet. This process is same as illustrated in Fig. S1 in the Supporting information.

We provided two examples here:
First one is the peptide used in the manuscript (see results in Figure 4D-c, all parameters are based on theoretical equation). And you can find a output with name of 'ep18.out' in TxtSet. Second one is the conformation2 of KSI as shown in Fig. 5D(2). It is quite important to note that the EM field distribution is totally different in different system, which you need to change the expression of EM field in line 148 of 'Positioner.py'. In the current version, we provided several field distribution equation corresponding to different system, which can be used to reproduce the spectra in the manuscript. Both of two examples will cost second to produce output files.

###### Intruction of using Multiwfn ######
To note that, the usage of Multiwfn is to do the Lorentzian expansion, which means it is not the only way to obtain spectra. If you also use this software to draw the spectra, please cite Tian Lu, Feiwu Chen, J. Comput. Chem., 33, 580-592 (2012). And it can be downloaded at http://sobereva.com/multiwfn/.
When the user obtained the file in TxtSet, they can directly drag the file into 'GUI' of Multiwfn:
"Enter" ——> Choose "11" ——> Choose "2" (Then you can modified the spectra by yourself, including range, width...) ——> Choose "2"
By this way, you can obtain the spectra in the folder of Multiwfn (You can also see the spectra in advance by choose "0" at the last step).

###### Maintainer and Contributors ######

Dr. Hao Ma, Xinyu Lu, Prof. Bin Ren

For more information or instruction, please contact us: oaham@xmu.edu.cn (H.M); bren@xmu.edu.cn (B.R).

###### License ######
To use SPARC, you are required to read and agree the following terms:
(a) Currently version of SPARC is free of charge and open-source for both academic and commerical usages, anyone is allowed to freely distribute the original or their modified codes to others.
(b) SPARC can be distributed as a free component of commercial code. Selling modified version of SPARC may also be granted, however, obtaining prior consent from the original author of SPARC (Dr. Hao Ma, Xinyu Lu, Prof. Bin Ren) is needed.

Thank you very much for using SPARC.
