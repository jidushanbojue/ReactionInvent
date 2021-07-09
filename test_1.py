from rdkit import Chem
m = Chem.MolFromSmiles('[CH3:1][S:2](=[O:3])(=[O:4])[c:5]1[cH:6][cH:7][c:8](-[c:9]2[cH:10][n:11][n:12](-[c:30]3[cH:29][cH:28][c:27]([S:26][CH3:25])[cH:32][cH:31]3)[c:13](=[O:14])[c:15]2-[c:16]2[cH:17][cH:18][c:19]([F:20])[cH:21][cH:22]2)[cH:23][cH:24]1')
nm = Chem.FragmentOnBonds(m, (11,12), addDummies=True)
Chem.MolToSmiles(nm)
mol_list = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=False)
Chem.MolToSmiles(mol_list[0])
Chem.MolToSmiles(mol_list[1])