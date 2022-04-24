import os
import urllib3
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import Draw



cnames = ['C00138', 'C00139']
savedir = "./KEGG_Compound_MOLiles/"

def get_molfile(cname, savepath='./'):
    url = "https://www.kegg.jp/entry/-f+m+" + cname
    if not os.path.exists(savepath):os.mkdir(savepath)
    http = urllib3.PoolManager()
    head = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.109 Safari/537.36"
    }
    mol = http.request('GET', url, headers=head)
    moldata = mol.data.decode('UTF-8')
    with open(savepath + cname + '.mol', 'w', encoding='utf-8') as f:
        f.write(moldata)


mols = []
for cname in cnames:
    get_molfile(cname, savedir)
    mol = Chem.MolFromMolFile(savedir + cname + '.mol')
    mols.append(mol)
    print(Chem.MolToSmiles(mol))

fps = [Chem.RDKFingerprint(x) for x in mols]
sm01 = DataStructs.FingerprintSimilarity(fps[0], fps[1])
sm02 = DataStructs.FingerprintSimilarity(fps[0], fps[2])
sm12 = DataStructs.FingerprintSimilarity(fps[1], fps[2])

print("similarity between mol1 and mol2: %.2f" %
      sm01)  # similarity between mol1 and mol2: 0.93
print("similarity between mol1 and mol3: %.2f" %
      sm02)  # similarity between mol1 and mol3: 0.87
print("similarity between mol2 and mol3: %.2f" %
      sm12)  # similarity between mol2 and mol3: 0.93