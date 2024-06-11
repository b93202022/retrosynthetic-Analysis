import re
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import Chem, RDLogger
from itertools import chain, permutations
from multiprocessing import Pool, freeze_support
from collections import defaultdict

def mol_list_to_str(mols):
    '''List of RDKit molecules to string separated by ++'''
    inchis = [Chem.MolToSmiles(mol, allHsExplicit=True,allBondsExplicit=True) for mol in mols]
    return ' ++ '.join(inchis)

def mol_list_from_str(inchis):
    '''string separated by ++ to list of RDKit molecules'''
    return [Chem.MolFromSmiles(inchi.strip()) for inchi in inchis.split('++')]

def convert_to_retro(transform):
    '''This function takes a forward synthesis and converts it to a
    retrosynthesis. Only transforms with a single product are kept, since
    retrosyntheses should have a single reactant (and split it up accordingly).'''    

    # Split up original transform
    reactants = transform.split('>>')[0]
    products  = transform.split('>>')[1]

    # Don't force products to be from different molecules (?)
    # -> any reaction template can be intramolecular (might remove later)
    #products = products[1:-1].replace(').(', '.')

    # Don't force the "products" of a retrosynthesis to be two different molecules!
    #reactants = reactants[1:-1].replace(').(', '.')

    return '>>'.join([products, reactants])

def read(line):
    prod_to_reacs = defaultdict(set)
    prod_to_noreacs = defaultdict(set)
# rule, prod, reac are strings not lists
    rule, prod, reac = line.strip().split('\t')
#    prod_to_reacs[prod].add(reac)
    #   print(rule)
    retro_canonical = convert_to_retro(rule)
    #   print(retro_canonical)
    rxn = AllChem.ReactionFromSmarts(retro_canonical)
    try:
        outcomes = rxn.RunReactants(mol_list_from_str(reac))
        for outcome in outcomes:
            for product in outcome:
                    
                try:
                    Chem.SanitizeMol(product)
                    product.UpdatePropertyCache()
                    #create product or reactant using molfromsmarts+sanitizemol is sometimes better than molfromsmiles, but still using molfromsmiles as possible as you can
                    product=Chem.MolFromSmiles(Chem.MolToSmiles(product,allHsExplicit=True,allBondsExplicit=True))
                except Exception as e:
    #                   print('warning1: {}'.format(e))
                    #use pass is not good behavior, however i have validation finally
                    continue
                if not product:
                    continue
                prodsmi=Chem.MolToSmiles(product,allHsExplicit=True,allBondsExplicit=True)
                if  prodsmi != prod:
                    prod_to_noreacs[prodsmi].add(reac)
                    continue
                if  prodsmi == prod:
                    prod_to_reacs[prodsmi].add(reac)
                #tolri doesnt work well due to repeated same prodsmi but different product   
    #               Tolri+=1
                 
    except Exception as e:
        print('error: {}'.format(e))
        print('rxn: {}'.format(reac))
    
    return prod_to_reacs,prod_to_noreacs

if __name__ == '__main__':
    print('Loading data...')
    prod_to_reacs = defaultdict(set)
    prod_to_noreacs = defaultdict(set)
    Tolri=0
    Tolpos=0
    Tolneg=0

    with open('data/templates_expansion1.dat', 'r') as f:
    #    for l in tqdm(f, desc='products'):
        with Pool() as p:
            for reacs,noreacs in tqdm(p.imap(read, f)):
                for reac, values in reacs.items():
                    for value in values:
                        prod_to_reacs[reac].add(value)
                
                for reac, values in noreacs.items():
                    for value in values:
                        prod_to_noreacs[reac].add(value)
                    
    transforms=[]  
    for prod, reacs in prod_to_reacs.items():
        for reac in reacs:
            transforms.append((prod, reac, '1'))
            Tolpos+=1

    for prod, noreacs in prod_to_noreacs.items():
        for noreac in noreacs:
            transforms.append((prod, noreac, '0'))             
            Tolneg+=1
    with open('data/inscopedata.dat', 'w') as f:
        f.write('\n'.join(['\t'.join(rxn_prod) for rxn_prod in transforms]))

    print('total positive Expansion rules examples:',Tolpos) 
    print('total negative Expansion rules examples:',Tolneg)