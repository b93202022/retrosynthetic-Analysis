import re
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import Chem, RDLogger
from itertools import chain, permutations
from multiprocessing import Pool, freeze_support
from collections import defaultdict
import random
from functools import partial

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

def clean(lines):
    a,b,c =lines.strip().split('\t')
    return [a,b,c]
def seedgen(seed):
    expansion_rules = []
    with open('data/expansion_expansion.dat', 'r') as f:
        for i, l in tqdm(enumerate(f), desc='expansion'):
            rule = l.strip()
            expansion_rules.append(rule)
    random.seed(seed)
    while True:
#        random.seed(seed)
##        random.shuffle(expansion_rules)
        rules=random.sample(expansion_rules, 500)
#        seed +=1
##        yield expansion_rules[:500]
        yield rules
def read(combo):
    rule = combo[0][0]
    prod = combo[0][1]
    reac = combo[0][2]
    rules = combo[1]
#    rules =combo[1]
#    expansion_rules = []
#    with open('data/expansion_expansion.dat', 'r') as f:
#        for i, l in tqdm(enumerate(f), desc='expansion'):
#            rule = l.strip()
#            expansion_rules.append(rule)
 #   rules = expansion_rules       
    prod_to_reacs = defaultdict(set)
    prod_to_noreacs = defaultdict(set)
# rule, prod, reac are strings not lists
#    rule, prod, reac = line.strip().split('\t')
#    prod_to_reacs[prod].add(reac)
    #   print(rule)
#    seed = rule+prod+reac
#    random.seed(seed)
#    random.shuffle(rules)
    rulesrad=rules
    for r in rulesrad:
        if r== rule: continue
#        print(r)
        retro_canonical = convert_to_retro(r)
        #   print(retro_canonical)
        rxn = AllChem.ReactionFromSmarts(retro_canonical)
        rcts= mol_list_from_str(reac)
        if len(rcts) != rxn.GetNumReactantTemplates(): continue
        try:
            outcomes = rxn.RunReactants(rcts)
            if not outcomes: continue
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
    expansion_rules = []
    tem_simp = set()
    seed =0
    random.seed(seed)
    
    with open('data/expansion_expansion.dat', 'r') as f:
        for i, l in tqdm(enumerate(f), desc='expansion'):
            rule = l.strip()
            expansion_rules.append(rule)
#    random.shuffle(expansion_rules)
#    with open('data/templates_expansion1.dat', 'r') as f:
#        for l in tqdm(f, desc='products'):
#                tem_simp.add(l.strip())
#    with open('data/templates_expansion3.dat', 'w') as f:
#        f.write('\n'.join(tem_simp))                           
    combo=[]
    with open('data/templates_expansion3.dat', 'r') as f:
#        for i, l in tqdm(enumerate(f), desc='expansion'):
#            combo.append(clean(l))

    #    for l in tqdm(f, desc='products'):
#        seeds=[a for a in range(len(f.readlines()))]
#        seeds=[expansion_rules[:55608] for a in range(len(f.readlines()))]
#        print('seeds length',len(seeds))
##        seeds= [expansion_rules[:55608] for a in range(len(combo))]
#        combo=[]
#            combo.append(clean(f))
#        combo[0]=f
#        combo[1]=seeds
##        comboo=zip(combo,seeds)
#        comboo=zip(combo,seedgen(expansion_rules))
#        comboo=zip(map(clean, f),seedgen(expansion_rules))
#        comboo=zip(map(clean, f),seeds)
#        comboo.append(combo)
#        comboo.append(seeds)
        with Pool() as p:
#            for reacs,noreacs in tqdm(p.imap(partial(read, rules=expansion_rules), zip(map(clean, f),seeds))):
#          for reacs,noreacs in tqdm(p.imap(read, comboo)):  
            for reacs,noreacs in tqdm(p.imap(read,zip(map(clean, f),seedgen(seed)))):
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
    with open('data/inscopedata3.dat', 'w') as f:
        f.write('\n'.join(['\t'.join(rxn_prod) for rxn_prod in transforms]))

    print('total positive Expansion rules examples:',Tolpos) 
    print('total negative Expansion rules examples:',Tolneg)