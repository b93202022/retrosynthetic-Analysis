    from tqdm import tqdm, trange
    import random
    import pickle
    
    seed = 0
    #spilt = 0.1
    spilt = 0
    
    print('Loading data...')
    tem_simp = set()
    tem_val_pos = []
    tem_val_neg = []
    tem_val = []
    
    
    prods = []
    reacs = []
    labels = []
    prodsv = []
    reacsv = []
    labelsv = []    
    countp =-1
    countn =-1
    with open('data/inscopedata.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            
            tem_simp.add(l.strip())
            if l.strip().split('\t')[2] =='1':
                tem_val_pos.append(l.strip())
            if l.strip().split('\t')[2] =='0':
                tem_val_neg.append(l.strip()) 
    for l in tem_val_pos:
        countp +=1
        if countp%10 ==0:
            tem_val.append(l)
    for l in tem_val_neg:
        countn +=1
        if countn%2 ==1:
            tem_val.append(l)    

    with open('data/inscopedata2_noise_2-7.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            tem_simp.add(l.strip()) 
    #print('check data:', tem_simp)
    
    #with open('data/inscopedata4all.dat', 'r') as f:
    #    for l in tqdm(f, desc='inscopedata'):
    #        tem_simp.add(l.strip()) 
    
    #with open('data/inscopedata4.dat', 'r') as f:
    #    for l in tqdm(f, desc='inscopedata'):
    #        tem_simp.add(l.strip()) 
    print('total pos samples:', len(tem_val_pos))
    print('total neg samples:', len(tem_val_neg))
    print('total samples:', len(tem_simp))
    print('total val samples:', len(tem_val))
    tem_simp.difference_update(tem_val)
    print('total train samples:', len(tem_simp))
    
    for item in tem_simp:
        prod,reac,label = item.split('\t')
        prods.append(prod)
        reacs.append(reac)
        #labels.append(label)
        labels.append(float(label))
    #print('check samples:', labels[1000000:1000010])
    print('total train samples:', len(prods))    
    # Shuffle
    xyz = list(zip(prods, reacs, labels))
    xyz.sort()
    random.seed(seed+1)
    random.shuffle(xyz)
   
    
    #with open('data/inscopedatatest.dat', 'w') as f:
    #    f.write('\n'.join(['\t'.join(rxn_prod) for rxn_prod in xyz]))
    
    
    #'''
    prods, reacs, labels = zip(*xyz)
    
    for item in tem_val:
        prod,reac,label = item.split('\t')
        prodsv.append(prod)
        reacsv.append(reac)
        #labels.append(label)
        labelsv.append(float(label))
    #print('check samples:', labels[1000000:1000010])
    print('total val samples:', len(prodsv))    
    # Shuffle
    xyz = list(zip(prodsv, reacsv, labelsv))
    xyz.sort()
    random.seed(seed)
    random.shuffle(xyz)
   
    
    #with open('data/inscopedatatest.dat', 'w') as f:
    #    f.write('\n'.join(['\t'.join(rxn_prod) for rxn_prod in xyz]))
    
    
    #'''
    prodsv, reacsv, labelsv = zip(*xyz)
    
    #data_spilt= round(len(prods)*(1-spilt))
    x_train = prods
    x_test = prodsv
    y_train = reacs
    y_test = reacsv
    z_train = labels
    z_test = labelsv
    
    with open('data/x_train1-n2-7.pickle', 'wb') as f:
        pickle.dump(x_train, f)
    with open('data/x_test1-n2-7.pickle', 'wb') as f:
        pickle.dump(x_test, f)
    with open('data/y_train1-n2-7.pickle', 'wb') as f:
        pickle.dump(y_train, f)
    with open('data/y_test1-n2-7.pickle', 'wb') as f:
        pickle.dump(y_test, f) 
    with open('data/z_train1-n2-7.pickle', 'wb') as f:
        pickle.dump(z_train, f)
    with open('data/z_test1-n2-7.pickle', 'wb') as f:
        pickle.dump(z_test, f)   
    #'''    
    print('total val samples:', len(xyz))