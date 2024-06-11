#Inscope改善方案:Learning rate, batch size1024 , 檔案先shuttle好不要放在porgram裡, fold重啟，+dropout，增worker降quene,改善FP反應物字典給值做法，改用簡單的reaction rules FP(rdkit內建版)
import os
import random
import numpy as np
import tensorflow as tf
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm, trange
from collections import defaultdict
from highway_layer import Highway
#匯入深度學習的框架函式庫：keras
import keras
from keras import backend as K
from keras.initializers import Constant
from keras.utils import plot_model
#keras用以建立模型架構的函數
from keras.models import Sequential, load_model, Model

#keras中建立深度學習layer的函數

from keras.layers import Dense, Dropout, BatchNormalization, Activation, Multiply, Add, Lambda, Input
from keras import metrics, losses
#keras訓練演算法函數
from keras import regularizers
from keras.optimizers import Adam

#keras提早判停的函數
from keras.callbacks import EarlyStopping, ModelCheckpoint

#it's hard to reproduce results, so close all seeds
#os.environ['PYTHONHASHSEED'] = '0'
#np.random.seed(0)
#tf.set_random_seed(0)
#random.seed(0)

#to solve problem:Blas GEMM launch failed
from keras.backend.tensorflow_backend import set_session
config = tf.ConfigProto()
#config = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
config.gpu_options.allocator_type = 'BFC' #A "Best-fit with coalescing" algorithm, simplified from a version of dlmalloc.
config.gpu_options.per_process_gpu_memory_fraction = 0.95
config.gpu_options.allow_growth = True
set_session(tf.Session(config=config)) 

import gc
import pickle

def fps_to_arr(fps):
    """Faster conversion to ndarray"""
    arrs = []
    for fp, info in zip(fps[0],fps[1]):
        onbits = list(fp.GetOnBits())
        arr = np.zeros(fp.GetNumBits())
        for onbit in onbits:
            arr[onbit] = len(info[onbit])
        arrs.append(arr)
    arrs = np.array(arrs)
    return arrs




def fingerprint_mols(mols, fp_dim):
    fps = []
    infos = []
    for mol in mols:
        mol = Chem.MolFromSmiles(mol)
        info={}
        # Necessary for fingerprinting
        # Chem.GetSymmSSSR(mol)

        # "When comparing the ECFP/FCFP fingerprints and
        # the Morgan fingerprints generated by the RDKit,
        # remember that the 4 in ECFP4 corresponds to the
        # diameter of the atom environments considered,
        # while the Morgan fingerprints take a radius parameter.
        # So the examples above, with radius=2, are roughly
        # equivalent to ECFP4 and FCFP4."
        # <http://www.rdkit.org/docs/GettingStartedInPython.html>
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=int(fp_dim), useChirality=1, bitInfo=info)
        # fold_factor = fp.GetNumBits()//fp_dim
        # fp = DataStructs.FoldFingerprint(fp, fold_factor)
        fps.append(fp)
        infos.append(info)
    return fps, infos

def preprocess(X, fp_dim):
    # Compute fingerprints
    FPs = fps_to_arr(fingerprint_mols(X, fp_dim))
    # Apply variance threshold
    # return np.log(X[:,self.idx] + 1) 
    #FPs = np.log(dataX[:,idx]+1)
#    FPs = np.log(dataX+1)
    return FPs
def smi_list_from_str(inchis):
    '''string separated by ++ to list of RDKit molecules'''
    return [inchi.strip() for inchi in inchis.split('++')]

class DataGenerator(keras.utils.Sequence):
    
    def __init__(self, X, y, z, batch_size=1, shuffle=True, fp_dim=16384, recfp_dim=2048):
        self.batch_size = batch_size
        self.X = X
        self.y = y
        self.z = z
        self.indexes = np.arange(len(self.X))
        self.shuffle = shuffle
        self.fp_dim = fp_dim
        self.recfp_dim = recfp_dim

    def __len__(self):
        #计算每一个epoch的迭代次数
        return int(np.floor(len(self.X) / int(self.batch_size)))

    def __getitem__(self, index):
        #生成每个batch数据，这里就根据自己对数据的读取方式进行发挥了
        # 生成batch_size个索引
        batch_indexs = self.indexes[index*self.batch_size:(index+1)*self.batch_size]
        # 根据索引获取datas集合中的数据
        batch_datasX = [self.X[k] for k in batch_indexs]
        batch_datasy = [self.y[k] for k in batch_indexs]
        batch_datasz = [self.z[k] for k in batch_indexs]
        # 生成数据
        X = preprocess(batch_datasX, self.fp_dim)
        y = np.zeros((len(batch_datasy),self.recfp_dim))
        for i,a in enumerate(batch_datasy):
            n = np.zeros((1,self.recfp_dim))
            for b in smi_list_from_str(a):
                n += preprocess([b], self.recfp_dim)
            p = X[i].reshape((-1,self.recfp_dim))    
            y[i] = np.sum(p, 0, keepdims=True)- n
        z = np.array(batch_datasz)
#        y = y.astype(np.int64)
        return [X, y], [z]

    def on_epoch_end(self):
        #在每一次epoch结束是否需要进行一次随机，重新随机一下index
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

def fold(x):
    z=tf.subtract(x[0], x[1])
#    z_shape=tf.Tensor.shape(z)

#    z_shape=z.get_shape().as_list()
    zv=tf.reshape(z,[-1,8,2048])
    return tf.reduce_sum(zv, 1) 

def cosine(x):
    prod_net = x[0]
    react_net = x[1]
#    prod_norm = tf.nn.l2_normalize(prod_net, axis=-1)
#    react_norm = tf.nn.l2_normalize(react_net, axis=-1)
    cosine_sim = tf.reduce_sum(tf.multiply(prod_net, react_net), axis=-1,keepdims=True)
#    cosine_sim = tf.squeeze(cosine_sim,[1])
#    return tf.nn.sigmoid(cosine_sim)
    return tf.nn.sigmoid(cosine_sim)
# get average auc between different batches over the epoch, so don't use. otherwise validation process  always get wrong results
def auc2(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc

# AUC for a binary classifier, this AUC is a little underestimated due to minimum areas.
def auc1(y_true, y_pred):
    ptas = tf.stack([binary_PTA(y_true,y_pred,k) for k in np.linspace(0, 1, 1000)],axis=0)
    pfas = tf.stack([binary_PFA(y_true,y_pred,k) for k in np.linspace(0, 1, 1000)],axis=0)
    pfas = tf.concat([tf.ones((1,)) ,pfas],axis=0)
    binSizes = -(pfas[1:]-pfas[:-1])
    s = ptas*binSizes*1000000
    return K.sum(s, axis=0)
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# PFA, prob false alert for binary classifier(FPR)
def binary_PFA(y_true, y_pred, threshold=K.variable(value=0.5)):
    y_pred = K.cast(y_pred >= threshold, 'float32')
    # N = total number of negative labels
    N = K.sum(1 - y_true)
    # FP = total number of false alerts, alerts from the negative class labels
    FP = K.sum(y_pred - y_pred * y_true)
    return FP/N
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# P_TA prob true alerts for binary classifier(TPR)
def binary_PTA(y_true, y_pred, threshold=K.variable(value=0.5)):
    y_pred = K.cast(y_pred >= threshold, 'float32')
    # P = total number of positive labels
    P = K.sum(y_true)
    # TP = total number of correct alerts, alerts from the positive class labels
    TP = K.sum(y_pred * y_true)
    return TP/P
# PFA, prob false alert for binary classifier(FPR)
def FPR(y_true, y_pred):
    y_pred = K.cast(y_pred >= 0.5, 'float32')
    # N = total number of negative labels
    N = K.sum(1 - y_true)
    # FP = total number of false alerts, alerts from the negative class labels
    FP = K.sum(y_pred - y_pred * y_true)
    return FP/N
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# P_TA prob true alerts for binary classifier(TPR)
def TPR(y_true, y_pred):
    y_pred = K.cast(y_pred >= 0.5, 'float32')
    # P = total number of positive labels
    P = K.sum(y_true)
    # TP = total number of correct alerts, alerts from the positive class labels
    TP = K.sum(y_pred * y_true)
    return TP/P

# ACC= (TP + TN) / (P + N)
def ACCR(y_true, y_pred):
    y_pred = K.cast(y_pred >= 0.9, 'float32')
    # P = total number of positive labels
    P = K.sum(y_true)
    # N = total number of negative labels
    N = K.sum(1 - y_true)    
    # TP = total number of correct alerts, alerts from the positive class labels
    TP = K.sum(y_pred * y_true)
    # TN = total number of correct alerts, alerts from the negtive class labels
    TN = K.sum((1-y_pred) * (1-y_true))    
    return (TP+TN)/(P+N)

def aucn(y_true, y_pred):
    return -1*auc1(y_true, y_pred)

if __name__ == '__main__':
                #設定訓練參數和訓練模型存放路徑
    #batch_size = 3
    batch_size = 512
    #num_classes = 6
    #epochs = 2000
    epochs = 100
    seed=0
    #validation spilt
    #spilt=0
    spilt=0.1
    #for variance threshold
    #fp_dim=1e6
    fp_dim=16384
    recfp_dim=2048
    #recfp_dim=16384
    model_name = 'trained_model_inscope_'+str(seed)
    model_name1 = 'trained_model_inscope_'+str(seed)
    save_dir = os.path.join(os.getcwd(), 'saved_models')




    print('Loading data...')
    tem_simp = set()
    prods = []
    reacs = []
    labels = []
    '''
    with open('data/inscopedata.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            tem_simp.add(l.strip())

    with open('data/inscopedata2.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            tem_simp.add(l.strip()) 
    #print('check data:', tem_simp)

    with open('data/inscopedata4.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            tem_simp.add(l.strip()) 
        
    for item in tem_simp:
        prod,reac,label = item.split('\t')
        prods.append(prod)
        reacs.append(reac)
        labels.append(int(label))
    #print('check samples:', labels[1000000:1000010])
    print('total samples:', len(tem_simp))    
    # Shuffle
    xyz = list(zip(prods, reacs, labels))
    xyz.sort()
    random.seed(seed)
    random.shuffle(xyz)
    
    prods, reacs, labels = zip(*xyz)
    '''
    '''
    with open('data/inscopedatatest.dat', 'r') as f:
        for l in tqdm(f, desc='inscopedata'):
            tem_simp.add(l.strip()) 
        
    for item in tem_simp:
        prod,reac,label = item.split('\t')
        prods.append(prod)
        reacs.append(reac)
        labels.append(int(label))
        
    data_spilt= round(len(prods)*(1-spilt))
    x_train = prods[:data_spilt]
    x_test = prods[data_spilt:]
    y_train = reacs[:data_spilt]
    y_test = reacs[data_spilt:]
    z_train = labels[:data_spilt]
    z_test = labels[data_spilt:]
    '''
    #print('traindata:',x_train[:2],y_train[:2],z_train[:2])
    #print('testdata:',x_test[:2],y_test[:2],z_test[:2])
    
    with open('data/x_train0-smaO.pickle', 'rb') as f:
        x_train = pickle.load(f)
    with open('data/x_train0-smaO.pickle', 'rb') as f:
        x_trainM = pickle.load(f)
        
#    with open('data/x_test0.pickle', 'rb') as f:
#        x_test = pickle.load(f)
    with open('data/x_train0-smaO.pickle', 'rb') as f:
        x_test = pickle.load(f)

    with open('data/y_train0-smaO.pickle', 'rb') as f:
        y_train = pickle.load(f)
    with open('data/y_train0-smaO.pickle', 'rb') as f:
        y_trainM = pickle.load(f)
        
    with open('data/y_train0-smaO.pickle', 'rb') as f:
        y_test = pickle.load(f)
    with open('data/z_train0-smaO.pickle', 'rb') as f:
        z_train = pickle.load(f)
    with open('data/z_train0-smaMXO.pickle', 'rb') as f:
        z_trainM = pickle.load(f)        
    with open('data/z_train0-smaO.pickle', 'rb') as f:
        z_test = pickle.load(f) 
        
    with open('data/x_train0-smaO.pickle', 'rb') as f:
        x_trainM2 = pickle.load(f)
    with open('data/y_train0-smaO.pickle', 'rb') as f:
        y_trainM2 = pickle.load(f)
    with open('data/z_train0-smaO.pickle', 'rb') as f:
        z_trainM2 = pickle.load(f)  
        
    with open('data/x_train0-smaO.pickle', 'rb') as f:
        x_trainM3 = pickle.load(f)
    with open('data/y_train0-smaO.pickle', 'rb') as f:
        y_trainM3 = pickle.load(f)
    with open('data/z_train0-smaO.pickle', 'rb') as f:
        z_trainM3 = pickle.load(f)        
    
    with open('data/x_train0-smaM.pickle', 'rb') as f:
        x_trainM4 = pickle.load(f)
    with open('data/y_train0-smaM.pickle', 'rb') as f:
        y_trainM4 = pickle.load(f)
    with open('data/z_train0-smaM.pickle', 'rb') as f:
        z_trainM4 = pickle.load(f)        
    
    x_train=list(x_train)
    y_train=list(y_train)
    z_train=list(z_train)
    



    
    
    '''
    x_train1=x_train
    y_train1=y_train
    z_train1=z_train
    print('length1:',len(x_train1))
    x_train.extend(x_train1)
    y_train.extend(y_train1)
    z_train.extend(z_train1)
    print('length2:',len(x_train1))
    x_train.extend(x_train1)
    y_train.extend(y_train1)
    z_train.extend(z_train1)
    print('length3:',len(x_train1))
    x_train.extend(x_train1)
    y_train.extend(y_train1)
    z_train.extend(z_train1)
    print('length4:',len(x_train1))
    x_train.extend(x_train1)
    y_train.extend(y_train1)
    z_train.extend(z_train1)
    print('length5:',len(x_train1))
    #x_train.extend(x_train1)
    #y_train.extend(y_train1)
    #z_train.extend(z_train1)
    #print('length6:',len(x_train1))

    
    
    '''

    x_trainM=list(x_trainM)
    y_trainM=list(y_trainM)
    z_trainM=list(z_trainM)
    x_trainM2=list(x_trainM2)
    y_trainM2=list(y_trainM2)
    z_trainM2=list(z_trainM2)
    x_trainM3=list(x_trainM3)
    y_trainM3=list(y_trainM3)
    z_trainM3=list(z_trainM3)
    
    x_trainM4=list(x_trainM4)
    y_trainM4=list(y_trainM4)
    z_trainM4=list(z_trainM4)
 
    
    x_trainM4=x_trainM4[:240000]
    y_trainM4=y_trainM4[:240000]
    z_trainM4=z_trainM4[:240000]
    
    #x_trainM2=x_trainM2[12000:25200]
    #y_trainM2=y_trainM2[12000:25200]
    #z_trainM2=z_trainM2[12000:25200]
    
    #x_trainM=x_trainM[:120000]
    #y_trainM=y_trainM[:120000]
    #z_trainM=z_trainM[:120000]


#    x_trainM=x_trainM[120000:132000]
#    y_trainM=y_trainM[120000:132000]
#    z_trainM=z_trainM[120000:132000]    
    
    '''
    x_trainM.extend(x_trainM)
    y_trainM.extend(y_trainM)
    z_trainM.extend(z_trainM)
    x_trainM.extend(x_trainM)
    y_trainM.extend(y_trainM)
    z_trainM.extend(z_trainM)
    x_trainM.extend(x_trainM)
    y_trainM.extend(y_trainM)
    z_trainM.extend(z_trainM)



    
    
    '''

    '''
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
    x_trainM2.extend(x_trainM2)
    y_trainM2.extend(y_trainM2)
    z_trainM2.extend(z_trainM2)
      



    
    
    '''
    

    x_train.extend(x_trainM)
    y_train.extend(y_trainM)
    z_train.extend(z_trainM)

#    x_train.extend(x_trainM2)
#    y_train.extend(y_trainM2)
#    z_train.extend(z_trainM2) 
#    x_train.extend(x_trainM3)
#    y_train.extend(y_trainM3)
#    z_train.extend(z_trainM3)
#    x_train.extend(x_trainM4)
#    y_train.extend(y_trainM4)
#    z_train.extend(z_trainM4)
    
   
    xyz = list(zip(x_train, y_train, z_train))
#    xyz.sort()
    random.seed(seed)
    random.shuffle(xyz)
    x_train, y_train, z_train= zip(*xyz) 
    '''     
    x_train=list(x_train)
    y_train=list(y_train)
    z_train=list(z_train)
    x_test=list(x_test)
    y_test=list(y_test)
    z_test=list(z_test)
    prods=x_train
    prods.extend(x_test)
    reacs=y_train
    reacs.extend(y_test)
    labels=z_train
    labels.extend(z_test)
    xyz = list(zip(prods, reacs, labels))
    xyz.sort()
    random.seed(seed)
    random.shuffle(xyz)
    prods, reacs, labels = zip(*xyz) 
    data_spilt= round(len(prods)*(1-spilt))
    x_train = prods[:data_spilt]
    x_test = prods[data_spilt:]
    y_train = reacs[:data_spilt]
    y_test = reacs[data_spilt:]
    z_train = labels[:data_spilt]
    z_test = labels[data_spilt:]
    
    #with open('data/z_testeval.pickle', 'wb') as f:
        #pickle.dump(z_test,f) 
    #with open('data/y_testeval.pickle', 'wb') as f:
        #pickle.dump(y_test,f) 
    #with open('data/x_testeval.pickle', 'wb') as f:
        #pickle.dump(x_test,f)         

    '''
    
    print('shuffle is over...')
    print(z_train[:20])
    print(z_test[:20])
    #build model
    visible = Input(shape=(fp_dim,))
    hidden = Lambda(lambda x: tf.math.log(x+1))(visible)
    hidden = Dense(1024, activation='elu')(hidden)
    hidden = Dropout(0.3)(hidden)

    # only for expansion rule policynet
    for _ in range(5):
        hidden = Highway()(hidden)
    #    hidden = Dropout(0.4)(hidden)
    #another branch
    #visible1 = Input(shape=(fp_dim,))
    visible2 = Input(shape=(recfp_dim,))
    #hidden1 = Lambda(fold)([visible, visible2])
    hidden1 = Dense(1024, activation='elu')(visible2)

    output = Lambda(cosine)([hidden, hidden1])
    #,output_shape=(1,)
    
    model = Model(inputs=[visible,visible2], outputs=output)
    # summarize layers
    print(model.summary())
    # plot graph
    #plot_model(model, to_file='expansionpolicynet_graph.png')
    # 初始化Adam optimizer
    opt = keras.optimizers.Adam(lr=0.0001)

    # 設定訓練方式，包含loss、optimizer..)
    loss1=losses.binary_crossentropy
    #loss='binary_crossentropy'
    model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  metrics=[metrics.binary_accuracy, ACCR, auc1, TPR, FPR, loss1])
    #metrics.binary_accuracy, ACCR, auc2, auc1, TPR, FPR
    # early stop存放模型設置


    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    model_path = os.path.join(save_dir, model_name)
    model_path1 = os.path.join(save_dir, model_name1)
    checkpoint = ModelCheckpoint(model_path, monitor='val_auc1', save_best_only=1, verbose=1, mode='max',period=1)

    # early stop參數設定
    earlystop = EarlyStopping(monitor='val_auc1', patience=15, verbose=1, mode='max')

    #continue training

    #K.clear_session()
    #gc.collect()
    #del model  # 删掉存在的模型

    #返回一个编译好的模型
    #与删掉的模型相同
    #model = load_model(model_path1, custom_objects={'ACCR': ACCR,'auc2': auc2,'auc1': auc1,'TPR': TPR, 'FPR': FPR,'Highway': Highway,'fold': fold,'cosine': cosine, 'tf': tf, 'loss1': loss1})
    #model.compile(loss='binary_crossentropy',
    #           optimizer=opt,
    #            metrics=[metrics.binary_accuracy, ACCR, auc1, loss1])
    #metrics.binary_accuracy, ACCR, auc2, auc1, TPR, FPR
    # 開始訓練
    #x_train1 = x_train[:500000]
    #y_train1 = y_train[:500000]
    #z_train1 = z_train[:500000]
    training_generator = DataGenerator(X=x_train, y=y_train, z=z_train, batch_size=batch_size, shuffle=True, fp_dim=fp_dim, recfp_dim=recfp_dim)
    validation_gen = DataGenerator(X=x_test, y=y_test, z=z_test, batch_size=batch_size, shuffle=True, fp_dim=fp_dim, recfp_dim=recfp_dim)    

    model_history = model.fit_generator( 
                    generator=training_generator,
                    epochs=70,
#                    class_weight = {1:1., 0:5.},
                    validation_data=validation_gen,
                    verbose=1,
                    initial_epoch=0,
                    workers=8, 
                    use_multiprocessing=1, 
#                    shuffle=False,
#                    max_queue_size = 3, 
#                    callbacks=[checkpoint]
                    callbacks=[earlystop, checkpoint]
                    )