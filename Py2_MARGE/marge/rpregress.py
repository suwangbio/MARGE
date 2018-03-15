import argparse,math,os,sys,tables
import numpy as np
from operator import itemgetter, attrgetter, methodcaller
from sklearn import linear_model, decomposition, datasets
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn import cross_validation
from sklearn import metrics
import regions
import random


UNCHANGED,DOWN,UP,TARGET = 0,1,2,3
statdict = { UNCHANGED:'.', DOWN:'DOWN', UP:'UP', TARGET:'TARGET' }

class dataset(object):
    def __init__(self):
        self.index = None
        self.info  = None # include chrom and TSS and gene symbol
        self.x     = None
        self.y     = None
        self.bwdir = None 
        self.bwfiles = None
        self.rpfiles = None 

def gene_sym(symfile):
    """
    One representative of each gene symbol is chosen. 
    """

    fp = open(symfile)
    symdict = {}
    for line in fp:
        f = line.strip().split('\t')
        #print f
        g = f[3]
        IDs = g.split(':')
        #print IDs
        if IDs[1] not in symdict:
            symdict[IDs[1]] = [IDs[0], f[0]]
 
    rsymdict = {}
    for elem in symdict:
        rsymdict[symdict[elem][0]] = [elem,symdict[elem][1]]
    fp.close()

    return rsymdict

def scalartransform(x):
    pcount = 1
    x = np.log2(x+pcount)
    return x

def transform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount 
    med = np.median( xt,0 )
    x = np.log2(xt) - np.log2(med)
    return x

def sqrttansform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount
    med = np.median( xt,0 )
    x = np.sqrt(xt) - np.sqrt(med)
    return x
#====================================================================

def read_hdf5( file_name, sample_names ):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 

    h5file = tables.open_file( file_name, driver="H5FD_CORE")

    X = None

    for elem in sample_names:
        a = h5file.get_node("/", elem )
        m = a.read()
        if X is None:
            X = m
        else:
            X = numpy.vstack((X,m))

    X = X.transpose()

    h5file.close()
    return X

def getSampleNames_hdf5(rp_hdf5):
    h5file = tables.open_file( rp_hdf5, driver="H5FD_CORE")
    samplenames = []
    for array in h5file.walk_nodes("/","Array"):
        if array.name not in samplenames:
            samplenames.append(array.name)
        else:
            continue
    return samplenames

def readregpotfiles(sym,genome,symfile,hdf5s,samplenames):

    # make list of regulatory potential file
    
    # read in regulatory potential from files in directory 

    index = None
    x = None
    nsamples = len(samplenames)


    for k,name in enumerate(samplenames):
       
        bed = regions.interval(genome)
        bed.read_bed(symfile,scorefield=4)
        if index == None:
            print name
            index = {}
            info = {}
            i = 0
            for j,geneids in enumerate(bed.name):
                #print geneid
                geneids = geneids.split(':')
                geneid = geneids[0]
                #print geneid
                if geneid in sym:
                    symid = sym[geneid][0]
                    #print symid
                    if symid not in index:
                        index[symid] = i
                        info[symid] = [bed.chrom[j],bed.start[j]]
                        i += 1
            ngenes = len(index)
            x = np.zeros((ngenes,nsamples))
            print np.shape(x)

        RP = None
        for hdf5 in hdf5s:
            if name in getSampleNames_hdf5(hdf5):
                RP = read_hdf5( hdf5, [name] )
            else:
                continue
        if RP is None:
            print "Please input hdf5 files with ',' as the separator"

        for i,geneids in enumerate(bed.name):
            geneids = geneids.split(':')
            geneid = geneids[0]
            if geneid in sym:
                symid = sym[geneid][0]
                rp = RP[i]
                try:
                    x[index[symid],k] = rp
                    #print x[index[symid],k]
                except:
                    pass
    #tempRP = x[:,1:10]
    #outf = open('tempRP.txt','w')
    #for t in range(np.shape(tempRP)[0]):
    #    for tt in tempRP[t,:]:
    #        outf.write(str(tt) + '\t')
    #    outf.write('\n')
    #outf.close()
    z         = dataset()  
    z.rpfiles = samplenames 
    z.x = x
    z.index = index
    
    z.info  = info
    return z

def read_genelistResponse(sym, fname, index, gname2):
    ID,RESPONSE = range(2)
    status = np.zeros( len(index) )
    train_index = []
    test_index = []

    fp = open(fname)
    
    for line in fp:
        f = line.split()
        if gname2:
            try:
                i = index[f[ID]]
            except:
                continue
        else:
            try:
                i = index[sym[f[ID]][0]]
            except:
                continue

        if f[RESPONSE] == '0':
            status[i] = UNCHANGED
        if f[RESPONSE] == '1':
            status[i] = TARGET
    print 'file: %s\ttarget: %d\tunchanged: %d\n' % ( fname, sum( status == TARGET ), sum( status == UNCHANGED ) )
    return (status,train_index,test_index) 

def read_genelistOnly(sym, fname, index, gname2, sepby_chrom=True):
    
    status = np.zeros( len(index) )
    print index.keys()[0:20]

    train_chroms = ['chr1','chr3','chr5','chr7','chr9','chr11','chr13','chr15','chr17','chr19','chr21']
    test_chroms = ['chr2','chr4','chr6','chr8','chr10','chr12','chr14','chr16','chr18','chr20','chr22']
    train_index = []
    test_index = []

    fp = open(fname).readlines()
    genes = [g.strip() for g in fp]

    allgenes = sym.keys()
    print allgenes[0:20]

    for ag in allgenes:
        if gname2:
            try:
                i = index[sym[ag][0]]
                if sym[ag][1] in train_chroms:
                    train_index.append(i)
                elif sym[ag][1] in test_chroms:
                    test_index.append(i)
                else:
                    pass
                #print i
                if sym[ag][0] in genes:
                    #print sym[ag][0]
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
            except:
                continue
        else:
            try:
                i = index[sym[ag][0]]
                if sym[ag][1] in train_chroms:
                    train_index.append(i)
                elif sym[ag][1] in test_chroms:
                    test_index.append(i)
                else:
                    pass
                if ag in genes:
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
            except:
                continue

    print 'file: %s\ttarget: %d\tunchanged: %d\n' % ( fname, sum( status == TARGET ), sum( status == UNCHANGED ) )
    return (status,train_index,test_index)  

def dataset_annotation(annotationf):
    #get the cell annotation for each datasetID
    inf = open(annotationf,'rU')
    ann = {}
    for line in inf:
        if line.startswith('datasetID'):
            pass
        else:
            line = line.strip().split('\t')
            ID = line[0]
            info = [line[4],line[5],line[7]]
            try:
                ann[ID] = info
            except:
                ann[ID] = 'NA'
    return ann


def regress(x, y, train_index, test_index, samplenames, name, change, maxsamples,sepby_chrom, ann):

    print x.shape
    print y.shape

    cross_validate_flag = True
    record   = []
    cvrecord = []

    col_list   = samplenames
    maxsamples = min( len(col_list), maxsamples )
    print maxsamples

    # forward step-wise regression - include the samples that produce the best cross validation results
    for i in range(maxsamples):
        #print 'first i = %s'%str(i)
        if cross_validate_flag:
            cvscore = []

            for i,elem in enumerate(col_list):
                if i not in record:
                    #print 'second i = %s'%str(i)
                    trial = record + [i]
                    rand_idx = range(x.shape[0])
                    if not sepby_chrom:
                        # select random subset for trial and testing
                        random.shuffle( rand_idx )
                        #xt = x[:,trial]
                        xt = x[:,trial]
                        xt = xt[rand_idx,:]
                        yt = y[rand_idx]
                        xt.reshape( (x.shape[0],len(trial)) )
                        X_train, X_test, y_train, y_test = cross_validation.train_test_split( xt, yt )
                        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
                        LR_l1.fit(X_train,y_train)
                        cvscore += [ np.mean( cross_validation.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ) ) ]
                    else:
                        xt = x[:,trial]
                        yt = y
                        X_train = xt[train_index,:]
                        y_train = yt[train_index]
                        
                        X_train.reshape((X_train.shape[0],len(trial)))
                        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
                        LR_l1.fit(X_train,y_train)
                        cvscore += [ np.mean( cross_validation.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ) ) ]
                else:
                    cvscore += [0]
                    continue
            k = np.argmax(cvscore) 
            cvrecord += [max(cvscore)]
            record += [k]
            #print record

    if sepby_chrom:
        xt = x[:,record]
        X_test = xt[test_index,:] 
        y_test = y[test_index]
        LR_l1.fit(X_test,y_test)
        y_test_hat = LR_l1.predict_log_proba(X_test)
        yhat = LR_l1.predict_log_proba(xt)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_test_hat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        print "=============="
    else:
        xt = x[:,record]
        xt.reshape( (xt.shape[0],len(record)) )
        LR_l1.fit(xt,y)
        yhat = LR_l1.predict_log_proba(xt)
        fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        print "=============="
    outf = open(name + '_' + change + '_regressionInfo.txt','w')
    for k,elem in enumerate(record):
        dataID = col_list[elem].split('_')[0]
        if dataID in ann.keys():
            annInfo = ann[dataID]
        else:
            annInfo = ['NA','NA','NA']

        print dataID, cvrecord[k], LR_l1.coef_[0,k], annInfo
        outf.write(dataID + '\t' + str(cvrecord[k]) + '\t' + str(LR_l1.coef_[0,k]) + '\t' + '\t'.join(annInfo) + '\n')
    outf.write('AUC = %s'%(str(round(auc,3))))
    outf.close()   

    return LR_l1.coef_, [ col_list[i]  for i in record ], yhat, record


def main( genome, rp_hdf5, gxfile, symfile, name, change, maxsamples, logtrans, sqrttrans, exptype, gname2, sepby_chrom, annotation):
    
    sym = gene_sym(symfile)
    hdf5s = rp_hdf5.strip().split(',')
    samplenames = []
    for hdf5 in hdf5s:
        sns = getSampleNames_hdf5(hdf5)
        samplenames += sns

    z   = readregpotfiles(sym,genome,symfile,hdf5s,samplenames)
  
    if logtrans:
        #z.x = scalartransform(z.x)
        z.x = transform(z.x)
    if sqrttrans:
        z.x = sqrttansform(z.x)

    rsymid = {}
    for elem in sym:
        rsymid[ sym[elem][0] ] = [elem,sym[elem][1]]  # reverse lookup 

    m = np.median(z.x,1)
    if exptype == 'LIMMA':
        (z.y,train_index,test_index) = read_limma_file(sym, gxfile, z.index, minfc=minfc, max_padj=adjp)
    elif exptype == 'DESeq':
        (z.y,train_index,test_index) = read_deseq_file(sym, gxfile, z.index, minfc=minfc, max_padj=adjp)
    elif exptype == 'Gene_Response':
        (z.y,train_index,test_index) = read_genelistResponse(sym, gxfile, z.index, gname2)
    else:
        (z.y,train_index,test_index) = read_genelistOnly(sym, gxfile, z.index, gname2, sepby_chrom)

    if change == 'down':
        print 'Do regrssion with DOWN genes'
        y = 1*( z.y == DOWN )
    elif change == 'up':
        print 'Do regrssion with UP genes'
        y = 1*( z.y == UP )
    elif change == 'target':
        print 'Do regrssion with TARGET genes'
        y = 1*( z.y == TARGET )
    else:
        print "Please choose the specfic direction, UP or DOWN or TARGET."
        sys.exit()

    x = z.x
    ann = dataset_annotation( annotation )
    coef,rprecord,yhat,record = regress( x, y, train_index, test_index, z.rpfiles, name, change, maxsamples, sepby_chrom, ann )

if __name__ == "__main__":

    try:
        parser = argparse.ArgumentParser(description="""Regression of regulatory potential to gene expression changes.""")
        parser.add_argument( '-e','--expr', dest='expr', required = True, type = str, help = 'The related differential expression file')
        parser.add_argument( '--exptype', dest='exptype',required = True, choices=['Gene_Response','Gene_Only'], type = str, help = 'The type of the expression file,\
                                             Gene_Response includes 2 columns, one is the geneID, and the other\
                                             is 1/0, 1 for target and 0 for un-target; Gene_Only includes 1 column, only the gene list of the targets. Only official gene symbol or \
                                             refseqID are allowd for the geneID.')
        parser.add_argument( '--gname2', dest='gname2', default=False, action='store_true', help = 'If this switch is on, gene or transcript IDs in files given through -e will be considered as official gene symbols, otherwise, it is RefSeqID, DEFAULT=FALSE')
        parser.add_argument( '-r','--historicalRP', dest='histRP', required = True, type = str, help = 'The file with hdf5 format which contain the H3K27ac RP information')
        parser.add_argument( '-n','--name', dest='name',required = True, type = str, help = 'The prefix of the output names')
        parser.add_argument( '-g','--genome', dest="genome", type=str, default='hg38', choices=['mm9','hg19','hg38','mm10'], required=False, help='genome')
        parser.add_argument( '--maxsamples', dest='maxsamples',  type=int, default=10, required=False, help='Maximum number of samples to include in regression model.' )
        parser.add_argument( '-l', '--logtrans', dest='logtrans',  default=False,  action='store_true', help='Use log transformation on regulatory potential (default is do nothing), this is complementary with -s, --sqrttrans.' )
        parser.add_argument( '-s', '--sqrttrans', dest='sqrttrans',  default=False,  action='store_true', help='Use sqrt transformation on regulatory potential (default is do nothing). this is complementary with -l, --logtrans.' )
        parser.add_argument( '-m', dest='sym', type=str, required=True, help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>')
        parser.add_argument( '--sc', dest='sepby_chrom',  default=False,  action='store_true', help='Whether to seperate chroms to do cross validation, default = False' )
        parser.add_argument( '-a', dest='annotation', required=True,  type=str, help='The annotation file for each dataset' )

        args = parser.parse_args()

 
        main( args.genome, args.histRP, args.expr, args.sym, args.name, 'target', args.maxsamples, args.logtrans, args.sqrttrans, args.exptype, args.gname2, args.sepby_chrom, args.annotation )

    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)

