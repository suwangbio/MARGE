import argparse,math,numpy,os,sys,tables
import time
from sklearn import linear_model, decomposition, datasets, cluster, cross_validation, metrics


# read in list
def read_list(fname):
    fp = open(fname)
    l = []
    for line in fp.readlines():
        if not ( line[0] == '#' ):
            f = line.strip().split()
            if not f:
                continue
            else:
                sym = f[0]
                l  += [sym]
    fp.close()
    return l


def gene_sym(symfile):
    """
    One representative of each gene symbol is chosen. 
    """
    fp = open(symfile)
    symdict = {}
    refrp = []
    for line in fp.readlines():
        f = line.split('\t')
        g = f[3]
        IDs = g.split(':')
        refrp.append(IDs[0])
        if IDs[1] not in symdict:
            symdict[IDs[1]] = IDs[0]
 
    rsymdict = {}
    for elem in symdict:
        rsymdict[symdict[elem]] = elem
    fp.close()
    return (rsymdict,refrp)

# read in H3K27ac sample IDs
def read_sample_list(fname):
    fp = open(fname)
    sl = []
    for line in fp:
        line  = line.strip().split('\t')
        if not line[0].startswith('AUC'):
            sl += [line[0]]
    fp.close()
    return sl

# read in H3K27ac signal on union DHS sites for sample
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
    return X

# rank numpy vector 
def rank(a):
    """ Rank from low to high """
    temp = a.argsort()
    rank = numpy.empty(len(a), int)
    rank[temp] = numpy.arange(len(a))
    return rank

def median_norm(X):
    col_median = numpy.median(X,0)
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - col_median[i]
    return X

def quantile_norm(X):
    """ input: numpy array data_points x features """
    Y = numpy.copy(X)
    Y.sort(0)
    row_mu = Y.mean(1)
    for i in range(Y.shape[1]):
        r = rank(X[:,i])
        X[r,i] = row_mu
    return X

def rowcenter(X):
    X = (X.transpose()-X.mean(1)).transpose()
    return X

def read_table(tffile,col=1,ctype=int):
    fp = open(tffile)
    l = []
    for line in fp.readlines():
        f = line.split()
        l += [ctype(f[col])]
    return l

def getSampleNames_hdf5(rp_hdf5):
    h5file = tables.open_file( rp_hdf5, driver="H5FD_CORE")
    samplenames = []
    for array in h5file.walk_nodes("/","Array"):
        if array.name not in samplenames:
            samplenames.append(array.name)
        else:
            continue
    return samplenames

def main(samplefile,genefile,casename,genome,RP_HDF5s,UDHS_H3K27ac_HDF5s,REF_SYM,tffile=None):

    # TODO allow choice of genomes
    #refrp = REFGENE_INDEX
    
    NCLUSTERS     =  7

    fpo    = open(casename+'_enhancer_prediction.txt','w')
    log_fp = open(casename+'_enhancer_prediction.log','w')
    log2_fp = open(casename+'_centroids.txt','w')

    # read in refseq to sym dictionary
    (refdict,refrp) = gene_sym(REF_SYM)
    symdict = {}

    # read tf file
    if tffile:
        tf_in_UDHS = read_table(tffile,col=1,ctype=int)
    else:
        tf_in_UDHS = None
    # generate sym to refseq dictionary
    for elem in refdict:
        symdict[ refdict[elem] ] = elem 

    sample_names = read_sample_list(samplefile)
    RP_sample_names = ['%s_all_RP' % elem for elem in sample_names]

    RP_HDF5s = RP_HDF5s.strip().split(',')

    RP = None
    for rp_samplename in RP_sample_names:
        rp = None
        for rp_hdf5 in RP_HDF5s:
            if rp_samplename in getSampleNames_hdf5(rp_hdf5):
                rp = read_hdf5( rp_hdf5, [rp_samplename] )
            else:
                continue
        if rp is None:
            print "Please input all the hdf5 file used in regression step, with ',' as the separator"
            sys.exit()
        else:
            if RP is None:
                RP = rp
            else:
                RP =  numpy.vstack((RP,rp))
    RP = RP.transpose()
    
    # read in list of refseqs for indexing RPs
    idx_rp  = {}
    idx_all = []

    i = 0
    for j,elem in enumerate(refrp):
        if elem in refdict:
            idx_all += [j]
            idx_rp[elem] = i
            i += 1

    # remove genes that are not in RP table
    for elem in refdict.keys():
        if elem not in idx_rp:
            del(symdict[refdict[elem]])
            del(refdict[elem])

    # filter RPs using unique gene symbols
    RP = RP[idx_all,:]
    RP = numpy.sqrt(RP)
    RP = median_norm(RP)
    RP = rowcenter(RP)
    RP_centroid,RP_label,RP_inertia = cluster.k_means(RP, NCLUSTERS, n_init=10, max_iter=300,random_state=100 )

    # read gene list
    genelist = read_list(genefile)
    reflist = [ symdict[elem] for elem in genelist if elem in symdict ]
    if len(reflist) == 0:
        reflist = [ elem for elem in genelist if elem in refdict ]
 
    # indices of genes in gene list
    idx_t = [ idx_rp[elem] for elem in reflist ]
    genestat = len(idx_all)*[False]
    for i in idx_t:
        genestat[i] = True
    genestat = numpy.array(genestat)

    # determine RP cluster most highly enriched in gene list
    p = numpy.array( [ sum( 1.0*genestat[RP_label==i])/sum(1.0*(RP_label==i) ) for i in range(NCLUSTERS) ] )
    idx_max = numpy.argmax(p)

    DHS_sample_names = [ '%s_Strength' % elem for elem in sample_names ]


    UDHS_H3K27ac_HDF5s = UDHS_H3K27ac_HDF5s.split(',')
    DHS = None
    for dhs_samplename in DHS_sample_names:
        dhs = None
        for dhs_hdf5 in UDHS_H3K27ac_HDF5s:
            if dhs_samplename in getSampleNames_hdf5(dhs_hdf5):
                dhs = read_hdf5( dhs_hdf5, [dhs_samplename] )
            else:
                continue

        if dhs is None:
            print "Please input all the hdf5 file used in regression step, with ',' as the separator"
            sys.exit()
        else:
            if DHS is None:
                DHS = dhs
            else:
                DHS =  numpy.vstack((DHS,dhs))

    DHS = DHS.transpose()
    

    DHS = numpy.sqrt(DHS)
    DHS = median_norm(DHS)
    DHS = rowcenter(DHS)

    #======== learning from TF binding data (best case) =====
    if tf_in_UDHS:
        LR_TF = linear_model.LogisticRegression(penalty='l1', tol=0.01)
        yt = 1.0*numpy.array(tf_in_UDHS)
        xt = DHS
        LR_TF.fit(xt,yt)
        DHS_pred = LR_TF.predict_proba(DHS)
        fpr, tpr, thresholds = metrics.roc_curve( tf_in_UDHS, DHS_pred[:,1], pos_label=1 )
        auc = metrics.auc(fpr,tpr)
        print >> log_fp, 'LR TF prediction',auc
    #============================================

        T = 1.0*numpy.array(tf_in_UDHS)
        T = T.reshape((T.shape[0],1))

    DHS_centroid,DHS_label,DHS_inertia = cluster.k_means(DHS, NCLUSTERS, n_init=10, max_iter=300,random_state=100 )
    print len(DHS_label)
    print DHS_label[0:10]
    #for lb in range(len(DHS_label)):
    #    print >> log3_fp, '\t'.join([str(lb+1),str(DHS_label[lb])])

    RP_dot_DHS_centroid  = numpy.dot( RP,  DHS_centroid.transpose() )
    DHS_dot_DHS_centroid = numpy.dot( DHS, DHS_centroid.transpose() )
       
    auc_gene = [] 
    print >> log2_fp, 'NCLUSTER' + '\t' + '\t'.join(sample_names)
    for i in range(NCLUSTERS):
        fpr, tpr, thresholds = metrics.roc_curve( genestat, RP_dot_DHS_centroid[:,i], pos_label=1 )
        auc = metrics.auc(fpr,tpr)
        auc_gene += [auc]
        #print >> log2_fp, 'Centroid of DHS cluster',i,'of',NCLUSTERS,':',DHS_centroid[i,:]
        print >> log_fp, 'RP AUC DHS cluster',i,'of',NCLUSTERS,':',auc
        centroids = DHS_centroid[i,:]
        centroids = ["%.3f" % x for x in centroids]
        print >> log2_fp, str(i) + '\t' + '\t'.join(centroids)
        if tf_in_UDHS:
            fpr, tpr, thresholds = metrics.roc_curve( tf_in_UDHS, DHS_dot_DHS_centroid[:,i], pos_label=1 )
            auc = metrics.auc(fpr,tpr)
            print >> log_fp, 'DHS AUC DHS cluster',i,'of',NCLUSTERS,':',auc

    auc_gene = numpy.array( auc_gene )
    k_max = numpy.argmax( auc_gene )

    T = DHS_dot_DHS_centroid[:,k_max]

    h5file = UDHS_H3K27ac_HDF5s[0]
    chrom = read_hdf5( h5file, ["chrom"] )
    start = read_hdf5( h5file, ["start"] )
    end = read_hdf5( h5file, ["end"] )
    ID = read_hdf5( h5file, ["end"] )

    out_res = []
    for i,elem in enumerate(T):
        out_res.append((chrom[i].decode('utf-8'),start[i],end[i],str(i+1),elem))
    sorted_out_res = sorted(out_res, key=lambda x: float(x[4]),reverse=True)

    print >> fpo, '%s\t%s\t%s\t%s\t%s' % ("chromosom",'start','end','UDHSID',"Score")
    for line in sorted_out_res:
        print >> fpo, '%s\t%d\t%d\t%s\t%3.2f' % (line[0],line[1],line[2],line[3],line[4])
    fpo.close()
    log_fp.close()
    log2_fp.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Predict enhancer elements from gene list and H3K27ac in union DNase-seq peaks.""")
    parser.add_argument( '-s', dest='samplefile', type=str, required=True, help='File that lists informative H3K27ac samples. One sample ID per line.' )
    parser.add_argument( '-i', dest='inputfile', type=str, required=True, help='Input gene list file. Gene symbols of refseq IDs.' )
    parser.add_argument( '-r', dest='refgeneAnn', type=str, required=True, help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>' )
    parser.add_argument( '--rp', dest='rphdf5', type=str, required=True, help='Path for the the hdf5 format rp files, if there are multiple HDF5 files, use "," as the separator' )
    parser.add_argument( '--k27ac', dest='k27achdf5', type=str, required=True, help='Path for the hdf5 format H3K27ac reads counts in UDHS regions, if there are multiple HDF5 files, use "," as the separator' )
    parser.add_argument( '-n', dest='name', type=str, required=True, help='Name of study, for output file naming.' )
    parser.add_argument( '-g', dest='genome', default='hg38', choices=['mm9','mm10','hg19','hg38'], required=False, help='genome' )
    parser.add_argument( '-t', dest='tffile', default=None, required=False, help='Indicators of TF binding at each UDHS sites 0 or 1. For performance evaluation.' )
    
    args = parser.parse_args()
    main( args.samplefile, args.inputfile, args.name, args.genome, args.rphdf5, args.k27achdf5, args.refgeneAnn, tffile=args.tffile )
