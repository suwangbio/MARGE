import argparse,math,numpy,os,sys,tables
import time

def write_hdf5( inputfiles, motif_file_name, idx=1 ):
    """ Save all motif data in current directory into hdf5 file."""    

    h5file = tables.open_file( motif_file_name, "a", driver="H5FD_CORE")

    for inputfile in inputfiles:
        fname = os.path.split(inputfile)[1]
        name = fname.split('.')[0]
        ext  = fname.split('.')[1]
        if ext != 'txt':
            pass 
        l = []
        name = name.replace('::','_')
        try:
            a = h5file.get_node("/", name )
        except:
            fp = open(inputfile)
            for i,line in enumerate(fp.readlines()):
                f  = line.split()
                l += [float(f[idx])]
            A = numpy.array(l)
            h5a = h5file.create_array( h5file.root, name, A )
            h5a.flush()

    h5file.close()
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Compute motif enrichment from pre-calculated data e.g. DHS peaks.""")
    parser.add_argument( '-i', dest='colindex', type=int, default=1, required=False, help='Column with values (starting at 0) for generating hdf5 file.' )
    parser.add_argument( '-n', dest='name', required=True, help='The prefix of the output hdf5 file.' )
    parser.add_argument( '-f', dest='files', required=True, nargs="+", help='The files that want to create the hdf5 format file.' )
    
    args = parser.parse_args()
    output = args.name
    print args.files
    write_hdf5( args.files, output, idx=args.colindex )


