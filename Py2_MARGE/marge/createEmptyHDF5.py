import argparse,math,numpy,os,sys,tables
import time

def write_hdf5( motif_file_name ):
    """ Save all motif data in current directory into hdf5 file."""    

    h5file = tables.open_file( motif_file_name, "a", driver="H5FD_CORE")
    h5file.close()
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Create an empty hdf5 file.""")
    parser.add_argument( '-n', dest='name', required=True, help='The name of the output hdf5 file.' )
    
    args = parser.parse_args()
    output = args.name
    write_hdf5( output)


