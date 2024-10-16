from aart_func import *
from params import * 
import argparse


parser = argparse.ArgumentParser(description="Magnetic Angle")
parser.add_argument('--lband', default="0", type=str)   # PATH_DIR
parser.add_argument('--rtray', default="0", type=str)   # PATH_DIR
args = parser.parse_args()

rtray = args.rtray
lband = args.lband

print("Magnetic Angle")
if lband == "0":
	fnbands = path + "LensingBands_a_%s_i_%s.h5" % (spin_case, i_case)
	print("using default lband")
else:
	print("using altered lband")
	fnbands = lband
print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
mask0=h5f['mask0'][:]
N0=int(h5f["N0"][0])

supergrid1=h5f['grid1'][:]
mask1=h5f['mask1'][:]
N1=int(h5f["N1"][0])

supergrid2=h5f['grid2'][:]
mask2=h5f['mask2'][:]
N2=int(h5f["N2"][0])

h5f.close()



if rtray == "0":
	print("using default rtray")
	fnbands = path + "Rays_a_%s_i_%s.h5"%(spin_case,i_case)
else:
	print("using altered rtray")
	fnbands = rtray
print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

rs0=h5f['rs0'][:]
sign0=h5f['sign0'][:]

rs1=h5f['rs1'][:]
sign1=h5f['sign1'][:]

rs2=h5f['rs2'][:]
sign2=h5f['sign2'][:]

h5f.close()

filename=path+"MagneticAngle_a_%s_i_%s.h5"%(spin_case,i_case)

h5f = h5py.File(filename, 'w')

h5f.create_dataset('cos2angB_n0',data=magneticf.cos2angB_f(supergrid0,mask0,N0,rs0,sign0,spin_case,thetao))
h5f.create_dataset('cos2angB_n1',data=magneticf.cos2angB_f(supergrid1,mask1,N1,rs1,sign1,spin_case,thetao))
h5f.create_dataset('cos2angB_n2',data=magneticf.cos2angB_f(supergrid2,mask2,N2,rs2,sign2,spin_case,thetao))

h5f.close()

print("File ",filename," created.")