from aart_func import *
from params import *
import argparse


parser = argparse.ArgumentParser(description='Intensity as function of radial distance from black hole')
parser.add_argument('--nu', default=ilp.kw_nu0.value, type=float)
parser.add_argument('--mass', default=ilp.kw_mass.value, type=float)
parser.add_argument('--scaleh', default=ilp.kw_scale_height, type=float)
parser.add_argument('--thetab', default=ilp.kw_theta_b.value, type=float)
parser.add_argument('--beta', default=ilp.kw_beta, type=float)
parser.add_argument('--rie', default=ilp.kw_beta, type=float)
parser.add_argument('--rb0', default=ilp.kw_rb_0, type=float)
parser.add_argument('--nth0', default=ilp.kw_n_th0.value, type=float)
parser.add_argument('--te0', default=ilp.kw_t_e0.value, type=float)
parser.add_argument('--b0', default=ilp.kw_bv_0.value, type=float)
parser.add_argument('--pdens', default=ilp.kw_p_dens, type=float)
parser.add_argument('--ptemp', default=ilp.kw_p_temp, type=float)
parser.add_argument('--pmag', default=ilp.kw_p_bv, type=float)
parser.add_argument('--nscale', default=ilp.kw_nscale, type=float)
parser.add_argument('--emodelkey', default=0, type=int)
parser.add_argument('--bkey', default=0, type=int)
parser.add_argument('--nnoisykey', default=0, type=int)
parser.add_argument('--tnoisykey', default=0, type=int)
parser.add_argument('--bnoisykey', default=0, type=int)
parser.add_argument('--thetabkey', default=0, type=int)
parser.add_argument('--lband', default="0", type=str)   # PATH_DIR
parser.add_argument('--rtray', default="0", type=str)   # PATH_DIR
parser.add_argument('--magang', default="0", type=str)  # PATH_DIR


args = parser.parse_args()
brightparams = {
	"nu0": args.nu*ilp.Hz,				# 0
	"mass": args.mass*ilp.grams,  		# 1
	"scale_height": args.scaleh,  		# 2
	"theta_b": args.thetab*ilp.rads,  	# 3
	"beta": args.beta,  				# 4
	"r_ie": args.rie, 					# 5
	"rb_0": args.rb0,  					# 6
	"n_th0": args.nth0*ilp.cmcubed, 	# 7
	"t_e0": args.te0*ilp.kelv,  		# 8
	"b_0": args.b0*ilp.gauss,		# 9
	"p_dens": args.pdens,				# 10
	"p_temp": args.ptemp, 				# 11
	"p_mag": args.pmag,  				# 12
	"nscale": args.nscale, 				# 13
}


funckeys = {
	"emodelkey": args.emodelkey, 		# 0
	"bkey": args.bkey, 					# 1
	"nnoisykey": args.nnoisykey, 		# 2
	"tnoisykey": args.tnoisykey, 		# 3
	"bnoisykey": args.bnoisykey, 		# 4
	"theta_bkey": args.thetabkey        # 0 for variable theta b, 1 for fixed
}

# Getting angles
# fnrays="./Results/Rays_a_%s_i_%s.h5"%(spin_case,i_case)
rtray = args.rtray
lband = args.lband
magAng = args.magang
if rtray == "0":
	print("using default rtray")
	fnrays = path + "Rays_a_%s_i_%s.h5"%(spin_case,i_case)
else:
	print("using altered rtray")
	fnrays = rtray

print("Reading file: ",fnrays)

h5f = h5py.File(fnrays,'r')
phi012 = [
	h5f['phi0'][:],
	h5f['phi0'][:],
	h5f['phi0'][:]
	]

h5f.close()
#-------------------------
if magAng == "0":
	print("using default magAng")
	fn = path + "MagneticAngle_a_%s_i_%s.h5"%(spin_case,i_case)
else:
	print("using altered rtray")
	fn = magAng

h5f = h5py.File(fn,'r')

#Points for the boundary of the BH shadow
anglen0=h5f['cos2angB_n0'][:]
anglen1=h5f['cos2angB_n1'][:]
anglen2=h5f['cos2angB_n2'][:]

h5f.close()

#-------------------------
print("Intensity")


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

if bvapp!=1:

	supergrid1=h5f['grid1'][:]
	mask1=h5f['mask1'][:]
	N1=int(h5f["N1"][0])

	supergrid2=h5f['grid2'][:]
	mask2=h5f['mask2'][:]
	N2=int(h5f["N2"][0])

	fnbands=fnrays

	print("Reading file: ",fnbands)

	h5f = h5py.File(fnbands,'r')

	rs0=h5f['rs0'][:]
	sign0=h5f['sign0'][:]
	rs1=h5f['rs1'][:]
	sign1=h5f['sign1'][:]
	rs2=h5f['rs2'][:]
	sign2=h5f['sign2'][:]
	h5f.close()

	obsint.br(supergrid0,mask0,N0,rs0,sign0,anglen0,supergrid1,mask1,N1,rs1,sign1,anglen1,
			  supergrid2,mask2,N2,rs2,sign2,anglen2,brightparams,funckeys,phi012)
else:

	h5f.close()

	fnrays=rtray
	print("Reading file: ",fnrays)

	h5f = h5py.File(fnrays,'r')

	rs0_bv=h5f['rs0_bv'][:]
	sign0_bv=h5f['sign0_bv'][:]

	h5f.close()

	obsint.br_bv(supergrid0,mask0,N0,rs0_bv,sign0_bv)

