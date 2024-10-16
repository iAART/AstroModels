from aart_func import *
import params
from params import *

funckeys = {
        "emodelkey" : 0, # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law (1 is WIP)
        "bkey" : 2,      # Type of magnetic field profile , 0 = true function from brodrick and loeb eq. 3, 1 = power law with lmfit values of 0, 2 = power law from values set in brightparams
        "nnoisykey" : 0, # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey" : 0, # tnoisykey Inoisy temperature
        "bnoisykey" : 0, # bnoisykey Inoisy magnetic field
        "theta_bkey": 0 # Variable impact parameter, 0 for varied, 1 for fixed
}

def createIntensityArgs(brightparams,lband="0",rtray="0",magAng="0",funckeys=funckeys):
    args = ' '
    cmd1_args = {
        "nu0": '--nu ',
        "mass": '--mass ',
        "scale_height": '--scaleh ',
        "theta_b": '--thetab ',
        "beta": '--beta ',
        "r_ie": '--rie ',
        "rb_0": '--rb0 ',
        "n_th0": '--nth0 ',
        "t_e0": '--te0 ',
        "b_0": '--b0 ',
        "p_dens": '--pdens ',
        "p_temp": '--ptemp ',
        "p_mag": '--pmag ',
        "nscale": '--nscale ',
    }
    cmd2_args = {
        "emodelkey": '--emodelkey ',
        "bkey": '--bkey ',
        "nnoisykey": '--nnoisykey ',
        "tnoisykey": '--tnoisykey ',
        "bnoisykey": '--bnoisykey ',
        "theta_bkey": '--thetabkey '
    }


    # brightparams = fpp.bp_steeperT
    # funckeys = fpp.fk_fiducial

    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '

    args += "--lband " + lband + " --rtray " + rtray + " --magang " + magAng

    return 'python3 radialintensity.py' + args


def intensityNameWrite(brightparams,funckeys=funckeys):
    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}_magkey_{}.h5').format(
        params.spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"].value),
        "{:.5e}".format(brightparams["mass"].value),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"].value),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"].value),
        "{:.1e}".format(brightparams["t_e0"].value),
        "{:.3e}".format(brightparams["b_0"].value),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        float(brightparams["p_mag"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"],
        funckeys["theta_bkey"]
    )

    return filename


def intensityNameNoUnits(brightparams,funckeys=funckeys):

    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'b0_{}_pdens_{}_ptemp_{}_pmag_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}_magkey_{}.h5').format(
        params.spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"]),
        "{:.5e}".format(brightparams["mass"]),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"]),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"]),
        "{:.1e}".format(brightparams["t_e0"]),
        "{:.3e}".format(brightparams["b_0"]),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        float(brightparams["p_mag"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"],
        funckeys["theta_bkey"]
    )

    return filename

