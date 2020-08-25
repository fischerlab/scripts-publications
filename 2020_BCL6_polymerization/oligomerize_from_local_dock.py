__author__ = "Shourya S. Roy Burman"
__email__ = "ssrburman@gmail.com"


import os
from pymol import cmd
import argparse

def fibre_formation_AC(pdb_dir, model, oligo_dir, n_dimers, chains):
    """
    Aligns chain A to C and B to D for successive tetramers.
    """

    # loading and aligning tetramers
    cmd.load(os.path.join(pdb_dir, model), "{}_00".format(model[:-7]))
    for ii in range(1, n_dimers):
        cmd.load(os.path.join(pdb_dir, model), "{}_{:02d}".format(model[:-7], ii))
        cmd.cealign("{}_{:02d} and chain C".format(model[:-7], ii - 1), "{}_{:02d} and chain A".format(model[:-7], ii))


    # making fibre structure
    cmd.remove("Chain C+D")

    for ii in range(1, n_dimers):
        cmd.alter("{}_{:02d} and chain A".format(model[:-7], ii), "chain='{}'".format(chains[2 * ii]))
        cmd.alter("{}_{:02d} and chain B".format(model[:-7], ii), "chain='{}'".format(chains[2 * ii + 1]))
    cmd.select("fibre", "all")
    cmd.save(os.path.join(oligo_dir, "fibre_AC_{}".format(model)), "fibre")
    cmd.delete("all")


def fibre_formation_AD(pdb_dir, model, oligo_dir, n_dimers, chains):
    """
    Aligns chain A to D and B to C for successive tetramers.
    """

    # loading and aligning tetramers
    cmd.load(os.path.join(pdb_dir, model), "{}_00".format(model[:-7]))
    for ii in range(1, n_dimers):
        cmd.load(os.path.join(pdb_dir, model), "{}_{:02d}".format(model[:-7], ii))
        cmd.cealign("{}_{:02d} and chain D".format(model[:-7], ii - 1), "{}_{:02d} and chain A".format(model[:-7], ii))


    # making fibre structure
    cmd.remove("Chain C+D")

    for ii in range(1, n_dimers):
        cmd.alter("{}_{:02d} and chain B".format(model[:-7], ii), "chain='{}'".format(chains[2 * ii]))
        cmd.alter("{}_{:02d} and chain A".format(model[:-7], ii), "chain='{}'".format(chains[2 * ii + 1]))
    cmd.select("fibre", "all")
    cmd.save(os.path.join(oligo_dir, "fibre_AD_{}".format(model)), "fibre")
    cmd.delete("all")



def main():

    # Parsing all arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="Path to directory of local docking results", required=True)
    parser.add_argument("-p", "--prefix", help="What the filenames start with", required=True)
    parser.add_argument("-s", "--suffix", help="Extension the filenames end with", required=True)
    parser.add_argument("-o", "--oligo_dir", help="Path to docking directory where fibres are stored", required=True)
    parser.add_argument("-n", "--n_dimers", help="Path to docking directory", type=int, default=26)

    args = parser.parse_args()
    pdb_dir = args.input_dir
    prefix = args.prefix
    suffix = args.suffix
    oligo_dir = args.oligo_dir
    n_dimers = args.n_dimers

    # That's the most number of chains a PDB will hold
    chain_indices = [ii for ii in range(65, 91)] + [jj for jj in range(97, 123)]
    chains = [chr(ii) for ii in chain_indices]

    # doing this for all PDBs
    for model in os.listdir(pdb_dir):
        if not model.startswith("._"):
            if model.startswith(prefix) and model.endswith(suffix):
                fibre_formation_AC(pdb_dir, model, oligo_dir, n_dimers, chains)
                fibre_formation_AD(pdb_dir, model, oligo_dir, n_dimers, chains)


main()
