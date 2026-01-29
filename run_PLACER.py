#!/usr/bin/env python3

import sys, os
import warnings
warnings.filterwarnings("ignore")
import time
import glob
import argparse
import json
import torch
from openbabel import openbabel
openbabel.obErrorLog.SetOutputLevel(0)

DIR = os.path.dirname(__file__)
sys.path.insert(0, DIR)
import PLACER
import modules.cdr_utils as cdr_utils


def main(args):
    ########################################################
    # 0. Load the network
    ########################################################
    placer = PLACER.PLACER(args.weights)
    if args.maxatoms is not None:
        placer._PLACER__params['DATALOADER']['featurizer']['maxatoms'] = args.maxatoms
    try:
        crop_size = placer._PLACER__params['DATALOADER']['featurizer']['maxatoms']
        print(f"Crop size set to: {crop_size}")
    except (KeyError, AttributeError, TypeError):
        print("Crop size: default")


    ########################################################
    # 1. Gather input files
    ########################################################
    if args.idir is not None:
        fnames = glob.glob(os.path.join(args.idir, '*.pdb'))
        if len(fnames) < 1:
            sys.exit(f"Error: no .pdb files found in '{args.idir}'")
    elif args.ifile is not None:
        if any([args.ifile.endswith(ext) for ext in ['.pdb', '.ent', '.cif', '.cif.gz']]):
            fnames = [args.ifile]
        else:
            with open(args.ifile) as f:
                fnames = [line.strip() for line in f.readlines()]
            if len(fnames) < 1:
                sys.exit(f"Error: no .pdb files found in '{args.ifile}'")
    else:
        sys.exit("Error: One of -i/--idir or -f/--ifile must be provided.")

    print(f"# number of PDB files to process: {len(fnames)}")


    ########################################################
    # 2. Parse input arguments into PLACER input object
    ########################################################
    placer_input = PLACER.PLACERinput()

    if args.exclude_common_ligands:
        placer_input.skip_ligands(PLACER.utils.get_common_ligands())

    if args.ligand_file is not None:
        ligand_ref = {lr.split(":")[0]: lr.split(":")[1] for lr in args.ligand_file}
        placer_input.ligand_reference(ligand_ref)

    if args.ignore_ligand_hydrogens:
        placer_input.ignore_ligand_hydrogens(True)

    if not args.use_sm:
        placer_input.exclude_sm(True)

    if args.poly_ligand_chains is not None:
        placer_input.poly_ligand_chains(args.poly_ligand_chains)

    if args.fixed_ligand_noise is not None:
        placer_input.fixed_ligand_noise(args.fixed_ligand_noise)

    if args.cdr_file is not None:
        cdr_def = cdr_utils.parse_cdr_definition(args.cdr_file)
        placer_input.cdr_residues(cdr_def)

    def evaluate_pred_fix_ligand_input(ligands):
        fixed_ligands = []
        for lig in ligands:
            parts = lig.split("-")
            if len(parts) == 2:
                fixed_ligands.append((parts[0], int(parts[1])))
            elif len(parts) == 3:
                fixed_ligands.append((parts[0], parts[1], int(parts[2])))
            else:
                sys.exit(f"Invalid fixed/predict ligand input: {lig}")
        return fixed_ligands

    if args.fixed_ligand is not None:
        placer_input.fixed_ligand(evaluate_pred_fix_ligand_input(args.fixed_ligand))

    if args.predict_ligand is not None:
        placer_input.predict_ligand(evaluate_pred_fix_ligand_input(args.predict_ligand))

    if args.predict_multi:
        placer_input.predict_multi(True)

    if args.target_res is not None:
        target_res = args.target_res.split("-")
        if len(target_res) == 2:
            placer_input.target_res((target_res[0], int(target_res[1])))
        elif len(target_res) == 3:
            placer_input.target_res((target_res[0], int(target_res[1]), target_res[2]))

    if args.bonds is not None:
        bonds = []
        for bond in args.bonds:
            a, b, bondlen = bond.split(':')
            a, b = a.split('-'), b.split('-')
            bonds.append([(a[0], int(a[1]), a[2], a[3]), (b[0], int(b[1]), b[2], b[3]), float(bondlen)])
        placer_input.bonds(bonds)

    if args.mutate is not None:
        mutate_dict = {}
        for mutres in args.mutate:
            pos, resn = mutres.split(':')
            resno = "".join([c for c in pos if c.isdigit()])
            chain = pos[len(resno):]
            mutate_dict[(chain, int(resno))] = resn
        placer_input.mutate(mutate_dict)

    if args.residue_json is not None:
        placer_input.add_custom_residues(json.load(open(args.residue_json)))

    if args.crop_centers is not None:
        _centers = [(c.split("-")[0], int(c.split("-")[1]), c.split("-")[2], c.split("-")[3]) for c in args.crop_centers]
        placer_input.crop_centers(_centers)

    if args.corruption_centers is not None:
        _centers = [(c.split("-")[0], int(c.split("-")[1]), c.split("-")[2], c.split("-")[3]) for c in args.corruption_centers]
        placer_input.corruption_centers(_centers)


    ########################################################
    # 3. Generate models
    ########################################################
    tic = time.time()
    for fname in fnames:

        label = os.path.basename(fname)
        for ext in [".cif.gz", ".cif", ".pdb"]:
            if label.endswith(ext):
                label = label.replace(ext, "")
        if args.suffix:
            label += f"_{args.suffix}"

        outfile_prefix = os.path.join(args.odir, label)
        if args.cautious and os.path.exists(outfile_prefix + ".csv"):
            print(f"{outfile_prefix}.csv already exists, skipping.")
            continue

        placer_input_iter = placer_input.copy()
        placer_input_iter.name(label)

        if fname.endswith(".pdb"):
            placer_input_iter.pdb(fname)
        else:
            placer_input_iter.cif(fname)

        # --- CDR mode: disable backbone noise for non-CDR residues ---
        if hasattr(placer_input_iter, 'cdr_residues') and placer_input_iter.cdr_residues() is not None:
            sigma_bb = 0.0
        else:
            sigma_bb = placer._PLACER__params.get('sigma_bb', 1.0)

        # Optional: print CDR vs fixed counts
        try:
            ligands_in_chains, fixed_ligands = PLACER.utils.parse_fixed_ligand_input(placer_input_iter, placer_input_iter.chains)
            print(f"{label}: Predicting {len(ligands_in_chains)} CDR residues; fixing {len(fixed_ligands)} residues")
        except Exception:
            pass

        # Run PLACER
        outputs = placer.run(placer_input_iter, args.nsamples, sigma_bb=sigma_bb)

        # Save outputs
        os.makedirs(args.odir, exist_ok=True)
        PLACER.protocol.dump_output(output_dict=outputs, filename=outfile_prefix, rerank=args.rerank)

    print(f"Finished predicting {len(fnames)} structures in {(time.time() - tic):.2f} seconds.")


if __name__ == "__main__":
    rank_options = ['prmsd', 'plddt', 'plddt_pde']

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--idir', type=str, help='Input folder with PDB/mmCIF files')
    parser.add_argument('-f','--ifile', type=str, help='Single PDB/mmCIF file or file containing list of files')
    parser.add_argument('-o','--odir', type=str, default="./", help='Output folder')
    parser.add_argument('-n','--nsamples', type=int, default=10, help='Number of samples to generate')
    parser.add_argument('--suffix', type=str, help='Suffix for output files')
    parser.add_argument('--cautious', action='store_true', default=False, help='Skip if output exists')
    parser.add_argument('--weights', type=str, default=f"{DIR}/weights/PLACER_model_1.pt", help='Weights file')
    parser.add_argument('--cdr-file', type=str, help='CDR definition file')
    parser.add_argument('--predict_multi', action='store_true', default=False)
    parser.add_argument('--fixed_ligand', nargs="+", type=str)
    parser.add_argument('--predict_ligand', nargs="+", type=str)
    parser.add_argument('--exclude_common_ligands', action='store_true', default=False)
    parser.add_argument('--ignore_ligand_hydrogens', action='store_true', default=False)
    parser.add_argument('--use_sm', action='store_true', default=True)
    parser.add_argument('--no-use_sm', dest='use_sm', action='store_false')
    parser.add_argument('--target_res', type=str)
    parser.add_argument('--bonds', nargs="+", type=str)
    parser.add_argument('--mutate', nargs="+", type=str)
    parser.add_argument('--crop_centers', nargs="+", type=str)
    parser.add_argument('--corruption_centers', nargs="+", type=str)
    parser.add_argument('--residue_json', type=str)
    parser.add_argument('--fixed_ligand_noise', type=float)
    parser.add_argument('--maxatoms', type=int)
    parser.add_argument('--ligand_file', nargs="+", type=str)
    parser.add_argument('--rerank', type=str, choices=rank_options)

    args = parser.parse_args()
    main(args)
