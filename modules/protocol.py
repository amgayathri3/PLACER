import os, sys
import numpy as np
import pandas as pd
import itertools
import random
from openbabel import openbabel
openbabel.obErrorLog.SetOutputLevel(0)
import utils


# =============================================================================================================
### Utility functions ###
def parse_ligand_reference(ligand_reference, ccd_ligands):
    """
    Parses user-provided ligand reference info
    Arguments:
        ligand_reference (dict)
        ccd_ligands (dict)   ::  dictionary of CCD ligands. Commonly retrieved from PLACER instance with mols() method.
    """
    ## Reading in user-defined reference ligand file(s)
    assert isinstance(ligand_reference, dict)

    reference_ligands = {}
    for lig_name in ligand_reference:
        lig_file = ligand_reference[lig_name]
        
        if lig_file == "CCD":
            assert lig_name in ccd_ligands.keys(), f"{lig_name} not in internal CCD database!"
            print(f"Reading ligand {lig_name} reference definition from internal CCD database.")
            lig_file = ccd_ligands[lig_name]["sdf"]
            ext = "sdf"
        else:
            if os.path.exists(lig_file):
                ext = lig_file.split(".")[-1]
                assert ext in ["sdf", "mol2"], f"Invalid reference ligand file type: {ext}"
                print(f"Reading ligand {lig_name} reference definition from file: {lig_file}")
                lig_file = open(lig_file, "r").read()
            else:
                print(f"Reading ligand {lig_name} reference definition input string.")
                if "@<TRIPOS>MOLECULE" in lig_file.split("\n")[0]:
                    ext = "mol2"
                elif "  0999 V2000" in lig_file.split("\n")[3]:
                    ext = "sdf"
                else:
                    sys.exit(f"Unable to figure out reference ligand {lig_name} reference type")

        reference_ligands[lig_name] = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat(ext)
        obConversion.ReadString(reference_ligands[lig_name], lig_file)
    return reference_ligands


def parse_input_structure(input_object, ligand_reference, pdbparser, cifparser):
    ### Parsing the input structure ###

    ## PDB parsing
    if input_object.pdb() is not None:

        # mutate positions
        # TODO: could this be done on the chains object?
        if input_object.mutate() is not None:
            residues = {''.join(k.split()[::-1]):list(g)
                        for k,g in itertools.groupby(input_object.pdb().split("\n"),key=lambda x : x[21:26])}

            for pos, rname in input_object.mutate().items():
                pos_str = f"{pos[1]}{pos[0]}"
                print(f'# mutation: position {pos} to residue {rname}')
                res = pdbparser.Parser.getRes(rname)
                if res is None:
                    sys.exit(f'Residue "{rname}" is not in the library')
                if pos_str not in residues.keys():
                    sys.exit(f'Residue "{pos}" is not in input')
                residues[pos_str] = utils.mutate(residues[pos_str], rname, res['res'].atoms)
            atom_lines = [l[:6] + '%-5d'%(i) + l[11:] for i,l in enumerate(itertools.chain.from_iterable(residues.values()))]
            chains = pdbparser.parseProtein(atom_lines)
        else:
            chains = pdbparser.parseProtein(input_object.pdb().split("\n"))


        # parse the ligand to OpenBabel object
        obmol = pdbparser.parse_ligand_from_pdb_to_obmol(input_object.pdb(),
                                                         ligand_reference,
                                                         ignore_hydrogens=input_object.ignore_ligand_hydrogens())

        # if ligand is found, then convert it to chains dictionary
        if obmol is not None:
            lig_chains = pdbparser.parseLigand(obmol, pdbstr=input_object.pdb().split("\n"))
            assert not any([k in chains.keys() for k in lig_chains.keys()]), "One or more of ligand chains already exist in parsed protein chains. Consider changing the chain letters of your ligands."
            chains.update(lig_chains)

    ## CIF parsing
    if input_object.cif() is not None:
        chains,asmb,covale,meta = cifparser.parse(input_object.cif())
        obmol = None
        # TODO: need to add ligand mol/sdf parsing here too?
        # Currently only works with native RCSB CIF files
        # AF3 and Rosetta CIF files are not parsed correctly.
    return chains, obmol


#def parse_fixed_ligand_input(input_object, chains):
    ### Parsing user choices about fixed ligands and to-be-prediced ligands ###
    # Ligand can be defined in 3 ways:
    # name3 - all ligands with this name will be fixed/predicted
    # (name3, resno) - ligands with this name and residue number will be fixed/predicted.
    #                  If there are multiple chains with the same residue numbering then all copies of the ligand will be fixed/predicted.
    # (chain, name3, resno) - ligand in this chain with this name and residue number will be fixed/predicted.
    #ligands_in_chains = []
    #for ch in chains:
    #    if chains[ch].type == "nonpoly" or ch in input_object.poly_ligand_chains():  # currently not supporting fixing side chains # Modify this line to check for user-specified polymer chains
    #        ligands_in_chains += list(set([(ch, at[2], int(at[1])) for at in chains[ch].atoms]))  # (str, str, int) // (chain, name3, resno)

    #fixed_ligands = []
    #if input_object.predict_ligand() is not None:
        ## Predicting only selected ligand. Everything else is fixed.
     #   for lig in ligands_in_chains:
      #      if lig[1] in input_object.predict_ligand():  # name3
       #         continue
        #    elif (lig[1], lig[2]) in input_object.predict_ligand():  # (name3, resno)
         #       continue
          #  elif lig[0] in input_object.predict_ligand():  # chain
           #     continue
            #elif lig in input_object.predict_ligand():  # (chain, name3, resno)
             #   continue
           # fixed_ligands.append(lig)

    #elif input_object.fixed_ligand() is not None:
        ## Fixing selected ligand(s). Everything else is predicted.
        #for lig in ligands_in_chains:
        #    if lig[1] in input_object.fixed_ligand():
          #      fixed_ligands.append(lig)
         #   elif (lig[1], lig[2]) in input_object.fixed_ligand():
         #       fixed_ligands.append(lig)
         #   elif lig in input_object.fixed_ligand():
         #       fixed_ligands.append(lig)
   # return ligands_in_chains, fixed_ligands

def parse_fixed_ligand_input(input_object, chains):
    """
    Parses user choices about fixed ligands and to-be-prediced ligands
    """
    ligands_in_chains = []
    for ch in chains:
        # Normal nonpoly chains (original behavior)
        if chains[ch].type == "nonpoly":  # currently not supporting fixing side chains
            ligands_in_chains += list(set([(ch, at[2], int(at[1])) for at in chains[ch].atoms]))  # (str, str, int) // (chain, name3, resno)
            
        # Antibody chains with CDR definitions
        elif hasattr(input_object, 'poly_ligand_chains') and \
             ch in input_object.poly_ligand_chains() and \
             hasattr(input_object, 'cdr_residues') and \
             input_object.cdr_residues() and \
             ch in input_object.cdr_residues():
            cdr_res_nums = input_object.cdr_residues()[ch]
            cdr_atoms = []
            for at in chains[ch].atoms:
                try:
                    res_num = int(at[1])
                    if res_num in cdr_res_nums:
                        cdr_atoms.append((ch, at[2], res_num))
                except ValueError:
                    continue
            
            ligands_in_chains += list(set(cdr_atoms))
            print(f"Added {len(set([a[2] for a in cdr_atoms]))} CDR residues from chain {ch} as ligands")
        
        # Handle case where chain is flagged as ligand but no CDR definition
        elif hasattr(input_object, 'poly_ligand_chains') and ch in input_object.poly_ligand_chains():
            print(f"Warning: Chain {ch} marked as ligand but no CDR definitions found. Adding entire chain.")
            ligands_in_chains += list(set([(ch, at[2], int(at[1])) for at in chains[ch].atoms]))

    # fixing
    fixed_ligands = []
    if input_object.predict_ligand() is not None:
        ## Predicting only selected ligand. Everything else is fixed.
        for lig in ligands_in_chains:
            if lig[1] in input_object.predict_ligand():  # name3
                continue
            elif (lig[1], lig[2]) in input_object.predict_ligand():  # (name3, resno)
                continue
            elif lig[0] in input_object.predict_ligand():  # chain
                continue
            elif lig in input_object.predict_ligand():  # (chain, name3, resno)
                continue
            fixed_ligands.append(lig)

    elif input_object.fixed_ligand() is not None:
        ## Fixing selected ligand(s). Everything else is predicted.
        for lig in ligands_in_chains:
            if lig[1] in input_object.fixed_ligand():
                fixed_ligands.append(lig)
            elif (lig[1], lig[2]) in input_object.fixed_ligand():
                fixed_ligands.append(lig)
            elif lig in input_object.fixed_ligand():
                fixed_ligands.append(lig)
    return ligands_in_chains, fixed_ligands

def build_crop(dataloader, input_object, chains, obmol, fixed_ligands):
    """
    Generates a set of cropped atoms.
    If CDR residues are provided, the crop is centered on ALL CDR atoms
    (multi-center crop) to support simultaneous sampling of all CDRs.
    """

    # ============================================================
    # CDR MULTI-CENTER CROP OVERRIDE
    # ============================================================
    if hasattr(input_object, "cdr_residues") and input_object.cdr_residues():
        print("CDR crop override active")

        cdr_atoms = []

        for ch, resnums in input_object.cdr_residues().items():
            if ch not in chains:
                print(f"Warning: chain {ch} not found in structure")
                continue

            for at_key in chains[ch].atoms:
                try:
                    resnum = int(at_key[1])
                except ValueError:
                    continue

                if resnum in resnums:
                    atom = chains[ch].atoms[at_key]
                    if atom.element > 1 and atom.occ != 0:
                        cdr_atoms.append(atom)

        print("Total CDR atoms used as crop centers:", len(cdr_atoms))
        print("Chains contributing CDRs:", list(input_object.cdr_residues().keys()))

        if len(cdr_atoms) == 0:
            sys.exit("ERROR: CDR override requested but no CDR atoms were found")

        # Perform multi-center crop
        cropped_atoms = dataloader.dataset.dataset.get_crop(
            chains,
            cdr_atoms,
            exclude=fixed_ligands
        )

        print("Total atoms in cropped region:", len(cropped_atoms))

        # Hard sanity checks
        if len(cropped_atoms) < 500:
            print("WARNING: Crop is very small; check CDR definitions or numbering")

        if hasattr(dataloader.dataset.dataset, "params"):
            maxatoms = dataloader.dataset.dataset.params.get("maxatoms", None)
            if maxatoms is not None and len(cropped_atoms) > maxatoms:
                print("WARNING: Crop exceeds maxatoms")

        return cropped_atoms, cdr_atoms

    # ============================================================
    # ORIGINAL PLACER LOGIC BELOW (UNCHANGED)
    # ============================================================

    user_defined_center = False
    if input_object.corruption_centers() is not None or input_object.crop_centers() is not None:
        user_defined_center = True

    if input_object.pdb() is not None and input_object.target_res() is None and user_defined_center is False:
        cropped_atoms, center = dataloader.dataset.dataset.get_crop_around_mol(
            chains,
            obmol,
            exclude=fixed_ligands,
            multicenter=input_object.predict_multi()
        )

    elif (
        input_object.cif() is not None or
        (input_object.pdb() is not None and input_object.target_res() is not None) or
        (input_object.pdb() is not None and user_defined_center is True)
    ):
        skip_chains = [ch for ch in chains if chains[ch].type != "nonpoly"]
        _fixed_ligands = []

        if input_object.target_res() is not None and "polypeptide" in chains[input_object.target_res()[0]].type:
            skip_chains = [ch for ch in chains if ch != input_object.target_res()[0]]
            residues = []
            for ch in chains:
                for atm in chains[ch].atoms:
                    if atm[:3] not in residues:
                        residues.append((atm[0], atm[2], atm[1]))

            for res in residues:
                if len(input_object.target_res()) == 2:
                    if (res[0], res[2]) != input_object.target_res():
                        _fixed_ligands.append(res)
                elif len(input_object.target_res()) == 3:
                    if res != input_object.target_res():
                        _fixed_ligands.append(res)

        if input_object.corruption_centers() is None and input_object.crop_centers() is None:
            center = dataloader.dataset.dataset.get_crop_center(
                chains,
                skip_chains=skip_chains,
                exclude=fixed_ligands + _fixed_ligands,
                multicenter=input_object.predict_multi()
            )
        else:
            ligands_in_chains = []
            for ch in chains:
                if chains[ch].type == "nonpoly":
                    ligands_in_chains += list(
                        set([(ch, at[2], int(at[1])) for at in chains[ch].atoms])
                    )

            ligands_to_predict = [lig for lig in ligands_in_chains if lig not in fixed_ligands]
            center = []

            for lig in ligands_to_predict:
                if input_object.crop_centers() is not None:
                    random_center = random.choice(input_object.crop_centers())
                else:
                    random_center = random.choice(input_object.corruption_centers())

                _lig_atoms = []
                for ch in chains:
                    if ch != lig[0]:
                        continue
                    for at in chains[ch].atoms:
                        atom = chains[ch].atoms[at]
                        if atom.occ == 0 or atom.element <= 1:
                            continue
                        if (ch, at[2], int(at[1])) == lig:
                            _lig_atoms.append(atom)

                center.append(random.choice(_lig_atoms))

        cropped_atoms = dataloader.dataset.dataset.get_crop(chains, center)
    else:
        sys.exit("build_crop :: No valid input provided")

    print("Total atoms in cropped region:", len(cropped_atoms))
    return cropped_atoms, center


def dump_output(output_dict, filename, rerank=None):
    """
    Dumps PLACER outputs into a CSV file (scores) and a multimodel PDB file ().
    Outputs are reranked based on requested scoreterm (prmsd, plddt, plddt_pde).
    Returns a pandas DataFrame of the scores.

    Parameters
    ----------
    output_dict : dict
        PLACER outputs dictionary.
    filename : str
        Name prefix for the produced files.
        Scores are saved as {filename}.csv
        Models are saved as {filename}_model.pdb
    rerank : str
        Scoreterm used to rerank the models. Must be some of these: ["prmsd", "plddt", "plddt_pde"]

    """
    if rerank is not None:
        assert rerank in ["prmsd", "plddt", "plddt_pde"]
        output_dict = utils.rank_outputs(output_dict, rerank)

    with open(f"{filename}_model.pdb", "w") as file:
        for n in output_dict.keys():
            file.write(output_dict[n]["model"])
    print(f"Wrote generated models to {filename}_model.pdb")

    ignore_csv_keys = ["item", "model", "center", "Xs", "Ds", "plDDTs", "pDEVs"]
    df = pd.DataFrame.from_dict({k: [output_dict[n][k] for n in output_dict] for k in output_dict[0].keys() if k not in ignore_csv_keys})
    df.to_csv(f"{filename}.csv", index=False)
    print(f"Wrote scores to {filename}.csv")

