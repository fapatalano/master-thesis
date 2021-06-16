#!/usr/bin/python3

"""
Author: Fabiana Patalano
Mail: fabiana.patalano97@gmail.com
Last updated: 30/04/2021

This script takes as input one or more trrosetta output file and discriminate interchain and intrachain connections
according to the tertiary structure of the monomer.
In case of files with more than 1 biological assembly, only the first one  is considered.
"""

import argparse
import itertools
import os
import tempfile
import seaborn as sns
import matplotlib as plt
import joblib
import numpy as np
from Bio import PDB

from retrieve_pdb import retrieve_pdb_file
import performance


def get_args():
    parser = argparse.ArgumentParser(description=" ".join(__doc__.splitlines()[4:]),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', type=str, required=False,
                        dest='native', help="true PDB file")

    parser.add_argument('-r', type=str, required=True,
                        dest='trrosetta', help="trRosetta prediction")
    parser.add_argument('-n', required=False, nargs='+', type=int, dest='relax_parameter', default=[0, 1, 2],
                        help='relax parameter')
    args = parser.parse_args()

    return args


def get_pdb_coord(pdbfile, pdb_code):
    """

    :param pdbfile:
    :param pdb_code:
    :return: coords:
             list for each chain of residues in the PDB file without Hetero atoms
    """

    parser = PDB.PDBParser()
    structure = parser.get_structure(pdb_code, pdbfile)
    chains = PDB.Selection.unfold_entities(structure[0], 'C')
    coords = [PDB.Selection.unfold_entities(chains[0], 'R'), PDB.Selection.unfold_entities(chains[1], 'R')]
    return coords


def trrosetta2maps(trrosetta, probcut=0.5, bin_step=0.5):
    """
    :param trrosetta:
    :param probcut:
    :param bin_step:
    :return:
    """

    a = trrosetta['dist']
    bins = np.array([2.25 + bin_step * i for i in range(36)])
    if len(a[0, 0, :]) != 37:
        raise Exception('ERROR! This is not a trRosetta prediction')
    D = np.full((len(a), len(a)), 20.0)  # res
    np.fill_diagonal(D, 0)
    for i in range(len(a)):
        for j in range(len(a[0])):
            maxprob_value = np.sum(a[i, j, 5:], axis=-1)
            if maxprob_value > probcut:
                D_slice = a[i, j, 1:]
                mean_dist = np.sum(np.multiply(bins, D_slice / np.sum(D_slice)))
                D[i, j] = mean_dist
    return D


def pdb2dmap(chains):
    """
    :param chains:
           list for each chain of residues in the PDB file without Heteroatoms
    :return: pdb_map:
             numpy array filled with the true intrachain and interchain distances
    """

    residues_tot = chains[0] + chains[1]
    L = len(chains[0]) * 2
    pdb_map = np.full((L, L), np.nan)
    for i in residues_tot:
        for j in residues_tot:
            if i.resname == 'GLY' and j.resname == 'GLY':
                if i.has_id('CA') and j.has_id('CA') and i.id[1] - 1 < L and j.id[1] - 1 < L:
                    pdb_map[i.id[1] - 1, j.id[1] - 1] = i["CA"] - j["CA"]
            elif i.resname == 'GLY' and j.has_id('CB'):
                if i.has_id('CA') and i.id[1] - 1 < L and j.id[1] - 1 < L: pdb_map[i.id[1] - 1, j.id[1] - 1] = i["CA"] - \
                                                                                                               j["CB"]
            elif j.resname == 'GLY' and i.has_id('CB'):
                if j.has_id('CA') and i.id[1] - 1 < L and j.id[1] - 1 < L: pdb_map[i.id[1] - 1, j.id[1] - 1] = i["CB"] - \
                                                                                                               j["CA"]
            elif i.has_id('CB') and j.has_id('CB') and i.id[1] - 1 < L and j.id[1] - 1 < L:
                pdb_map[i.id[1] - 1, j.id[1] - 1] = i["CB"] - j["CB"]
    # pdb_map[pdb_map > 20] = np.nan
    return pdb_map


def get_true_intrachain_dist(pdb_cmap, pdb_coords):
    """
    :param pdb_coords:
    :param pdb_cmap:
            PDB coordinates
    :return: distances: intra-chain distances obtained from the PDB file
        Intrachain contact between two residues i and j is said to exist if the Euclidean distance between the
        respective C beta
        (C alpha for glycine) atoms of residues i and j is less than or equal to 12.0 Å
    """

    L = len(pdb_coords[0])
    distance_matrix = pdb_cmap[:L, :L]
    distance_matrix[distance_matrix > 8] = np.nan
    return distance_matrix


def filter_contact(predicted_map, true_intrachain_c):
    """
    :param predicted_map: predicted contact map
    :param true_intrachain_c:
    :return: filtered_dist

    filter contact removing the matching true intrachain contact and remove contact with a
    sequence distance lower than 6
    """
    seq_dist = 6
    true_contact = np.where(~np.isnan(true_intrachain_c))
    listOfCoordinates = list(zip(true_contact[0], true_contact[1]))
    # Short sequence separation
    for row in range(len(predicted_map)):
        i = 0 if (row - seq_dist + 1) < 0 else row - seq_dist + 1
        j = row + seq_dist if (row + seq_dist) <= (len(predicted_map)) else (len(predicted_map))
        predicted_map[row, i:j] = np.nan
    # remove matching true contacts
    for coord in listOfCoordinates:
        if coord[0] < len(predicted_map) and coord[1] < len(predicted_map): predicted_map[coord] = np.nan
    predicted_map[predicted_map >= 12] = np.nan
    return predicted_map


def relax_removal(filtered_dist, intra_dist_pdb, n):
    """
    :param filtered_dist:
           a numpy array containing filtered inter-residue distances
    :param intra_dist_pdb:
           numpy array filled with the real intra-chain distances obtained from the pdb file
    :param n:
           relax parameter it could range from 0 to 2
    :return:  relax_removal
            If position (i,j) is a true intrachain contact, and the relax parameter is n (where, n=0, 1, and 2),
            then let X = [i-n,i+n] and Y= [j-n, j+n]. Remove all contacts (Xp,Yq) from the predicted contact map
            where Xp={i-n, i-n+1, ..., i+n} and
            Yq={j-n, j-n+1, ..., j+n}. This removes (sets to zero) a square matrix of dimension n x n centered at (i,
            j) from the predicted intrachain contact map
    """

    relaxed_d_matrix = np.pad(filtered_dist, ((n, n), (n, n)), 'constant', constant_values=np.nan)
    true_contact = np.where(~np.isnan(intra_dist_pdb))
    listOfCoordinates = list(zip(true_contact[0], true_contact[1]))
    for i,j in listOfCoordinates:
        relaxed_d_matrix[i - n:i + n + 1, j - n:j + n + 1] = np.nan
    if n != 0:
        return relaxed_d_matrix[n:-n, n:-n]
    else:
        return relaxed_d_matrix


def get_interacting_residues(model, r_cutoff=6):
    """
    Return residue-residue interactions between all chains in the model.

    Parameters
    ----------
    :param r_cutoff:
    :param model:

    Returns
    -------
    dict
        A dictionary of interactions between chains i (0..n-1) and j (i+1..n).
        Keys are (chain_idx, chain_id, residue_idx, residue_resnum, residue_amino_acid) tuples.
        (e.g. (0, 'A', 0, '0', 'M'), (0, 1, '2', 'K'), ...)
        Values are a list of tuples having the same format as the keys.

    Examples
    --------
    You can reverse the order of keys and values like this::

        complement = dict()
        for key, values in get_interacting_chains(model):
            for value in values:
                complement.setdefault(value, set()).add(key)


    You can get a list of all interacting chains using this command::

        {(key[0], value[0])
         for (key, values) in get_interacting_chains(model).items()
         for value in values}


    """
    interactions_between_chains = dict()
    atoms_c2 = PDB.Selection.unfold_entities(model[1], 'A')
    ns = PDB.NeighborSearch(atoms_c2)

    # Residue 1
    for residue_1 in model[0]:
        interacting_residues = set()
        for atom_1 in residue_1: interacting_residues.update(ns.search(atom_1.get_coord(), r_cutoff, 'R'))
        # ns search return a list of atoms at distance x from the target atom

        # Residue 2
        interacting_residue_ids = []
        for residue_2 in interacting_residues:
            interacting_residue_ids.append(residue_2)
            if interacting_residue_ids:
                interactions_between_chains.setdefault(residue_1, set()).update(interacting_residue_ids)
    return interactions_between_chains


def get_true_inter_distances(interactions_between_chains, coords):
    """
    :param coords:
    :param interactions_between_chains:
    :return: inter
    """
    L = len(coords[0])
    inter = np.full((L, L), np.nan)
    for key, values in interactions_between_chains.items():
        for value in values:
            if key.resname == 'GLY' and value.resname == 'GLY':
                i = key.id[1] - 1
                j = (value.id[1] - 1) - len(coords[0])
                if i < L and j < L: inter[i, j] = key['CA'] - value['CA']
            elif key.resname == 'GLY' and value.has_id('CB'):
                i = key.id[1] - 1
                j = (value.id[1] - 1) - len(coords[0])
                if i < L and j < L: inter[i, j] = key['CA'] - value['CB']
            elif value.resname == 'GLY' and key.has_id('CB'):
                i = key.id[1] - 1
                j = (value.id[1] - 1) - len(coords[0])
                if i < L and j < L: inter[i, j] = key['CB'] - value['CA']
            elif key.has_id('CB') and value.has_id('CB'):
                i = key.id[1] - 1
                j = (value.id[1] - 1) - len(coords[0])
                if i < L and j < L: inter[i, j] = key['CB'] - value['CB']

    return inter


def disteval_main(trrosetta, relax_parameter):
    tot_dict = {}
    tot_true_intra_distances = []
    tot_true_inter_distances = []
    tot_pred_intra_distances = []
    tot_pred_inter_distances_filter = []
    tot_pred_inter_distances_relax = {0: [], 1: [], 2: []}

    for entry in os.scandir(trrosetta):
        if entry.path.endswith(".npz") and entry.is_file():
            performance_dict = {}
            # Prepare the monomer structure and the prediction

            pdb_code = entry.name.split('-')[0][:-4].upper()
            print('\nStructure', pdb_code)
            pdbfile = f'/home/fapatalano/Desktop/tesi/dataset/pdb_file/{pdb_code}.pdb'
            residues = get_pdb_coord(pdbfile, pdb_code)
            pdb_cmap = pdb2dmap(residues)
            L = pdb_cmap.shape[0] // 2
            tot_true_intra_distances += pdb_cmap[:L, :L].flatten().tolist()
            tot_true_inter_distances += pdb_cmap[:L, L:].flatten().tolist()

            print('\nTrue dmap: ', pdb_cmap.shape)
            true_intrachain_contact = get_true_intrachain_dist(pdb_cmap, residues)

            print('\nLoad the input trRosetta prediction..')
            trrosetta_map = np.load(entry)
            D = trrosetta2maps(trrosetta_map)
            print('TrRosetta map: ', D.shape)
            tot_pred_intra_distances += D.flatten().tolist()

            interactions_between_chains = get_interacting_residues(residues, r_cutoff=6)
            true_interchain_contact = get_true_inter_distances(interactions_between_chains, residues)
            print('Total number of interchain contacts', len(np.where(~np.isnan(true_interchain_contact))[0]))

            # Filter short range contacts and true intrachain contacts
            filtered_dist = filter_contact(D, true_intrachain_contact)
            tot_pred_inter_distances_filter += filtered_dist.flatten().tolist()

            # Relax removal
            relax_parameter_list = list(itertools.product(relax_parameter, repeat=2))
            for n in relax_parameter_list:
                relax_removal_map = relax_removal(filtered_dist, true_intrachain_contact, n[0])
                tot_pred_inter_distances_relax[n[0]] += relax_removal_map.flatten().tolist()

                tp, tn, fp, fn = performance.relaxation(true_interchain_contact, relax_removal_map, n[1])
                performance_dict[n] = [tp, tn, fp, fn]
                # print(f'TP={tp}\nTN={tn}\nFP={fp}\nFN={fn}')
    return tot_true_intra_distances, tot_true_inter_distances, tot_pred_intra_distances, \
           tot_pred_inter_distances_filter, tot_pred_inter_distances_relax


if __name__ == "__main__":
    args = get_args()
    print(args)

    native = None
    pdb_code = None

    relax_parameter = args.relax_parameter
    trrosetta = os.path.abspath(args.trrosetta)

    if args.native is not None:
        native = os.path.abspath(args.native)
        pdb_code = os.path.basename(native).split('.')[0]

    disteval_main(trrosetta=trrosetta,
                  relax_parameter=relax_parameter)
