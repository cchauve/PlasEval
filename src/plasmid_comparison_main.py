__author__ = "amane"

import logging
import os

import pandas as pd
from bidict import bidict

import compare_sets
from log_errors_utils import check_file, create_directory


def get_plasmid_details(contigs_dict, filename, side, min_len):
    """Get plasmid details.

    Arguments
    ---------
    contigs_dict: Key: contig (str), Value: Nested dictionary: length (int), L_copies/R_copies (list of contig copies in plasmid set)
    path to input file
    side ('L' or 'R')

    Returns
    -------
    plasmids: list of list of contig ids
    updated contigs_dict
    plasmids_keys: bidict of plasmid indices <-> plasmid names/ids
    """
    plasmids = []
    plasmids_keys = bidict()
    count = 0
    pls_ctg_df = pd.read_csv(filename, sep="\t")

    for _, row in pls_ctg_df.iterrows():
        plasmid, contig, length = (
            f"{side}_{row['plasmid']}",
            str(row["contig"]),
            row["contig_len"],
        )
        if length >= min_len:
            if contig not in contigs_dict:
                contigs_dict[contig] = {
                    "length": length,
                    "L_copies": [],
                    "R_copies": [],
                }
            if plasmid not in plasmids_keys:
                plasmids_keys[plasmid] = count
                plasmids.append([])
                count += 1
            pls_index = plasmids_keys[plasmid]
            plasmids[pls_index].append(contig)
            contigs_dict[contig][f"{side}_copies"].append(
                [contig, pls_index, len(plasmids[pls_index])],
            )
    return contigs_dict, plasmids_keys


def comp_mode(
    left_plasmids_file,
    right_plasmids_file,
    p,
    min_len,
    max_calls,
    output_file,
    log_file,
):
    """
    Reads input files
    Initializes plasmid dicts and stores plasmid bins for both sides
    Calls compare_sets function to compute dissimilarity between the two sides
    """
    for in_file in [left_plasmids_file, right_plasmids_file]:
        check_file(in_file)
    output_dir = os.path.dirname(output_file)
    log_dir = os.path.dirname(log_file)
    create_directory([output_dir, log_dir])
    results_file = open(output_file, "w")
    # Initialize logging
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(name)s - %(levelname)s - %(message)s",
    )
    # Reading data and saving it to a dictionary with plasmids as keys and a nested dictionary of contigs as values
    contigs_dict = {}
    pls_ids_dict = {"L": {}, "R": {}}
    contigs_dict, pls_ids_dict["L"] = get_plasmid_details(
        contigs_dict,
        left_plasmids_file,
        "L",
        min_len,
    )
    contigs_dict, pls_ids_dict["R"] = get_plasmid_details(
        contigs_dict,
        right_plasmids_file,
        "R",
        min_len,
    )
    compare_sets.run_compare_plasmids(
        contigs_dict,
        pls_ids_dict,
        p,
        max_calls,
        results_file,
    )
