import pandas as pd
import numpy as np
import uuid
import os
import logging
import json as json_package
import re
from .utils import get_uniprot_seqs

from .constants import UNIPROT_REST_API

from .gget_diamond import diamond, tsv_to_df


from .gget_setup import (
    ELM_FILES,
    ELM_INSTANCES_FASTA,
    ELM_CLASSES_TSV,
    ELM_INSTANCES_TSV,
)


def motif_in_query(row):
    """
    Checks if motif is in the overlapping region with the query sequence

    Args:
    row     - row in dataframe

    Returns: True if the motif is in between the target start and end of sequence. False otherwise
    """
    return (
        True
        if (row["Start in ortholog"] >= row["target_start"])
        & (row["End in ortholog"] <= row["target_end"])
        else False
    )


def get_elm_instances(UniProtID, verbose=False):
    """
    Get ELM instances and their information from local ELM tsv files

    Args:

    UniProt ID - UniProt ID to search for in the accession column of ELM tsv files
    verbose    - If True, turns on logging. Default: False

    Returns:

    df_final - dataframe combining ELM instances and information (description, functional site...)

    """

    # Check if ELM files are present
    if os.path.exists(ELM_INSTANCES_FASTA):
        if verbose:
            logging.info(f"ELM fasta file installed succesfully.")

    else:
        logging.error("ELM fasta file download failed.")
        return

    if os.path.exists(ELM_CLASSES_TSV) and os.path.exists(ELM_INSTANCES_TSV):
        if verbose:
            logging.info("ELM tsv files installed successfully.")
    else:
        logging.error("ELM tsv files download failed.")
        return

    # return matching rows from elm_instances.tsv
    df_full_instances = tsv_to_df(ELM_INSTANCES_TSV)

    df_full_instances.rename(columns={"Primary_Acc": "UniProt ID"}, inplace=True)
    df_full_instances.rename(columns={"Start": "Start in ortholog"}, inplace=True)
    df_full_instances.rename(columns={"End": "End in ortholog"}, inplace=True)
   
    # print("Uniprot ID input", UniProtID)
    # print("Matching uniprot id from instances.tsv", df_full_instances["UniProt ID"])
    df_instances_matching = df_full_instances.loc[
        df_full_instances["UniProt ID"].str.contains(UniProtID)
    ]
    # return (df_instances_matching)

    # get class descriptions from elm_classes.tsv
    df_classes = tsv_to_df(ELM_CLASSES_TSV)
    df_classes.rename(columns={"Accession": "class_accession"}, inplace=True)

    # merge two dataframes using ELM Identifier
    df = df_instances_matching.merge(df_classes, how="left", on=["ELMIdentifier"])
    # print(f"df merged orthologs columns {df.columns}")
    # reorder columns
    change_column = [
        "UniProt ID",
        "class_accession",
        "ELMIdentifier",
        "FunctionalSiteName",
        "Description",
        "Regex",
        "Probability",
        "Start in ortholog",
        "End in ortholog",
        "Query Cover",
        "Per. Ident",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "ProteinName",
        "Organism",
        "References",
        "InstanceLogic",
        "PDB",
        "#Instances",
        "#Instances_in_PDB",
    ]
    df_final = df.reindex(columns=change_column)
    return df_final


def seq_workflow(
    sequences,
    sequence_lengths,
    reference=ELM_INSTANCES_FASTA,
    out=None,
    sensitivity="very-sensitive",
    json=False,
    verbose=True,
):
    """
    Alignment of sequence using DIAMOND to get UniProt ID. Use the UniProt ID to construct an ortholog dataframe similar to the UniProt workflow
    except for additional columns for start, end and whether the motif overlaps the target sequence.

    Args:
    sequences        - list of user input amino acid sequence
    sequence_lengths - list of lengths respective to each sequence
    input_file       - Set to fasta file path (include .fa) if input contains multiple sequences. Default: None
    reference        - Set to reference file path (include .dmnd). If not specified, the ELM instances tsv file is used to construct the reference database file.
    out              - Folder name to save output files. Default: None (output is converted and returned in dataframe format. The output temporary files is not saved)
    sensitivity      - Sensitivity level to do DIAMOND alignment. The sensitivity can be adjusted using the options --fast, --mid-sensitive, --sensitive, --more-sensitive, --very-sensitive and --ultra-sensitive. Default: very-sensitive
    json             - If True, returns results in json format instead of data frame. Default: False.
    verbose          - If True, turns on logging for INFO_level messages

    Returns:
    df              - dataframe consisting of ELM instances, class information, start, end in query, and if motif overlaps with target sequence

    """
    df = pd.DataFrame()

    df_diamond = diamond(
        sequences,
        reference=reference,
        sensitivity=sensitivity,
        json=json,
        verbose=verbose,
        out=out,
    )
    # print(df_diamond)
    seq_number = 1
    for sequence, seq_len in zip(sequences, sequence_lengths):
        sequence = str(sequence)

        # If no match found for sequence, raise error

        if len(df_diamond) == 0:
            # !!! TODO change to warning
            logging.info(
                f"Sequence #{seq_number}: No orthologous proteins found in ELM database."
            )
        else:
            logging.info(
                f"Sequence #{seq_number}: Found orthologous proteins in ELM database. Retrieving data about ELMs occurring in orthologs..."
            )

            # Construct df with elm instances from uniprot ID returned from diamond
            # TODO double check that this gets info if more than one UniProt ID matched
            uniprot_ids = str(df_diamond["target_accession"]).split("|")[1]
            logging.info(
                f"Pairwise sequence alignment with DIAMOND matched the following UniProt IDs {uniprot_ids}. Retrieving ELMs for each UniProt ID..."
            )
            # print(df_diamond.columns)

            for i, uniprot_id in enumerate(df_diamond["target_accession"].values):
                # print(f"UniProt ID {uniprot_id}")
                df_elm = get_elm_instances(str(uniprot_id).split("|")[1], verbose)
                # missing motifs other than the first one
                df_elm["Query Cover"] = df_diamond["length"].values[i] / seq_len * 100
                df_elm["Per. Ident"] = df_diamond["Per. Ident"].values[i]  
                df_elm["query_start"] = int(df_diamond["query_start"].values[i])
                df_elm["query_end"] = int(df_diamond["query_end"].values[i])
                print("Target start", int(df_diamond["target_start"].values[i]))
                df_elm["target_start"] = int(df_diamond["target_start"].values[i])
                # print(df_elm["target_start"])
                print("Target end", int(df_diamond["target_end"].values[i]))
                df_elm["target_end"] = int(df_diamond["target_end"].values[i])
                # print(f"df_seq_workflow: {df_elm.columns}")
                df_elm["motif_in_query"] = df_elm.apply(motif_in_query, axis=1)

                df = pd.concat([df, df_elm])

        seq_number += 1

    return df


def regex_match(sequence):
    """
    Compare ELM regex with input sequence and return all matching elms

    Args:
    sequence - user input sequence (can be either amino acid seq or UniProt ID)

    Returns:
    df_final - dataframe containing regex matches
    TODO: Make sure this returns empty dataframe if no matches were found
    """
    # Get all motif regex patterns from elm db local file
    df_elm_classes = tsv_to_df(ELM_CLASSES_TSV)
    df_full_instances = tsv_to_df(ELM_INSTANCES_TSV)

    elm_ids = df_elm_classes["Accession"]
    # print(f"all elm ids {elm_ids}")
    regex_patterns = df_elm_classes["Regex"]

    df_final = pd.DataFrame()

    # Compare ELM regex with input sequence and return all matching elms
    for elm_id, pattern in zip(elm_ids, regex_patterns):
        regex_matches = re.finditer(pattern, sequence)

        for match_string in regex_matches:
            # print(match_string)
            elm_row = df_elm_classes[df_elm_classes["Accession"] == elm_id]
            # print(f"ELM ID {elm_id} match with pattern {pattern}")
            elm_row.insert(
                loc=1,
                column="Instances (Matched Sequence)",
                value=match_string.group(0),
            )

            (start, end) = match_string.span()
            elm_row.insert(loc=2, column="Start in query", value=str(start))
            elm_row.insert(loc=3, column="End in query", value=str(end))

            elm_identifier = [str(x) for x in elm_row["ELMIdentifier"]][0]

            df_instances_matching = df_full_instances.loc[
                df_full_instances["ELMIdentifier"] == elm_identifier
            ]

            # merge two dataframes using ELM Identifier, since some Accessions are missing from elm_instances.tsv
            elm_row = elm_row.merge(
                df_instances_matching, how="left", on=["ELMIdentifier"]
            )

            df_final = pd.concat([df_final, elm_row])

    df_final.pop("Accession_y")
    df_final.pop("#Instances")
    df_final.pop("#Instances_in_PDB")
    df_final.pop("References")
    df_final.pop("InstanceLogic")

    df_final.rename(columns={"Accession_x": "instance_accession"}, inplace=True)

    change_column = [
        "instance_accession",
        "ELMIdentifier",
        "FunctionalSiteName",
        "ELMType",
        "Description",
        "Instances (Matched Sequence)",
        "Probability",
        "Start in query",
        "End in query",
        "Methods",
        "ProteinName",
        "Organism",
    ]
    df_final = df_final.reindex(columns=change_column)

    return df_final


def elm(
    sequence,
    uniprot=False,
    json=False,
    input_file=None,
    reference=ELM_INSTANCES_FASTA,
    out=None,
    sensitivity="very-sensitive",
    verbose=True,
):
    """
    Searches the Eukaryotic Linear Motif resource for Functional Sites in Proteins.

    Args:
     - sequence       Amino acid sequence or Uniprot ID (str).
                      If Uniprot ID, set 'uniprot==True'.
     - uniprot        Set to True if input is a Uniprot ID instead of amino acid sequence. Default: False.
     - json           If True, returns results in json format instead of data frame. Default: False.
     - input_file     Set to fasta file path (include .fa) if input contains multiple sequences. Default: None
     - reference      Set to reference file path (include .dmnd). If not specified, the ELM instances tsv file is used to construct the reference database file.
     - out            Folder name to save output files. Default: None (output is converted and returned in dataframe format. The output temporary files is not saved)
     - sensitivity    Sensitivity level to do DIAMOND alignment. The sensitivity can be adjusted using the options --fast, --mid-sensitive, --sensitive, --more-sensitive, --very-sensitive and --ultra-sensitive. Default: very-sensitive
     - verbose        True/False whether to print progress information. Default: True.

    Returns two data frames: orthologs and regex matches from ELM database.
    """
    # TODO: check if handle user input file correctly since seq workflow does not have input file arg

    if not uniprot:
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZXBJZ")
        # Convert input sequence to upper case letters
        sequence = sequence.upper()

        # If sequence is not a valid amino sequence, raise error
        if not set(sequence) <= amino_acids:
            logging.warning(
                f"Input amino acid sequence contains invalid characters. If the input is a UniProt ID, please specify `uniprot=True` (python: uniprot=True)."
            )

    df = pd.DataFrame()

    # building first ortholog dataframe
    if uniprot:
        df = get_elm_instances(sequence, verbose)
        df["Query Cover"] = np.nan
        df["Per. Ident"] = np.nan

        if len(df) == 0:
            logging.warning(
                "UniProt ID does not match UniProt IDs in the ELM database. Converting UniProt ID to amino acid sequence..."
            )
            df_uniprot = get_uniprot_seqs(server=UNIPROT_REST_API, ensembl_ids=sequence)

            try:
                # only grab sequences where id match exact input uniprot id
                aa_seqs = df_uniprot[df_uniprot["uniprot_id"] == sequence][
                    "sequence"
                ].values
                seq_lens = df_uniprot["sequence_length"].values
                # print(f"aa_seqs {aa_seqs}")
            except KeyError:
                raise ValueError(
                    f"No sequences found for UniProt ID {sequence} from searching the UniProt server. Please double check your UniProt ID and try again."
                )

    if len(df) == 0:
        # add input aa sequence and its length to list
        if not uniprot:
            aa_seqs = [sequence]
            seq_lens = [len(sequence)]
            if verbose:
                logging.info(
                    f"Performing pairwise sequence alignment against ELM database using DIAMOND for {len(aa_seqs)} sequence(s)..."
                )

        df = seq_workflow(
            sequences=aa_seqs,
            sequence_lengths=seq_lens,
            reference=reference,
            out=out,
            sensitivity=sensitivity,
            json=json,
            verbose=verbose,
        )

        if len(df) == 0:
            # TODO: change to warning
            logging.info(
                "No ELM database orthologs found for input sequence or UniProt ID."
            )

        if not uniprot and len(df) > 0:
            try:
                target_start_values = df["target_start"].values
                target_end_values = df["target_end"].values

                # # TO-DO: we do not want to drop these results, we just want to add a column "motif_in_overlap" True/False
                # if (df["Per. Ident"] is not None):
                #     # ignore nonoverlapping motifs
                #     df.drop(df[ (df['Start'] <= target_start[0]) | (df['End'] >= target_end[0]) ].index, inplace=True)

            except KeyError:
                logging.warning(
                    "No target start found for input sequence. If you entered a UniProt ID, please set 'uniprot' to True."
                )

    # building second data frame with regex motif match
    if uniprot:
        # use amino acid sequence associated with UniProt ID to do regex match
        df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, sequence)
        # sequences is an array
        sequences = df_uniprot[df_uniprot["uniprot_id"] == sequence]["sequence"].values

        # TODO What if no amino acid seqs are found for ID?
        if len(sequence) == 0:
            logging.warning(
                "No sequences found for UniProt ID from UniProt REST Server"
            )

        if len(sequences) > 1:
            logging.info(
                f"More than one UniProt amino acid sequence found for UniProt ID {sequence}. Using best match to find regex motifs."
            )
        sequence = sequences[0]

    df_regex_matches = regex_match(sequence)

    if len(df_regex_matches) == 0:
        logging.warning("No regex matches found for input sequence or UniProt ID.")

    if json:
        ortholog_dict = json_package.loads(df.to_json(orient="records"))
        regex_dict = json_package.loads(df_regex_matches.to_json(orient="records"))
        if out:
            with open("ortholog.json", "w", encoding="utf-8") as f:
                json_package.dump(ortholog_dict, f, ensure_ascii=False, indent=4)
            with open("regex.json", "w", encoding="utf-8") as f:
                json_package.dump(regex_dict, f, ensure_ascii=False, indent=4)

        return ortholog_dict, regex_dict

    else:
        ROOT_DIR = os.path.abspath(os.curdir)
        if out is None:
            # Create temporary results folder
            path = os.path.join(ROOT_DIR, "results")
        else:
            path = os.path.join(ROOT_DIR, out)
        try:
            if len(df) > 0 and len(df_regex_matches) > 0:
                df.to_csv(os.path.join(path, "ortholog"))
                df_regex_matches.to_csv(os.path.join(path, "regex_match"))
            elif len(df) > 0:
                df.to_csv(os.path.join(path, "ortholog"))
            elif len(df_regex_matches) > 0:
                df_regex_matches.to_csv(os.path.join(path, "regex_match"))

        except OSError:
            os.mkdir(path)

    return df, df_regex_matches
