import pandas as pd
import numpy as np
import os
import logging
import json as json_package
import re

from .utils import get_uniprot_seqs, tsv_to_df
from .constants import UNIPROT_REST_API
from .gget_diamond import diamond
from .gget_setup import (
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
        if (row["motif_start_in_target"] >= row["target_start"])
        & (row["motif_end_in_target"] <= row["target_end"])
        else False
    )


def get_elm_instances(UniProtID):
    """
    Get ELM instances and their information from local ELM tsv files.

    Args:
    - UniProtID   UniProt ID to search for in the accession column of ELM tsv files.

    Returns: dataframe combining ELM instances and information (description, functional site...)
    """
    # Get matching rows from elm_instances.tsv
    # ELM Instances.tsv file contains 5 lines before headers and data
    df_full_instances = tsv_to_df(ELM_INSTANCES_TSV, skiprows=5)
    df_instances_matching = df_full_instances[
        df_full_instances["Primary_Acc"] == UniProtID
    ]
    # Rename columns
    df_instances_matching = df_instances_matching.rename(
        columns={
            "Primary_Acc": "Ortholog_UniProt_ID",
            "Start": "motif_start_in_target",
            "End": "motif_end_in_target",
        }
    )

    # Get matching class descriptions from elm_classes.tsv
    df_classes = tsv_to_df(ELM_CLASSES_TSV, skiprows=5)
    df_classes.rename(columns={"Accession": "class_accession"}, inplace=True)

    # Merge dataframes using ELM Identifier
    df_final = df_instances_matching.merge(df_classes, how="left", on="ELMIdentifier")

    return df_final


def seq_workflow(
    sequences,
    sequence_lengths,
    reference,
    sensitivity,
    threads,
    verbose,
    diamond_binary,
):
    """
    Alignment of sequence using DIAMOND to get UniProt ID. Use the UniProt ID to construct an ortholog dataframe similar to the UniProt workflow
    except for additional columns for start, end and whether the motif overlaps the target sequence.

    Args:
    sequences        - list of user input amino acid sequence
    sequence_lengths - list of lengths respective to each sequence
    reference        - Path to reference FASTA file
    sensitivity      - Sensitivity of DIAMOND alignment.
                       One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive or ultra-sensitive.
    threads          - Number of threads used for DIAMOND alignment
    verbose          - If True, turns on logging for INFO_level messages
    diamond_binary   - Path to DIAMOND binary

    Returns: data frame consisting of ELM instances, class information, start, end in query, and if motif overlaps with target sequence
    """
    if verbose:
        logging.info(
            f"Performing pairwise sequence alignment against ELM database using DIAMOND for {len(sequences)} sequence(s)..."
        )

    df = pd.DataFrame()
    seq_number = 1
    for sequence, seq_len in zip(sequences, sequence_lengths):
        df_diamond = diamond(
            query=sequence,
            reference=reference,
            sensitivity=sensitivity,
            threads=threads,
            verbose=verbose,
            diamond_binary=diamond_binary,
        )

        if len(df_diamond) == 0:
            logging.warning(
                f"Sequence #{seq_number}: No orthologous proteins found in ELM database."
            )

        else:
            # Construct df with elm instances from UniProt ID returned from diamond
            # TODO double check that this gets info if more than one UniProt ID matched
            uniprot_ids = str(df_diamond["target_accession"]).split("|")[1]
            if verbose:
                logging.info(
                    f"Sequence #{seq_number}: DIAMOND found the following orthologous proteins: {', '.join(uniprot_ids)}. Retrieving ELMs for each UniProt ID..."
                )

            for i, uniprot_id in enumerate(df_diamond["target_accession"].values):
                # print(f"UniProt ID {uniprot_id}")
                df_elm = get_elm_instances(str(uniprot_id).split("|")[1])
                # missing motifs other than the first one
                df_elm["Query Cover"] = df_diamond["length"].values[i] / seq_len * 100
                df_elm["Per. Ident"] = df_diamond["Per. Ident"].values[i]
                df_elm["query_start"] = int(df_diamond["query_start"].values[i])
                df_elm["query_end"] = int(df_diamond["query_end"].values[i])
                df_elm["target_start"] = int(df_diamond["target_start"].values[i])
                df_elm["target_end"] = int(df_diamond["target_end"].values[i])
                df_elm["Motif_in_query"] = df_elm.apply(motif_in_query, axis=1)

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
    df_elm_classes = tsv_to_df(ELM_CLASSES_TSV, skiprows=5)
    df_full_instances = tsv_to_df(ELM_INSTANCES_TSV, skiprows=5)

    elm_ids = df_elm_classes["Accession"]
    regex_patterns = df_elm_classes["Regex"]

    df_final = pd.DataFrame()

    # Compare ELM regex with input sequence and return all matching elms
    for elm_id, pattern in zip(elm_ids, regex_patterns):
        regex_matches = re.finditer(pattern, sequence)

        for match_string in regex_matches:
            elm_row = df_elm_classes[df_elm_classes["Accession"] == elm_id]
            elm_row.insert(
                loc=1,
                column="Instances (Matched Sequence)",
                value=match_string.group(0),
            )

            (start, end) = match_string.span()
            elm_row.insert(loc=2, column="motif_start_in_query", value=str(start))
            elm_row.insert(loc=3, column="motif_end_in_query", value=str(end))

            elm_identifier = [str(x) for x in elm_row["ELMIdentifier"]][0]

            # df_instances_matching = df_full_instances.loc[
            #     df_full_instances["ELMIdentifier"] == elm_identifier
            # ]

            # merge two dataframes using ELM Identifier, since some Accessions are missing from elm_instances.tsv
            elm_row = elm_row.merge(
                df_full_instances, how="left", on="ELMIdentifier"
            )

            df_final = pd.concat([df_final, elm_row])

    # df_final.pop("Accession_y")
    # df_final.pop("#Instances")
    # df_final.pop("#Instances_in_PDB")
    # df_final.pop("References")
    # df_final.pop("InstanceLogic")

    if len(df_final) > 0:
        df_final.rename(columns={"Accession_x": "Instance_accession"}, inplace=True)

        # Reorder columns
        change_column = [
            "Instance_accession",
            "ELMIdentifier",
            "FunctionalSiteName",
            "ELMType",
            "Description",
            "Regex",
            "Instances (Matched Sequence)",
            "Probability",
            "motif_start_in_query",
            "motif_end_in_query",
            "Methods",
            "ProteinName",
            "Organism",
            "References",
            "InstanceLogic",
            "#Instances",
            "#Instances_in_PDB",
        ]

        for col in change_column:
            if col not in df_final.columns:
                df_final[col] = np.NaN

        df_final = df_final[change_column]

    return df_final


def elm(
    sequence,
    uniprot=False,
    sensitivity="very-sensitive",
    threads=1,
    diamond_binary=None,
    verbose=True,
    json=False,
    out=None,
):
    """
    Searches the Eukaryotic Linear Motif resource for Functional Sites in Proteins.

    Args:
    - sequence         Amino acid sequence or Uniprot ID (str).
                       If Uniprot ID, set 'uniprot==True'.
    - uniprot          Set to True if input is a Uniprot ID instead of amino acid sequence. Default: False.
    - sensitivity      Sensitivity of DIAMOND alignment.
                       One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive or ultra-sensitive.
                       Default: "very-sensitive"
    - threads          Number of threads used in DIAMOND alignment. Default: 1.
    - diamond_binary   Path to DIAMOND binary. Default: None -> Uses DIAMOND binary installed with gget.
    - verbose          True/False whether to print progress information. Default: True.
    - json             If True, returns results in json format instead of data frame. Default: False.
    - out              Path to folder to save DIAMOND results in. Default: Standard out, temporary files are deleted.

    Returns two data frames: orthologs and regex matches from ELM database.
    """
    # Check if ELM files were downloaded
    if (
        not os.path.exists(ELM_INSTANCES_FASTA)
        or not os.path.exists(ELM_CLASSES_TSV)
        or not os.path.exists(ELM_INSTANCES_TSV)
    ):
        raise FileNotFoundError(
            f"Some or all ELM database files are missing. Please run 'gget setup elm' (Python: gget.setup('elm')) once to download the necessary files."
        )

    # Let user know when local ELM was last updated
    lines_number = 2
    with open(ELM_CLASSES_TSV) as input_file:
        head = [next(input_file) for _ in range(lines_number)]
    if verbose:
        logging.info(head)
    with open(ELM_INSTANCES_TSV) as input_file:
        head = [next(input_file) for _ in range(lines_number)]
    if verbose:
        logging.info(head)

    if not uniprot:
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZXBJZ")
        # Convert input sequence to upper case letters
        sequence = sequence.upper()

        # If sequence is not a valid amino sequence, raise error
        if not set(sequence) <= amino_acids:
            logging.warning(
                f"Input amino acid sequence contains invalid characters. If the input is a UniProt ID, please use flag --uniprot (Python: uniprot=True)."
            )

    # Build ortholog dataframe
    ortho_df = pd.DataFrame()
    if uniprot:
        ortho_df = get_elm_instances(sequence)

        if len(ortho_df) == 0:
            logging.warning(
                "UniProt ID does not match UniProt IDs in the ELM database. Fetching amino acid sequence from UniProt..."
            )
            df_uniprot = get_uniprot_seqs(server=UNIPROT_REST_API, ensembl_ids=sequence)

            try:
                # Only grab sequences where IDs match exactly
                aa_seqs = df_uniprot[df_uniprot["uniprot_id"] == sequence][
                    "sequence"
                ].values
                seq_lens = df_uniprot["sequence_length"].values

            except KeyError:
                raise ValueError(
                    f"No sequences found for UniProt ID {sequence} from the UniProt server. Please double check your UniProt ID and try again."
                )

    if len(ortho_df) == 0:
        # Add input aa sequence and its length to list
        if not uniprot:
            aa_seqs = [sequence]
            seq_lens = [len(sequence)]

        ortho_df = seq_workflow(
            sequences=aa_seqs,
            sequence_lengths=seq_lens,
            reference=ELM_INSTANCES_FASTA,
            sensitivity=sensitivity,
            threads=threads,
            verbose=verbose,
            diamond_binary=diamond_binary,
        )

        if len(ortho_df) == 0:
            logging.warning(
                "No ELM database orthologs found for input sequence or UniProt ID."
            )

        # if not uniprot and len(ortho_df) > 0:
        #     try:
        #         target_start_values = ortho_df["target_start"].values
        #         target_end_values = ortho_df["target_end"].values

        #         # # TO-DO: we do not want to drop these results, we just want to add a column "motif_in_overlap" True/False
        #         # if (ortho_df["Per. Ident"] is not None):
        #         #     # ignore nonoverlapping motifs
        #         #     ortho_df.drop(ortho_df[ (ortho_df['Start'] <= target_start[0]) | (ortho_df['End'] >= target_end[0]) ].index, inplace=True)

        #     except KeyError:
        #         logging.warning(
        #             "No target start found for input sequence. If you entered a UniProt ID, please set 'uniprot' to True."
        #         )

    # Reorder columns of ortholog data frame
    final_cols = [
        "Ortholog_UniProt_ID",
        "ProteinName",
        "class_accession",
        "ELMIdentifier",
        "FunctionalSiteName",
        "Description",
        "Regex",
        "Probability",
        "Query Cover",
        "Per. Ident",
        "Motif_in_query",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "motif_start_in_target",
        "motif_end_in_target",
        "Organism",
        "References",
        "InstanceLogic",
        "PDB",
        "#Instances",
        "#Instances_in_PDB",
    ]
    for col in final_cols:
        if col not in ortho_df.columns:
            ortho_df[col] = np.NaN

    ortho_df = ortho_df[final_cols]

    # Build data frame containing regex motif matches
    fetch_aa_failed = False
    if uniprot:
        # use amino acid sequence associated with UniProt ID to do regex match

        # do not fetch sequence again if already done above
        if not "df_uniprot" in locals():
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, sequence)

        try:
            # Only grab sequences where IDs match exactly
            sequences = df_uniprot[df_uniprot["uniprot_id"] == sequence][
                "sequence"
            ].values

            if len(sequences) > 1:
                logging.warning(
                    f"More than one amino acid sequence found for UniProt ID {sequence}. Using best match to find regex motifs."
                )

            sequence = sequences[0]

        except KeyError:
            logging.warning(
                "No sequences found for UniProt ID {sequence} from the UniProt server."
            )
            fetch_aa_failed = True

    df_regex_matches = pd.DataFrame()
    if not fetch_aa_failed:
        df_regex_matches = regex_match(sequence)

    if len(df_regex_matches) == 0:
        logging.warning("No regex matches found for input sequence or UniProt ID.")

    # Create out folder if it does not exist
    if out:
        directory = "/".join(out.split("/")[:-1])
        if directory != "":
            os.makedirs(directory, exist_ok=True)

    if json:
        ortholog_dict = json_package.loads(ortho_df.to_json(orient="records"))
        regex_dict = json_package.loads(df_regex_matches.to_json(orient="records"))

        if out:
            with open(f"{out}/ELM_ortho_results.json", "w", encoding="utf-8") as f:
                json_package.dump(ortholog_dict, f, ensure_ascii=False, indent=4)
            with open(f"{out}/ELM_regex_results.json", "w", encoding="utf-8") as f:
                json_package.dump(regex_dict, f, ensure_ascii=False, indent=4)

        return ortholog_dict, regex_dict

    else:
        if out:
            if len(ortho_df) > 0 and len(df_regex_matches) > 0:
                ortho_df.to_csv(f"{out}/ELM_ortho_results.csv", index=False)
                df_regex_matches.to_csv(f"{out}/ELM_regex_results.csv", index=False)

            elif len(ortho_df) > 0:
                ortho_df.to_csv(f"{out}/ELM_ortho_results.csv", index=False)

            elif len(df_regex_matches) > 0:
                df_regex_matches.to_csv(f"{out}/ELM_regex_results.csv", index=False)

    return ortho_df, df_regex_matches
