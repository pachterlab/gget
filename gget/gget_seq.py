import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

import numpy as np

# Custom functions
from .utils import rest_query, get_uniprot_seqs
from .gget_info import info

# Constants
from .constants import ENSEMBL_REST_API, UNIPROT_REST_API


def seq(
    ens_ids,
    transcribe=False,
    seqtype=None,
    isoforms=False,
    save=False,
):
    """
    Fetch nucleotide or amino acid sequence (FASTA) of a gene
    (and all its isoforms) or transcript by Ensembl, WormBase or FlyBase ID.

    Args:
    - ens_ids       One or more Ensembl IDs (passed as string or list of strings).
                    Also supports WormBase and FlyBase IDs.
    - transcribe    True/False (default: False -> returns nucleotide sequences).
                    Defines whether nucleotide or amino acid sequences are returned.
                    Nucleotide sequences are fetched from the Ensembl REST API server.
                    Amino acid sequences are fetched from the UniProt REST API server.
    - isoforms      If True, returns the sequences of all known transcripts (default: False).
                    (Only for gene IDs.)
    - save          If True, saves output FASTA to current directory (default: False).

    Returns a list (or FASTA file if 'save=True') containing the requested sequences.

    Deprecated arguments: 'seqtype' (use True/False flag 'transcribe' instead.)
    """
    # Handle deprecated arguments
    if seqtype:
        logging.error(
            "'seqtype' argument deprecated! Please use True/False argument 'transcribe' instead."
        )
        return

    ## Clean up arguments
    # Clean up Ensembl IDs
    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]
    # Remove Ensembl ID version if passed
    ens_ids_clean = []
    temp = 0
    for ensembl_ID in ens_ids:
        # But only for Ensembl ID (and not for flybase/wormbase IDs)
        if ensembl_ID.startswith("ENS"):
            ens_ids_clean.append(ensembl_ID.split(".")[0])

            if "." in ensembl_ID and temp == 0:
                logging.info(
                    "We noticed that you may have passed a version number with your Ensembl ID.\n"
                    "Please note that gget seq will return information linked to the latest Ensembl ID version."
                )
                temp = +1

        else:
            ens_ids_clean.append(ensembl_ID)

    # Initiate empty 'fasta'
    fasta = []

    ## Fetch nucleotide sequece
    if transcribe is False:
        # Define Ensembl REST API server
        server = ENSEMBL_REST_API
        # Define type of returned content from REST
        content_type = "application/json"

        # Initiate dictionary to save results for all IDs in
        master_dict = {}

        # Query REST APIs from https://rest.ensembl.org/
        for ensembl_ID in ens_ids_clean:
            # Create dict to save query results
            results_dict = {ensembl_ID: {}}

            # If isoforms False, just fetch sequences of passed Ensembl ID
            if isoforms == False:
                # sequence/id/ query: Request sequence by stable identifier
                query = "sequence/id/" + ensembl_ID + "?"

                # Try if query valid
                try:
                    # Submit query; this will throw RuntimeError if ID not found
                    df_temp = rest_query(server, query, content_type)

                    # Delete superfluous entries
                    keys_to_delete = ["query", "id", "version", "molecule"]
                    for key in keys_to_delete:
                        # Pop keys, None -> do not raise an error if key to delete not found
                        df_temp.pop(key, None)

                    # Add results to main dict
                    results_dict[ensembl_ID].update({"seq": df_temp})

                    logging.info(
                        f"Requesting nucleotide sequence of {ensembl_ID} from Ensembl."
                    )

                except RuntimeError:
                    logging.error(
                        f"ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )

            # If isoforms true, fetch sequences of isoforms instead
            if isoforms == True:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID, verbose=False)

                # Check if Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.warning(
                        f"ID '{ensembl_ID}' not found. Please double-check spelling/arguments."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, get the IDs of all its transcripts
                if ens_ID_type == "Gene":
                    logging.info(
                        f"Requesting nucleotide sequences of all transcripts of {ensembl_ID} from Ensembl."
                    )

                    for transcipt_id in info_df.loc[ensembl_ID]["all_transcripts"]:
                        # Remove version number for Ensembl IDs (not for flybase/wormbase IDs)
                        if transcipt_id.startswith("ENS"):
                            transcipt_id = transcipt_id.split(".")[0]

                        # Try if query is valid
                        try:
                            # Define the REST query
                            query = "sequence/id/" + transcipt_id + "?"
                            # Submit query
                            df_temp = rest_query(server, query, content_type)

                            # Delete superfluous entries
                            keys_to_delete = ["query", "version", "molecule"]
                            for key in keys_to_delete:
                                # Pop keys, None -> do not raise an error if key to delete not found
                                df_temp.pop(key, None)

                            # Add results to main dict
                            results_dict[ensembl_ID].update(
                                {f"{transcipt_id}": df_temp}
                            )

                        except RuntimeError:
                            logging.error(
                                f"ID {transcipt_id} not found. "
                                "Please double-check spelling/arguments and try again."
                            )

                # If isoform true, but ID is not a gene; ignore the isoform parameter
                else:
                    # Try if query is valid
                    try:
                        # Define the REST query
                        query = "sequence/id/" + ensembl_ID + "?"

                        # Submit query
                        df_temp = rest_query(server, query, content_type)

                        # Delete superfluous entries
                        keys_to_delete = ["query", "id", "version", "molecule"]
                        for key in keys_to_delete:
                            # Pop keys, None -> do not raise an error if key to delete not found
                            df_temp.pop(key, None)

                        # Add results to main dict
                        results_dict[ensembl_ID].update({"seq": df_temp})

                        logging.info(
                            f"Requesting nucleotide sequence of {ensembl_ID} from Ensembl."
                        )
                        logging.warning("The isoform option only applies to gene IDs.")

                    except RuntimeError:
                        logging.error(
                            f"ID {ensembl_ID} not found. "
                            "Please double-check spelling/arguments and try again."
                        )

            # Add results to master dict
            master_dict.update(results_dict)

        # Build FASTA file
        for ens_ID in master_dict:
            for key in master_dict[ens_ID].keys():
                if key == "seq":
                    fasta.append(">" + ens_ID + " " + master_dict[ens_ID][key]["desc"])
                    fasta.append(master_dict[ens_ID][key]["seq"])
                else:
                    fasta.append(
                        ">"
                        + master_dict[ens_ID][key]["id"]
                        + " "
                        + master_dict[ens_ID][key]["desc"]
                    )
                    fasta.append(master_dict[ens_ID][key]["seq"])

    ## Fetch amino acid sequences from UniProt
    if transcribe:
        if isoforms is False:
            # List to collect transcript IDs
            trans_ids = []

            for ensembl_ID in ens_ids_clean:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID, verbose=False)

                # Check that Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.warning(
                        f"ID '{ensembl_ID}' not found. Please double-check spelling/arguments."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, use the ID of its canonical transcript
                if ens_ID_type == "Gene":
                    # Get ID of canonical transcript
                    can_trans = info_df.loc[ensembl_ID]["canonical_transcript"]

                    if ensembl_ID.startswith("ENS"):
                        # Remove Ensembl ID version from transcript IDs and append to transcript IDs list
                        temp_trans_id = can_trans.split(".")[0]
                        trans_ids.append(temp_trans_id)

                    elif ensembl_ID.startswith("WB"):
                        # Remove added "." at the end of transcript IDs
                        temp_trans_id1 = ".".join(can_trans.split(".")[:-1])
                        # For WormBase transcript IDs, also remove the version number for submission to UniProt API
                        temp_trans_id = ".".join(temp_trans_id1.split(".")[:-1])
                        trans_ids.append(temp_trans_id)

                    else:
                        # Remove added "." at the end of other transcript IDs
                        temp_trans_id = ".".join(can_trans.split(".")[:-1])
                        trans_ids.append(temp_trans_id)

                    logging.info(
                        f"Requesting amino acid sequence of the canonical transcript {temp_trans_id} of gene {ensembl_ID} from UniProt."
                    )

                # If the ID is a transcript, append the ID directly
                elif ens_ID_type == "Transcript":
                    # For WormBase transcript IDs, remove the version number for submission to UniProt API
                    if ensembl_ID.startswith("T"):
                        trans_ids.append(".".join(ensembl_ID.split(".")[:-1]))
                    else:
                        trans_ids.append(ensembl_ID)

                    logging.info(
                        f"Requesting amino acid sequence of {ensembl_ID} from UniProt."
                    )

                else:
                    logging.warning(
                        f"{ensembl_ID} not recognized as either a gene or transcript ID. It will not be included in the UniProt query."
                    )

            # Check if this is a Wrombase ID:
            if ensembl_ID.startswith("WB") or ensembl_ID.startswith("T"):
                id_type = "wormbase"
            # Check if this is a flybase ID:
            elif ensembl_ID.startswith("FB"):
                id_type = "flybase"
            else:
                id_type = "ensembl"

            # Fetch the amino acid sequences of the transcript Ensembl IDs
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, trans_ids, id_type=id_type)

        if isoforms == True:
            # List to collect transcript IDs
            trans_ids = []

            for ensembl_ID in ens_ids_clean:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID, verbose=False)

                # Check that Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.warning(
                        f"ID '{ensembl_ID}' not found. Please double-check spelling/arguments."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, get the IDs of all isoforms
                if ens_ID_type == "Gene":
                    # Get the IDs of all transcripts from the gget info results
                    for transcipt_id in info_df.loc[ensembl_ID]["all_transcripts"]:
                        if ensembl_ID.startswith("ENS"):
                            # Append transcript ID (without Ensembl version number) to list of transcripts to fetch
                            trans_ids.append(transcipt_id.split(".")[0])

                        elif ensembl_ID.startswith("WB"):
                            # For WormBase transcript IDs, remove the version number for submission to UniProt API
                            temp_trans_id = ".".join(transcipt_id.split(".")[:-1])
                            trans_ids.append(temp_trans_id)

                        else:
                            # Note: No need to remove the added "." at the end of unversioned transcripts here, because "all_transcripts" are returned without it
                            trans_ids.append(transcipt_id)

                    logging.info(
                        f"Requesting amino acid sequences of all transcripts of gene {ensembl_ID} from UniProt."
                    )

                elif ens_ID_type == "Transcript":
                    # For WormBase transcript IDs, remove the version number for submission to UniProt API
                    if ensembl_ID.startswith("T"):
                        trans_ids.append(".".join(ensembl_ID.split(".")[:-1]))

                    else:
                        trans_ids.append(ensembl_ID)

                    logging.info(
                        f"Requesting amino acid sequence of {ensembl_ID} from UniProt."
                    )
                    logging.warning("The isoform option only applies to gene IDs.")

                else:
                    logging.warning(
                        f"{ensembl_ID} not recognized as either a gene or transcript ID. It will not be included in the UniProt query."
                    )

            # Check if this is a Wormbase ID:
            if ensembl_ID.startswith("WB") or ensembl_ID.startswith("T"):
                id_type = "wormbase"
            # Check if this is a flybase ID:
            elif ensembl_ID.startswith("FB"):
                id_type = "flybase"
            else:
                id_type = "ensembl"

            # Fetch amino acid sequences of all isoforms from the UniProt REST API
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, trans_ids, id_type=id_type)

        # Check if less results were found than IDs put in
        if len(df_uniprot) != len(trans_ids) and len(df_uniprot) > 0:
            logging.warning(
                "The number of results does not match the number of IDs requested. "
                "It is possible that UniProt amino acid sequences were not found for some of the IDs."
            )

        # Check if no results were found
        if len(df_uniprot) < 1:
            logging.error("No UniProt amino acid sequences were found for these ID(s).")

        else:
            # Build UniProt results FASTA file
            for (
                uniprot_id,
                query_ensembl_id,
                gene_name,
                organism,
                sequence_length,
                uniprot_seq,
            ) in zip(
                df_uniprot["uniprot_id"].values,
                df_uniprot["query"].values,
                df_uniprot["gene_name"].values,
                df_uniprot["organism"].values,
                df_uniprot["sequence_length"].values,
                df_uniprot["sequence"].values,
            ):
                fasta.append(
                    ">"
                    + query_ensembl_id
                    + " uniprot_id: "
                    + uniprot_id
                    + " ensembl_id: "
                    + query_ensembl_id
                    + " gene_name(s): "
                    + gene_name
                    + " organism: "
                    + organism
                    + " sequence_length: "
                    + str(sequence_length)
                )
                fasta.append(uniprot_seq)

    # Save
    if save:
        file = open("gget_seq_results.fa", "w")
        for element in fasta:
            file.write(element + "\n")
        file.close()

    return fasta
