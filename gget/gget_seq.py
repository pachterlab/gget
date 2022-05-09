import logging
# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
import numpy as np

# Custom functions
from .utils import rest_query, get_uniprot_seqs
from .gget_info import info

# Constants
from .constants import ENSEMBL_REST_API, UNIPROT_REST_API


def seq(
    ens_ids,
    seqtype="gene",
    isoforms=False,
    save=False,
):
    """
    Fetch nucleotide or amino acid sequence (FASTA) of a gene
    (and all its isoforms) or transcript by Ensembl ID.

    Args:
    - ens_ids   One or more Ensembl IDs (passed as string or list of strings).
    - seqtype   'gene' (default) or 'transcript'.
                Defines whether nucleotide or amino acid sequences are returned.
                Nucleotide sequences are fetched from the Ensembl REST API server.
                Amino acid sequences are fetched from the UniProt REST API server.
    - isoforms  If True, returns the sequences of all known transcripts (default: False).
                (Only for gene IDs in combination with seqtype='transcript'.)
    - save      If True: Save output FASTA to current directory.

    Returns a list (or FASTA file if 'save=True') containing the requested sequences.
    """

    ## Clean up arguments
    # Check if seqtype is valid
    seqtypes = ["gene", "transcript"]
    seqtype = seqtype.lower()
    if seqtype not in seqtypes:
        raise ValueError(
            f"Argument 'seqtype' was specified as '{seqtype}'. Expected one of {', '.join(seqtypes)}"
        )

    # Clean up Ensembl IDs
    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]
    # Remove Ensembl ID version if passed
    ens_ids_clean = []
    for ensembl_ID in ens_ids:
        ens_ids_clean.append(ensembl_ID.split(".")[0])

    # Initiate empty 'fasta'
    fasta = []

    ## Fetch nucleotide sequece
    if seqtype == "gene":
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
                        f"Ensembl ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )

            # If isoforms true, fetch sequences of isoforms instead
            if isoforms == True:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID, expand=True)

                # Check that Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.error(
                        f"Ensembl ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, get the IDs of all its transcripts
                if ens_ID_type == "Gene":
                    logging.info(
                        f"Requesting nucleotide sequences of all transcripts of {ensembl_ID} from Ensembl."
                    )

                    for transcipt_id in info_df.loc[ensembl_ID]["all_transcripts"]:

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
                                f"Ensembl ID {transcipt_id} not found. "
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
                            f"Ensembl ID {ensembl_ID} not found. "
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
    if seqtype == "transcript":
        if isoforms == False:
            # List to collect transcript IDs
            trans_ids = []

            for ensembl_ID in ens_ids_clean:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID)

                # Check that Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.error(
                        f"Ensembl ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, use the ID of its canonical transcript
                if ens_ID_type == "Gene":
                    # Get ID of canonical transcript
                    can_trans = info_df.loc[ensembl_ID]["canonical_transcript"]
                    # Remove Ensembl ID version from transcript IDs and append to transcript IDs list
                    trans_ids.append(can_trans.split(".")[0])
                    logging.info(
                        f"Requesting amino acid sequence of the canonical transcript {can_trans.split('.')[0]} of gene {ensembl_ID} from UniProt."
                    )

                # If the ID is a transcript, append the ID directly
                elif ens_ID_type == "Transcript":
                    trans_ids.append(ensembl_ID)
                    logging.info(
                        f"Requesting amino acid sequence of {ensembl_ID} from UniProt."
                    )

                else:
                    logging.warning(
                        f"{ensembl_ID} not recognized as either a gene or transcript ID. It will not be included in the UniProt query."
                    )

            # Fetch the amino acid sequences of the transcript Ensembl IDs
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, trans_ids)

        if isoforms == True:
            # List to collect transcript IDs
            trans_ids = []

            for ensembl_ID in ens_ids_clean:
                # Get ID type (gene, transcript, ...) using gget info
                info_df = info(ensembl_ID, expand=True)

                # Check that Ensembl ID was found
                if isinstance(info_df, type(None)):
                    logging.error(
                        f"Ensembl ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )
                    continue

                ens_ID_type = info_df.loc[ensembl_ID]["object_type"]

                # If the ID is a gene, get the IDs of all isoforms
                if ens_ID_type == "Gene":
                    # Get the IDs of all transcripts from the gget info results
                    for transcipt_id in info_df.loc[ensembl_ID]["all_transcripts"]:
                        # Append transcript ID (without Ensembl version number) to list of transcripts to fetch
                        trans_ids.append(transcipt_id.split(".")[0])

                    logging.info(
                        f"Requesting amino acid sequences of all transcripts of gene {ensembl_ID} from UniProt."
                    )

                elif ens_ID_type == "Transcript":
                    # Append transcript ID to list of transcripts to fetch
                    trans_ids.append(ensembl_ID)
                    logging.info(
                        f"Requesting amino acid sequence of {ensembl_ID} from UniProt."
                    )
                    logging.warning("The isoform option only applies to gene IDs.")
                else:
                    logging.warning(
                        f"{ensembl_ID} not recognized as either a gene or transcript ID. It will not be included in the UniProt query."
                    )

            # Fetch amino acid sequences of all isoforms from the UniProt REST API
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, trans_ids)

        # Check if less results were found than IDs put in
        if len(df_uniprot) != len(trans_ids) and len(df_uniprot) > 0:
            logging.warning(
                "The number of results does not match the number of IDs requested. "
                "It is possible that UniProt amino acid sequences were not found for some of the IDs."
            )

        # Check if no results were found
        if len(df_uniprot) < 1:
            logging.error(
                "No UniProt amino acid sequences were found for these Ensembl ID(s)."
            )

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
