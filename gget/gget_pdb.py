from urllib.request import urlopen
from urllib.error import HTTPError
import json
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

from .constants import RCSB_PDB_API


def pdb(pdb_id, resource="pdb", identifier=None, save=False):
    """
    Query RCSB PDB for the protein structutre/metadata of a given PDB ID.

    Args:
    - pdb_id        PDB ID to be queried (str), e.g. "7S7U".
    - resource      Defines type of information to be returned.
                    "pdb": Returns the protein structure in PDB format (default).
                    "entry": Information about PDB structures at the top level of PDB structure hierarchical data organization.
                    "pubmed": Get PubMed annotations (data integrated from PubMed) for a given entry's primary citation.
                    "assembly": Information about PDB structures at the quaternary structure level.
                    "branched_entity": Get branched entity description (define entity ID as "identifier").
                    "nonpolymer_entity": Get non-polymer entity data (define entity ID as "identifier").
                    "polymer_entity": Get polymer entity data (define entity ID as "identifier").
                    "uniprot": Get UniProt annotations for a given macromolecular entity (define entity ID as "identifier").
                    "branched_entity_instance": Get branched entity instance description (define chain ID as "identifier").
                    "polymer_entity_instance": Get polymer entity instance (a.k.a chain) data (define chain ID as "identifier").
                    "nonpolymer_entity_instance": Get non-polymer entity instance description (define chain ID as "identifier").
    -  identifier   Can be used to define assembly, entity or chain ID if applicable (default: None).
                    Assembly/entity IDs are numbers (e.g. 1), and chain IDs are letters (e.g. "A").
    - save          True/False wether to save JSON/PDB with query results in the current working directory (default: False).

    Returns requested information in JSON format (except for resource="pdb" which returns protein structure in PDB format).
    """

    # Check if resource argument is valid
    resources = [
        "pdb",
        "entry",
        "pubmed",
        "assembly",
        "branched_entity",
        "nonpolymer_entity",
        "polymer_entity",
        "uniprot",
        "branched_entity_instance",
        "polymer_entity_instance",
        "nonpolymer_entity_instance",
    ]
    if resource not in resources:
        raise ValueError(
            f"'resource' argument specified as {resource}. Expected one of: {', '.join(resources)}"
        )

    # Check if required identifiers are present
    if resource == "assembly" and identifier is None:
        raise ValueError("Please define assembly ID (e.g. '1') as 'identifier'.")

    need_entitiy_id = [
        "branched_entity",
        "nonpolymer_entity",
        "polymer_entity",
        "uniprot",
    ]
    if resource in need_entitiy_id and identifier is None:
        raise ValueError("Please define entity ID (e.g. '1') as 'identifier'.")

    need_chain_id = [
        "branched_entity_instance",
        "nonpolymer_entity_instance",
        "polymer_entity_instance",
    ]
    if resource in need_chain_id and identifier is None:
        raise ValueError("Please define chain ID (e.g. 'A') as 'identifier'.")

    # Define URLs for HTTP request
    if resource != "pdb":
        # URLs to request resources other than PDB file
        if identifier is not None:
            url = f"{RCSB_PDB_API}{resource}/{pdb_id}/{identifier}"
        else:
            url = f"{RCSB_PDB_API}{resource}/{pdb_id}"

    else:
        # URL to request PDB file
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    # Submit URL request
    try:
        r = urlopen(url)
    except HTTPError:
        if resource == "assembly":
            logging.error(
                f"{resource} for {pdb_id} assembly {identifier} was not found. Please double-check arguments and try again."
            )
        elif resource in need_entitiy_id:
            logging.error(
                f"{resource} for {pdb_id} entity {identifier} was not found. Please double-check arguments and try again."
            )
        elif resource in need_chain_id:
            logging.error(
                f"{resource} for {pdb_id} chain {identifier} was not found. Please double-check arguments and try again."
            )
        else:
            logging.error(
                f"{resource} for {pdb_id} was not found. Please double-check arguments and try again."
            )
        return

    if r.status != 200:
        raise RuntimeError(
            f"The RCSB PDB server responded with status code: {r.status}. "
            "Please double-check arguments and try again.\n"
        )

    if resource != "pdb":
        # Read json formatted results
        results = json.load(r)
    else:
        # Read PDB file
        results = r.read().decode()

    if save:
        if resource != "pdb":
            # Save the results in json format
            if identifier is not None:
                out_name = f"{pdb_id}_{identifier}_{resource}.json"
            else:
                out_name = f"{pdb_id}_{resource}.json"

            with open(out_name, "w", encoding="utf-8") as f:
                json.dump(results, f, ensure_ascii=False, indent=4)

        else:
            # Save the PDB file
            with open(f"{pdb_id}.pdb", "w") as f:
                f.write(results)

    return results
