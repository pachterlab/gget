# Copyright 2022 Laura Luebbert

from bs4 import BeautifulSoup
import requests
import re
import numpy as np

def rest_query(server, query, content_type):
    """
    Function to query a 

    Parameters:
    - server
    Serve to query.
    - Query
    Query that is passed to server.
    - content_type
    Contect type requested from server.

    Returns server output.
    """

    r = requests.get(
        server + query, 
        headers={ "Content-Type" : content_type}
    )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text
    
def gget_species_options(release=105):
    """
    Function to find all available species core databases for gget.

    Parameters:
    - release
    Ensembl release for which the databases are fetched.

    Returns list of available core databases.
    """
    # Find all available databases
    url = f"http://ftp.ensembl.org/pub/release-{release}/mysql/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
        
    soup = BeautifulSoup(html.text, "html.parser")

    # Return list of all available databases
    databases = []
    for subsoup in soup.body.findAll('a'):
        if "core" in subsoup["href"]:
            databases.append(subsoup["href"].split("/")[0])
    
    return databases
    
    
def ref_species_options(which, release=None):
    """
    Function to find all available species for gget ref.

    Parameters:
    - release
    Ensembl release for which available species should be fetched.
    - which
    Which type of FTP. Possible entries: 'dna', 'cdna', 'gtf'.

    Returns list of available species.
    """
    ## Find latest Ensembl release
    url = "http://ftp.ensembl.org/pub/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
        
    soup = BeautifulSoup(html.text, "html.parser")
    # Find all releases
    releases = soup.body.findAll(text=re.compile("release-"))
    # Get release numbers
    rels = []
    for rel in releases:
        rels.append(rel.split("/")[0].split("-")[-1])

    # Find highest release number (= latest release)
    ENS_rel = np.array(rels).astype(int).max()
        
    # If release != None, use user-defined Ensembl release    
    if release != None:
        # Do not allow user-defined release if it is higher than the latest release
        if release > ENS_rel:
            raise ValueError("Defined Ensembl release number cannot be greater than latest release.")
        else:
            ENS_rel = release

    # Find all available species for this release and FTP type
    if which == "gtf":
        url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/gtf/"
    if which == "dna" or which == "cdna":
        url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")

    sps = []
    for subsoup in soup.body.findAll('a'):
        sps.append(subsoup["href"].split("/")[0])

    species_list = sps[1:]
    
    # Return list of all available species
    return species_list


def parse_blast_ref_page(handle):
    """
    Extract a tuple of RID, RTOE from the NCBI 'please wait' page.
    RTOE = 'Estimated time fo completion.'
    RID = 'Request ID'.

    Code adapted from the Biopython BLAST NCBIWWW project written
    by Jeffrey Chang (Copyright 1999), Brad Chapman, and Chris Wroe distributed under the
    Biopython License Agreement and BSD 3-Clause License
    https://github.com/biopython/biopython/blob/171697883aca6894f8367f8f20f1463ce7784d0c/LICENSE.rst
    """
    s = handle.read().decode()
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i + len("RID =") : j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i + len("RTOE =") : j].strip()

    if not rid and not rtoe:
        # Can we reliably extract the error message from the HTML page?
        # e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        # or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        # This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i + len('<div class="error msInf">') :].strip()
            msg = msg.split("</div>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i + len('<p class="error">') :].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # Generic search based on the way the error messages start:
        i = s.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = s[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError("Error message from NCBI: %s" % msg)
        # If we cannot recognise the error layout:
        raise ValueError(
            "No request ID and no estimated time to completion found in the NCBI 'please wait' page, "
            "there was probably an error in your request but we "
            "could not extract a helpful error message."
        )
    elif not rid:
        # Can this happen?
        raise ValueError(
            "No request ID found in the 'please wait' page. (Although estimated time to completion = %r)"
            % rtoe
        )
    elif not rtoe:
        # Can this happen?
        raise ValueError(
            "No estimated time to completion found in the 'please wait' page. (Although request ID = %r)"
            % rid
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            "A non-integer estimated time to completion found in the 'please wait' page, %r"
            % rtoe
        ) from None
