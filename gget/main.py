# Packages for use from terminal
import __init__ 
import argparse
import sys
import os

# Import gget code
from gget import *
import help_
import utils

def main():
    """
    Function containing argparse parsers and arguments to allow the use of gget from the terminal.
    """
    # Define parent parser 
    parent_parser = argparse.ArgumentParser(description=f"gget v{__init__.__version__}", add_help=False)
    # Initiate subparsers
    parent_subparsers = parent_parser.add_subparsers(dest="command")
    # Define parent (not sure why I need both parent parser and parent, but otherwise it does not work)
    parent = argparse.ArgumentParser(add_help=False)
    
    # Add custom help argument to parent parser
    parent_parser.add_argument(
            "-h","--help",
            action="store_true",
            help="Print manual. Recommendation: pipe into less by running 'gget -h | less'"
    )
    # Add custom version argument to parent parser
    parent_parser.add_argument(
            "-v","--version",
            action="store_true",
            help="Print version."
    )
    
    ## gget ref subparser
    parser_ref = parent_subparsers.add_parser(
        "ref",
        parents=[parent],
        description="Fetch FTP links for a specific species from Ensemble.",
        help="Fetch FTP links for a specific species from Ensemble.",
        add_help=False
        )
    # ref parser arguments
    parser_ref.add_argument(
        "--species", "-s", 
        default=None,
        type=str,
        metavar="/ -s",
        help="Species for which the FTPs will be fetched, e.g. homo_sapiens."
    )
    # ref parser arguments
    parser_ref.add_argument(
        "--list", "-l", 
        default=None, 
        action="store_true",
        required=False,
        help="List out all available species."
    )
    parser_ref.add_argument(
        "--which", "-w", 
        default="all", 
        type=str,
        nargs='+',
        required=False,
        metavar=(""),
        help=("Defines which results to return." 
              "Possible entries are:"
              "'all' - Returns GTF, cDNA, and DNA links and associated info (default)." 
              "Or one or a combination of the following:"  
              "'gtf' - Returns the GTF FTP link and associated info." 
              "'cdna' - Returns the cDNA FTP link and associated info."
              "'dna' - Returns the DNA FTP link and associated info."
             )
        )
    parser_ref.add_argument(
        "--release", "-r",
        default=None,  
        type=int, 
        required=False,
        metavar="/ -r",
        help="Ensemble release the FTPs will be fetched from, e.g. 104 (default: latest Ensembl release).")
    parser_ref.add_argument(
        "-ftp", "--ftp",  
        default=False, 
        action='store_true',
        required=False,
        help="If True: return only the FTP link instead of a json.")
    parser_ref.add_argument(
        "--out", "-o",
        type=str,
        required=False,
        metavar=("/ -o"),
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json." 
            "Default: None (just prints results)."
        )
    )

    ## gget search subparser
    parser_gget = parent_subparsers.add_parser(
        "search",
         parents=[parent],
         description="Query Ensembl for genes based on species and free form search terms.", 
         help="Query Ensembl for genes based on species and free form search terms.",
         add_help=False
         )
    # Search parser arguments
    parser_gget.add_argument(
        "--searchwords", "-sw", 
        type=str, 
        nargs="+",
        required=True, 
        metavar="",  
        help="One or more free form searchwords for the query, e.g. gaba, nmda."
    )
    parser_gget.add_argument(
        "--species", "-s",
        type=str,  
        required=True, 
        metavar="/ -s",
        help="Species to be queried, e.g. homo_sapiens."
    )
    parser_gget.add_argument(
        "--operator", "-op",
        type=str,
        choices=["or", "and"],
        default="or",
        required=False, 
        metavar="/ -op",
        help=("'or': Returns all genes that include at least one of the searchwords."
              "'and': Only returns genes that include all searchwords.")
    )
    parser_gget.add_argument(
        "--limit", "-l",
        type=int, 
        default=None,
        required=False,
        metavar="/ -l",
        help="Limits the number of results, e.g. 10 (default: None)."
    )
    parser_gget.add_argument(
        "--out", "-o",
        type=str,
        required=False,
        metavar="/ -o",
        help=(
            "Path to the csv file the results will be saved in, e.g. path/to/directory/results.csv." 
            "Default: None (just prints results)."
        )
    )
    
    ## gget info subparser
    parser_info = parent_subparsers.add_parser(
        "info",
        parents=[parent],
        description="Look up information about Ensembl IDs.", 
        help="Look up information about Ensembl IDs.",
        add_help=False
        )
    # info parser arguments
    parser_info.add_argument(
        "--ens_ids", "-id", 
        type=str,
        nargs="+",
        required=True, 
        metavar="",   
        help="One or more Ensembl IDs."
    )
    parser_info.add_argument(
        "--homology", "-H",
        default=False, 
        action='store_true',
        required=False, 
        help="Returns homology information of ID (default: False)."
    )
    parser_info.add_argument(
        "--xref", "-x",
        default=False, 
        action='store_true',
        required=False, 
        help="Returns information from external references (default: False)."
    )
    parser_info.add_argument(
        "--out", "-o",
        type=str,
        required=False,
        metavar="/ -o",
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json." 
            "Default: None (just prints results)."
        )
    )
    
    ## Show help when no arguments are given
    if len(sys.argv) == 1:
        parent_parser.print_help(sys.stderr)
        sys.exit(1)

    args = parent_parser.parse_args()

    ### Define return values
    ## Help return
    if args.help:
        help_.help_()
        
    ## Version return
    if args.version:        
        sys.stdout.write(f"gget version: {__init__.__version__}")
        
    ## ref return
    if args.command == "ref":
        # If list flag but no release passed, return all available species for latest release
        if args.list and args.release is None:
                # Find all available species for GTFs for this Ensembl release
                species_list_gtf = utils.ref_species_options('gtf')
                # Find all available species for FASTAs for this Ensembl release
                species_list_dna = utils.ref_species_options('dna') 

                # Find intersection of the two lists 
                # (Only species which have GTF and FASTAs available can continue)
                species_list = list(set(species_list_gtf) & set(species_list_dna))
                
                # Print available species list
                sys.stdout.write(f"{species_list}")
                
        # If list flag and release passed, return all available species for this release
        if args.list and args.release:
                # Find all available species for GTFs for this Ensembl release
                species_list_gtf = utils.ref_species_options('gtf', release=args.release)
                # Find all available species for FASTAs for this Ensembl release
                species_list_dna = utils.ref_species_options('dna', release=args.release) 

                # Find intersection of the two lists 
                # (Only species which have GTF and FASTAs available can continue)
                species_list = list(set(species_list_gtf) & set(species_list_dna))
                
                # Print available species list
                sys.stdout.write(f"{species_list}\n")
        
        # Raise error if neither species nor list flag passed
        if args.species is None and args.list is None:
            parser_ref.error("\n\nThe following arguments are required to fetch FTPs: -s/--species, e.g. '-s homo_sapiens'\n\n"
                             "gget ref --list -> lists out all available species. " 
                             "Combine with [-r] to define specific Ensembl release (default: latest release).")
        
        ## Clean up 'which' entry if passed
        if type(args.which) != str:
            which_clean = []
            # Split by comma (spaces are automatically split by nargs:"+")
            for which in args.which:
                which_clean.append(which.split(","))
            # Flatten which_clean
            which_clean_final = [item for sublist in which_clean for item in sublist]   
            # Remove empty strings resulting from split
            while("" in which_clean_final):
                which_clean_final.remove("")   
        else:
            which_clean_final = args.which

        if args.species:
            # Query Ensembl for requested FTPs using function ref
            ref_results = ref(args.species, which_clean_final, args.release, args.ftp)

            # Print or save list of URLs
            if args.ftp==True:
                if args.out:
                    os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
                    file = open(args.out, "w")
                    for element in ref_results:
                        file.write(element + "\n")
                    file.close()
                    sys.stderr.write(f"\nResults saved as {args.out}.\n")

                else:
                    results = " ".join(ref_results)
                    sys.stdout.write(f"{results}\n")

            # Print or save json file
            else:
                # Save in specified directory if -o specified
                if args.out:
                    os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
                    with open(args.out, 'w', encoding='utf-8') as f:
                        json.dump(ref_results, f, ensure_ascii=False, indent=4)
                    sys.stderr.write(f"\nResults saved as {args.out}.\n")
                # Print results if no directory specified
                else:
                    sys.stdout.write(f"{json.dumps(ref_results, ensure_ascii=False, indent=4)}\n")
                    sys.stderr.write("\nTo save these results, use flag '-o' in the format:\n" 
                                     "'-o path/to/directory/results.json'.\n")
        
    ## search return
    if args.command == "search":
        
        ## Clean up args.searchwords
        sw_clean = []
        # Split by comma (spaces are automatically split by nargs:"+")
        for sw in args.searchwords:
            sw_clean.append(sw.split(","))
        # Flatten which_clean
        sw_clean_final = [item for sublist in sw_clean for item in sublist]   
        # Remove empty strings resulting from split
        while("" in sw_clean_final) :
            sw_clean_final.remove("")  
        
        # Query Ensembl for genes based on species and searchwords using function search
        gget_results = search(sw_clean_final, args.species, operator=args.operator, limit=args.limit)
        
        # Save in specified directory if -o specified
        if args.out:
            os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
            gget_results.to_csv(args.out, index=False)
            sys.stderr.write(f"\nResults saved as {args.out}.\n")
        
        # Print results if no directory specified
        else:
            sys.stdout.write(f"{gget_results}")
            sys.stderr.write("\nTo save these results, use flag '-o' in the format:\n" 
                             "'-o path/to/directory/results.json'.\n")
            
    ## info return
    if args.command == "info":

        ## Clean up args.ens_ids
        ids_clean = []
        # Split by comma (spaces are automatically split by nargs:"+")
        for id_ in args.ens_ids:
            ids_clean.append(id_.split(","))
        # Flatten which_clean
        ids_clean_final = [item for sublist in ids_clean for item in sublist]   
        # Remove empty strings resulting from split
        while("" in ids_clean_final) :
            ids_clean_final.remove("")  

        # Look up requested Ensembl IDs
        info_results = info(ids_clean_final, homology=args.homology, xref=args.xref)

        # Print or save json file
        # Save in specified directory if -o specified
        if args.out:
            os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
            with open(args.out, 'w', encoding='utf-8') as f:
                json.dump(info_results, f, ensure_ascii=False, indent=4)
            sys.stderr.write(f"\nResults saved as {args.out}.\n")
        # Print results if no directory specified
        else:
            sys.stdout.write(f"{json.dumps(info_results, ensure_ascii=False, indent=4)}\n")
            sys.stderr.write("\nTo save these results, use flag '-o' in the format:\n" 
                             "'-o path/to/directory/results.json'.\n")