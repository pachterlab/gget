import os
import logging
import shutil
import sys
import subprocess
import platform
import uuid
from platform import python_version

from .utils import set_up_logger

logger = set_up_logger()

from .compile import PACKAGE_PATH
from .constants import (
    ELM_INSTANCES_FASTA_DOWNLOAD,
    ELM_CLASSES_TSV_DOWNLOAD,
    ELM_INSTANCES_TSV_DOWNLOAD,
    ELM_INTDOMAINS_TSV_DOWNLOAD,
)

## Variables for elm module
ELM_FILES = os.path.join(PACKAGE_PATH, "elm_files")
ELM_INSTANCES_FASTA = os.path.join(ELM_FILES, "elm_instances.fasta")
ELM_CLASSES_TSV = os.path.join(ELM_FILES, "elms_classes.tsv")
ELM_INSTANCES_TSV = os.path.join(ELM_FILES, "elm_instances.tsv")
ELM_INTDOMAINS_TSV = os.path.join(ELM_FILES, "elm_interaction_domains.tsv")

## Variables for alphafold module
ALPHAFOLD_GIT_REPO = "https://github.com/deepmind/alphafold"
ALPHAFOLD_GIT_REPO_VERSION = "main"  # Get version currently hosted on main branch
PDBFIXER_GIT_REPO = "https://github.com/openmm/pdbfixer.git"
# Unique ID to name temporary jackhmmer folder
UUID = "fcb45c67-8b27-4156-bbd8-9d11512babf2"
# # Path to temporary mounted disk (global)
# TMP_DISK = ""
# Model parameters
PARAMS_URL = (
    "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar"
)
PARAMS_DIR = os.path.join(PACKAGE_PATH, "bins/alphafold/")
PARAMS_PATH = os.path.join(PARAMS_DIR, "params_temp.tar")


def _install(package: str, import_name: str, verbose: bool = True):
    if verbose:
        logger.info(f"Installing {package} package (requires pip).")
    command = f"pip install -q -U {package}"
    with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
        stderr = process.stderr.read().decode("utf-8")
    # Exit system if the subprocess returned with an error
    if process.wait() != 0:
        if stderr:
            # Log the standard error if it is not empty
            sys.stderr.write(stderr)
        logger.error(
            f"{package} installation with pip (https://pypi.org/project/{package}) failed."
        )
        return

    # Test installation
    try:
        exec(f"import {import_name}")

        if verbose:
            logger.info(f"{import_name} installed succesfully.")
    except ImportError as e:
        logger.error(
            f"{package} installation with pip (https://pypi.org/project/{package}) failed. Import error:\n{e}"
        )
        return


def setup(module, verbose=True, out=None):
    """
    Function to install third-party dependencies for a specified gget module.
    Some modules require pip to be installed (https://pip.pypa.io/en/stable/installation).
    Some modules require curl to be installed (https://everything.curl.dev/get).

    Args:
    - module    (str) gget module for which dependencies should be installed, e.g. "alphafold", "cellxgene", "elm", "gpt", or "cbio".
    - verbose   True/False whether to print progress information. Default True.
    - out       (str) Path to directory to save downloaded files in (currently only applies when module='elm').
                NOTE: Do not use this argument when downloading the files for use with 'gget.elm'.
                Default None (files are saved in the gget installation directory).
    """
    supported_modules = ["alphafold", "cellxgene", "elm", "gpt", "cbio"]
    if module not in supported_modules:
        raise ValueError(
            f"'module' argument specified as {module}. Expected one of: {', '.join(supported_modules)}"
        )

    if module == "gpt":
        if verbose:
            logger.info("Installing openai version <=0.28.1 (requires pip).")
        command = "pip install -q -U 'openai<=0.28.1'"
        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            if stderr:
                # Log the standard error if it is not empty
                sys.stderr.write(stderr)
            logger.error(
                "Installation of openai version <=0.28.1 with pip (https://pypi.org/project/openai) failed."
            )
            return

        # Test installation
        try:
            import openai

            if verbose:
                logger.info(f"openai installed succesfully.")
        except ImportError as e:
            logger.error(
                f"openai installation with pip (https://pypi.org/project/openai) failed. Import error:\n{e}"
            )
            return

    elif module == "cellxgene":
        _install("cellxgene-census", "cellxgene_census", verbose=verbose)

    elif module == "elm":
        if verbose:
            logger.info(
                "ELM data can be downloaded & distributed for non-commercial use according to the following license: http://elm.eu.org/media/Elm_academic_license.pdf"
            )
            logger.info(
                "Downloading ELM database files (requires curl to be installed)..."
            )

        if out is not None:
            elm_files_out = os.path.abspath(out)
            elm_instances_fasta = os.path.join(elm_files_out, "elm_instances.fasta")
            elm_classes_tsv = os.path.join(elm_files_out, "elms_classes.tsv")
            elm_instances_tsv = os.path.join(elm_files_out, "elm_instances.tsv")
            elm_intdomains_tsv = os.path.join(
                elm_files_out, "elm_interaction_domains.tsv"
            )

            # Create folder for ELM files (if it does not exist)
            if not os.path.exists(elm_files_out):
                os.makedirs(elm_files_out)

        else:
            elm_instances_fasta = ELM_INSTANCES_FASTA
            elm_classes_tsv = ELM_CLASSES_TSV
            elm_instances_tsv = ELM_INSTANCES_TSV
            elm_intdomains_tsv = ELM_INTDOMAINS_TSV

            # Create folder for ELM files (if it does not exist)
            if not os.path.exists(ELM_FILES):
                os.makedirs(ELM_FILES)

        if platform.system() == "Windows":
            # The double-quotation marks allow white spaces in the path, but this does not work for Windows
            command = f"""
                curl -o {elm_instances_fasta} \"{ELM_INSTANCES_FASTA_DOWNLOAD}\" \\
                &&  curl -o {elm_classes_tsv} \"{ELM_CLASSES_TSV_DOWNLOAD}\" \\
                &&  curl -o {elm_instances_tsv} \"{ELM_INSTANCES_TSV_DOWNLOAD}\" \\
                &&  curl -o {elm_intdomains_tsv} \"{ELM_INTDOMAINS_TSV_DOWNLOAD}\"
                """

        else:
            command = f"""
                curl -o '{elm_instances_fasta}' {ELM_INSTANCES_FASTA_DOWNLOAD} \\
                &&  curl -o '{elm_classes_tsv}' {ELM_CLASSES_TSV_DOWNLOAD} \\
                &&  curl -o '{elm_instances_tsv}' {ELM_INSTANCES_TSV_DOWNLOAD} \\
                &&  curl -o '{elm_intdomains_tsv}' '{ELM_INTDOMAINS_TSV_DOWNLOAD}'
                """

        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
            # Log the standard error if it is not empty
            if stderr:
                sys.stderr.write(stderr)

        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            logger.error("ELM database files download failed.")
            return

        # Check if files are present
        if os.path.exists(elm_instances_fasta):
            if verbose:
                logger.info(f"ELM sequences file present.")
        else:
            logger.error("ELM FASTA file missing.")

        if os.path.exists(elm_classes_tsv):
            if verbose:
                logger.info("ELM classes file present.")
        else:
            logger.error("ELM classes file missing.")

        if os.path.exists(elm_instances_tsv):
            if verbose:
                logger.info("ELM instances file present.")
        else:
            logger.error("ELM instances file missing.")

        if os.path.exists(elm_intdomains_tsv):
            if verbose:
                logger.info("ELM interactions domains file present.")
        else:
            logger.error("ELM interactions domains file missing.")

    elif module == "alphafold":
        if platform.system() == "Windows":
            logger.error(
                "gget setup alphafold and gget alphafold are not supported on Windows OS."
            )

        ## Ask user to install openmm if not already installed
        try:
            import simtk.openmm as openmm

            # Silence openmm logger
            logging.getLogger("openmm").setLevel(logging.WARNING)

            # Commenting the following out because openmm v7.7.0 does not support __version__
            # # Check if correct version was installed
            # if openmm.__version__ != "7.5.1":
            #     raise ImportError()

            # if verbose:
            #   logger.info(f"openmm v{openmm.__version__} already installed.")

        except ImportError as e:
            raise ImportError(
                f"""
                Trying to import openmm resulted in the following error:
                {e}

                Please install AlphaFold third-party dependency openmm v7.5.1 (or v7.7.0 for Python >= 3.10) by running the following command from the command line: 
                'conda install -qy conda==4.13.0 && conda install -qy -c conda-forge openmm=7.5.1'   (or 'openmm=7.7.0' for Python >= 3.10)
                (Recommendation: Follow with 'conda update -qy conda' to update conda to the latest version afterwards.)
                """
            )

        ## Install py3Dmol
        if verbose:
            logger.info("Installing py3Dmol (requires pip).")
        command = "pip install -q py3Dmol"
        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            if stderr:
                # Log the standard error if it is not empty
                sys.stderr.write(stderr)
            logger.error(
                "py3Dmol installation with pip (https://pypi.org/project/py3Dmol) failed."
            )
            return

        # Test installation
        try:
            import py3Dmol

            if verbose:
                logger.info(f"py3Dmol installed succesfully.")
        except ImportError as e:
            logger.error(
                f"py3Dmol installation with pip (https://pypi.org/project/py3Dmol/) failed. Import error:\n{e}"
            )
            return

        ## Install Alphafold if not already installed
        if verbose:
            logger.info("Installing AlphaFold from source (requires pip and git).")

        ## Install AlphaFold and change jackhmmer directory where database chunks are saved in
        # Define AlphaFold folder name and location
        alphafold_folder = os.path.join(
            PACKAGE_PATH, "tmp_alphafold_" + str(uuid.uuid4())
        )

        # Clone AlphaFold github repo
        # Replace directory where jackhmmer database chunks will be saved
        # Insert "logging.set_verbosity(logging.WARNING)" to mute all info loggers
        # Pip install AlphaFold from local directory
        if platform.system() == "Darwin":
            command = """
                git clone --branch main -q --branch {} {} {} \\
                && sed -i '' 's/\\/tmp\\/ramdisk/{}/g' {}/alphafold/data/tools/jackhmmer.py \\
                && sed -i '' '/from absl import logging/a logging.set_verbosity(logging.WARNING)' {}/alphafold/data/tools/jackhmmer.py \\
                && pip install -q -r {}/requirements.txt \\
                && pip install -q --no-dependencies {}
                """.format(
                ALPHAFOLD_GIT_REPO_VERSION,
                ALPHAFOLD_GIT_REPO,
                alphafold_folder,
                os.path.expanduser(f"~/tmp/jackhmmer/{UUID}").replace(
                    "/", "\\/"
                ),  # Replace directory where jackhmmer database chunks will be saved
                alphafold_folder,
                alphafold_folder,
                alphafold_folder,
                alphafold_folder,
            )
        else:
            command = """
                git clone --branch main -q --branch {} {} {} \\
                && sed -i 's/\\/tmp\\/ramdisk/{}/g' {}/alphafold/data/tools/jackhmmer.py \\
                && sed -i 's/from absl import logging/from absl import logging\\\nlogging.set_verbosity(logging.WARNING)/g' {}/alphafold/data/tools/jackhmmer.py \\
                && pip install -q -r {}/requirements.txt \\
                && pip install -q --no-dependencies {}
                """.format(
                ALPHAFOLD_GIT_REPO_VERSION,
                ALPHAFOLD_GIT_REPO,
                alphafold_folder,
                os.path.expanduser(f"~/tmp/jackhmmer/{UUID}").replace(
                    "/", "\\/"
                ),  # Replace directory where jackhmmer database chunks will be saved
                alphafold_folder,
                alphafold_folder,
                alphafold_folder,
                alphafold_folder,
            )

        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            if stderr:
                # Log the standard error if it is not empty
                sys.stderr.write(stderr)
            logger.error("AlphaFold installation failed.")
            return

        # Remove cloned directory
        shutil.rmtree(alphafold_folder)

        try:
            import alphafold as AlphaFold

            if verbose:
                logger.info(f"AlphaFold installed succesfully.")
        except ImportError as e:
            logger.error(f"AlphaFold installation failed. Import error:\n{e}")
            return

        ## Append AlphaFold to path
        alphafold_path = os.path.abspath(os.path.dirname(AlphaFold.__file__))
        if alphafold_path not in sys.path:
            sys.path.append(alphafold_path)

        ## Install pdbfixer
        if verbose:
            logger.info("Installing pdbfixer from source (requires pip and git).")

        pdbfixer_folder = os.path.join(
            PACKAGE_PATH, "tmp_pdbfixer_" + str(uuid.uuid4())
        )

        try:
            if openmm.__version__ == "7.5.1":
                # Install pdbfixer version compatible with openmm v7.5.1
                PDBFIXER_VERSION = "v1.7"
            else:
                PDBFIXER_VERSION = "v1.8.1"
        except:
            PDBFIXER_VERSION = "v1.8.1"

        command = f"""
            git clone -q --branch {PDBFIXER_VERSION} {PDBFIXER_GIT_REPO} {pdbfixer_folder} \\
            && pip install -q {pdbfixer_folder}
            """

        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            if stderr:
                # Log the standard error if it is not empty
                sys.stderr.write(stderr)
            logger.error("pdbfixer installation failed.")
            return

        # Remove cloned directory
        shutil.rmtree(pdbfixer_folder)

        # Check if pdbfixer was installed successfully
        command = "pip list | grep pdbfixer"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        pdb_out, err = process.communicate()

        if pdb_out.decode() != "":
            logger.info(f"pdbfixer installed succesfully.")
        else:
            logger.error("pdbfixer installation failed.")
            return

        ## Download model parameters
        # Download parameters if the params directory is empty
        if not os.path.exists(os.path.join(PARAMS_DIR, "params/")):
            # Create folder to save parameter files
            os.makedirs(os.path.join(PARAMS_DIR, "params/"), exist_ok=True)

        if len(os.listdir(os.path.join(PARAMS_DIR, "params/"))) < 12:
            if verbose:
                logger.info(
                    "Downloading AlphaFold model parameters (requires 4.1 GB of storage). This might take a few minutes."
                )
            if platform.system() == "Windows":
                # The double-quotation marks allow white spaces in the path, but this does not work for Windows
                command = f"""
                    curl -# -o {PARAMS_PATH} {PARAMS_URL} \\
                    && tar --extract --file={PARAMS_PATH} --directory={PARAMS_DIR+'params/'} --preserve-permissions \\
                    && rm {PARAMS_PATH}
                    """
            else:
                command = f"""
                    curl -# -o '{PARAMS_PATH}' '{PARAMS_URL}' \\
                    && tar --extract --file='{PARAMS_PATH}' --directory='{PARAMS_DIR+'params/'}' --preserve-permissions \\
                    && rm '{PARAMS_PATH}'
                    """

            with subprocess.Popen(
                command, shell=True, stderr=subprocess.PIPE
            ) as process:
                stderr = process.stderr.read().decode("utf-8")
                # Log the standard error if it is not empty
                if stderr:
                    sys.stderr.write(stderr)
            # Exit system if the subprocess returned with an error
            if process.wait() != 0:
                logger.error("Model parameter download failed.")
                return
            else:
                if verbose:
                    logger.info("Model parameter download complete.")
        else:
            if verbose:
                logger.info("AlphaFold model parameters already downloaded.")

    elif module == "cbio":
        _install("bravado", "bravado", verbose=verbose)
