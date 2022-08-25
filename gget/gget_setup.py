import os
import shutil
import sys
import subprocess
import platform
import uuid
from platform import python_version

import logging

logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

from .compile import PACKAGE_PATH

## Variables for alphafold module
ALPHAFOLD_GIT_REPO = "https://github.com/deepmind/alphafold"
PDBFIXER_GIT_REPO = "https://github.com/openmm/pdbfixer.git"
# Unique ID to name temporary jackhmmer folder
UUID = "fcb45c67-8b27-4156-bbd8-9d11512babf2"
# Path to temporary mounted disk (global)
TMP_DISK = ""
# Model parameters
PARAMS_URL = (
    "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-03-02.tar"
)
PARAMS_DIR = os.path.join(PACKAGE_PATH, "bins/alphafold/")
PARAMS_PATH = os.path.join(PARAMS_DIR, "params_temp.tar")


def setup(module):
    """
    Function to install third-party dependencies for a specified gget module.

    Args:
    module      - (str) gget module for which dependencies should be installed, e.g. "alphafold"
    """
    supported_modules = ["alphafold"]
    if module not in supported_modules:
        raise ValueError(
            f"'module' argument specified as {module}. Expected one of: {', '.join(supported_modules)}"
        )

    if module == "alphafold":
        if platform.system() == "Windows":
            logging.warning(
                "gget setup alphafold and gget alphafold are not supported on Windows OS."
            )

        ## Make sure package paths are appended so openmm can be imported
        site_packages_path = os.__file__.split("os.py")[0] + "site-packages"
        if site_packages_path not in sys.path:
            sys.path.append(site_packages_path)

        site_packages_path_python = (
            "/".join(str(os.__file__.split("os.py")[0]).split("/")[:-2])
            + f"/python{'.'.join(python_version().split('.')[:2])}/site-packages"
        )
        if site_packages_path_python not in sys.path:
            sys.path.append(site_packages_path_python)

        site_packages_path_python37 = (
            "/".join(str(os.__file__.split("os.py")[0]).split("/")[:-2])
            + f"/python3.7/site-packages"
        )
        if site_packages_path_python37 not in sys.path:
            sys.path.append(site_packages_path_python37)

        conda_python_path = os.path.expanduser(
            f"~/opt/conda/lib/python{'.'.join(python_version().split('.')[:2])}/site-packages"
        )
        if conda_python_path not in sys.path:
            sys.path.append(conda_python_path)

        conda_python37_path = os.path.expanduser(
            f"~/opt/conda/lib/python3.7/site-packages"
        )
        if conda_python37_path not in sys.path:
            sys.path.append(conda_python37_path)

        # Global location of temporary disk
        global TMP_DISK

        ## Ask user to install openmm if not already installed
        try:
            import simtk.openmm as openmm

            # Silence openmm logger
            logging.getLogger("openmm").setLevel(logging.WARNING)

            # Check if correct version was installed
            if openmm.__version__ != "7.5.1":
                raise ImportError()

            logging.info(f"openmm v{openmm.__version__} already installed.")

        except ImportError:
            logging.error(
                """
        Please install AlphaFold third-party dependency openmm v7.5.1 by running the following command from the command line: 
        'conda install -qy conda==4.13.0 && conda install -qy -c conda-forge openmm=7.5.1' 
        (Recommendation: Follow with 'conda update -qy conda' to update conda to the latest version afterwards.)
        """
            )
            return

        ## Install Alphafold if not already installed
        logging.info("Installing AlphaFold from source (requires pip and git).")
        # Install AlphaFold and apply OpenMM patch.
        # command = f"""
        #     git clone {ALPHAFOLD_GIT_REPO} alphafold \
        #     && sed -i '' 's/\/tmp\/ramdisk/~\/tmp\/jackhmmer\/{UUID}/g' ./alphafold/alphafold/data/tools/jackhmmer.py \
        #     && pip install ./alphafold \
        #     && pushd {os.__file__.split('os.py')[0] + 'site-packages/'} \
        #     && patch -p0 < /content/alphafold/docker/openmm.patch \
        #     && popd \
        #     && rm -rf alphafold
        #     """

        ## Install AlphaFold and change jackhmmer directory where database chunks are saved in
        # Define AlphaFold folder name and location
        alphafold_folder = os.path.join(
            PACKAGE_PATH, "tmp_alphafold_" + str(uuid.uuid4())
        )
        # command = """
        #     git clone -q {} {} \
        #     && sed -i '' 's/\/tmp\/ramdisk/{}/g' {}/alphafold/data/tools/jackhmmer.py \
        #     && sed -i '' '/{}/d' {}/alphafold/data/tools/jackhmmer.py \
        #     && pip install -q {} \
        #     """.format(
        #     ALPHAFOLD_GIT_REPO,
        #     alphafold_folder,
        #     os.path.expanduser(f"~/tmp/jackhmmer/{UUID}").replace(
        #         "/", "\/"
        #     ),  # Replace directory where jackhmmer database chunks will be saved
        #     alphafold_folder,
        #     "logging.info",  # Delete all info loggers
        #     alphafold_folder,
        #     alphafold_folder,
        # )

        # Clone AlphaFold github repo
        # Replace directory where jackhmmer database chunks will be saved
        # Insert "logging.set_verbosity(logging.WARNING)" to mute all info loggers
        # Pip install AlphaFold from local directory
        if platform.system() == "Darwin":
            command = """
                git clone -q {} {} \
                && sed -i '' 's/\/tmp\/ramdisk/{}/g' {}/alphafold/data/tools/jackhmmer.py \
                && sed -i '' 's/from absl import logging/from absl import logging\\\nlogging.set_verbosity(logging.WARNING)/g' {}/alphafold/data/tools/jackhmmer.py \
                && pip install -q {} \
                """.format(
                ALPHAFOLD_GIT_REPO,
                alphafold_folder,
                os.path.expanduser(f"~/tmp/jackhmmer/{UUID}").replace(
                    "/", "\/"
                ),  # Replace directory where jackhmmer database chunks will be saved
                alphafold_folder,
                alphafold_folder,
                alphafold_folder,
            )
        else:
            command = """
                git clone -q {} {} \
                && sed -i 's/\/tmp\/ramdisk/{}/g' {}/alphafold/data/tools/jackhmmer.py \
                && sed -i 's/from absl import logging/from absl import logging\\\nlogging.set_verbosity(logging.WARNING)/g' {}/alphafold/data/tools/jackhmmer.py \
                && pip install -q {} \
                """.format(
                ALPHAFOLD_GIT_REPO,
                alphafold_folder,
                os.path.expanduser(f"~/tmp/jackhmmer/{UUID}").replace(
                    "/", "\/"
                ),  # Replace directory where jackhmmer database chunks will be saved
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
            logging.error("AlphaFold installation failed.")
            return

        # Remove cloned directory
        shutil.rmtree(alphafold_folder)

        try:
            import alphafold as AlphaFold

            logging.info(f"AlphaFold installed succesfully.")
        except ImportError:
            logging.error("AlphaFold installation failed.")
            return

        ## Append AlphaFold to path
        alphafold_path = os.path.abspath(os.path.dirname(AlphaFold.__file__))
        if alphafold_path not in sys.path:
            sys.path.append(alphafold_path)

        ## Install pdbfixer v1.7 (compadible with openmm v7.5.1)
        logging.info("Installing pdbfixer from source (requires pip and git).")

        pdbfixer_folder = os.path.join(
            PACKAGE_PATH, "tmp_pdbfixer_" + str(uuid.uuid4())
        )
        command = f"""
            git clone -q --branch v1.7 {PDBFIXER_GIT_REPO} {pdbfixer_folder} \
            && pip install -q {pdbfixer_folder} \
            """

        with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
            stderr = process.stderr.read().decode("utf-8")
        # Exit system if the subprocess returned with an error
        if process.wait() != 0:
            if stderr:
                # Log the standard error if it is not empty
                sys.stderr.write(stderr)
            logging.error("pdbfixer installation failed.")
            return

        # Remove cloned directory
        shutil.rmtree(pdbfixer_folder)

        # Check if pdbfixer was installed successfully
        command = "pip list | grep pdbfixer"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        pdb_out, err = process.communicate()

        if pdb_out.decode() != "":
            logging.info(f"pdbfixer installed succesfully.")
        else:
            logging.error("pdbfixer installation failed.")
            return

        # ## Create a temporary file system (TMPFS) to store a database chunk to make Jackhmmer run fast.
        # # TMPFS uses local memory for file system reads and writes, which is typically much faster than reads and writes in a UFS file system.
        # logging.info("Creating temporary file system (TMPFS) to store a database chunk and make Jackhmmer run faster.")
        # if platform.system() == "Linux":
        #   command = f"mkdir -m 777 -p /tmp/ramdisk && mount -t tmpfs -o size=9G ramdisk /tmp/ramdisk"

        #   with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
        #       stderr = process.stderr.read().decode("utf-8")

        #   # Exit system if the subprocess returned with an error
        #   if process.wait() != 0:
        #       # if stderr:
        #       # # Log the standard error if it is not empty
        #       # sys.stderr.write(stderr)
        #       logging.warning("Creating TMPFS failed. Jackhmmer will run slower.")

        # elif platform.system() == "Darwin":
        #   # Attach disk with 9GB
        #   command1 = "hdiutil attach -nomount ram://18432000"
        #   process = subprocess.Popen(command1, shell=True, stdout=subprocess.PIPE)
        #   hdi_out, err = process.communicate()

        #   # Record number of new disk
        #   TMP_DISK = hdi_out.decode("utf-8").strip()
        #   DISK_NUMBER = f"$({TMP_DISK} | tr -dc '0-9')"

        #   # Set up TMPFS
        #   command2 = f"newfs_hfs -v tmp /dev/rdisk{DISK_NUMBER}"
        #   command3 = f"diskutil eraseVolume HFS+ /tmp/ramdisk {TMP_DISK}"

        #   command = command2 + " && " + command3

        #   with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
        #     stderr = process.stderr.read().decode("utf-8")

        #   # Exit system if the subprocess returned with an error
        #   if process.wait() != 0:
        #       # if stderr:
        #       # # Log the standard error if it is not empty
        #       # sys.stderr.write(stderr)
        #       logging.warning("Creating TMPFS failed. Jackhmmer will run slower.")

        # else:
        # # Create folder to save temporary jackhmmer database chunks in
        # !!! Delete the duplicate line from get_msa function in gget_alphafold.py if this is uncommented
        # os.makedirs(f"~/tmp/jackhmmer/{UUID}", exist_ok=True)

        ## Download model parameters
        # Download parameters if the params directory is empty
        if not os.path.exists(os.path.join(PARAMS_DIR, "params/")):
            # Create folder to save parameter files
            os.makedirs(os.path.join(PARAMS_DIR, "params/"), exist_ok=True)

        if len(os.listdir(os.path.join(PARAMS_DIR, "params/"))) < 12:
            logging.info(
                "Downloading AlphaFold model parameters (requires 4.1 GB of storage). This might take a few minutes."
            )
            if platform.system() == "Windows":
                # The double-quotation marks allow white spaces in the path, but this does not work for Windows
                command = f"""
                    curl -# -o {PARAMS_PATH} {PARAMS_URL} \
                    && tar --extract --file={PARAMS_PATH} --directory={PARAMS_DIR+'params/'} --preserve-permissions \
                    && rm {PARAMS_PATH}
                    """
            else:
                command = f"""
                    curl -# -o '{PARAMS_PATH}' '{PARAMS_URL}' \
                    && tar --extract --file='{PARAMS_PATH}' --directory='{PARAMS_DIR+'params/'}' --preserve-permissions \
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
                logging.error("Model parameter download failed.")
                return
            else:
                logging.info("Model parameter download complete.")
        else:
            logging.info("AlphaFold model parameters already downloaded.")
