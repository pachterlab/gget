import hashlib
import gzip
import json
import math
import os
import subprocess
from typing import Literal, TypeVar

import pandas as pd
import numpy as np

import requests
from collections import defaultdict, OrderedDict

from .utils import set_up_logger

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm

from .constants import CBIO_CANCER_TYPE_TO_TISSUE_DICTIONARY

_K = TypeVar("_K")
_V = TypeVar("_V")

logger = set_up_logger()


def _ints_between(
    start: int, end: int, max_count: int, min_count: int, verbose: bool = False
) -> list[int]:
    """
    Generate a list of integers between start and end (inclusive) with a maximum count of max_count and a minimum count min_count.
    The list is guaranteed to contain start and end, and the spacing between the numbers will be as even as possible.
    If a perfect spacing is not possible, the spacing will omit a number rather than overcrowding.

    This method is designed to be used for plot labels
    """
    assert max_count >= 2, "max_count must be at least 2"
    assert start < end, "start must be less than end"

    if max_count == 2:
        return [start, end]
    else:
        step = int(math.ceil((end - start) / (max_count - 1)))
        if verbose:
            print(f"Original step: {step}, {(end - start) % step=}")
            print(f"{end - start=}")
        # check if it comes out even
        while ((end - start) / step) + 1 >= min_count and (end - start) % step != 0:
            step += 1
            if verbose:
                print(f"New step: {step}, {(end - start) % step=}")
        if (end - start) % step != 0:
            step = int(math.ceil((end - start) / (max_count - 1)))
            if verbose:
                print("Reverted step")

        out = [start]
        current = start
        while current + step * 2 <= end:
            current += step
            out.append(current)
        out.append(end)
        if verbose:
            print(f"{out=}")
        return out


def _describe_bytes(size: int) -> str:
    """
    Describe a size in bytes in human-readable format.

    :param size:    size in bytes
    """

    steps = ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"]

    unit = steps.pop(0)

    while size >= 1000 and steps:
        size /= 1000
        unit = steps.pop(0)

    if unit == "bytes":
        return f"{size} {unit}"
    return f"{size:.2f} {unit}"


def _download_file_from_git_lfs(
    target_path: str, oid: str, size: int, verbose: bool = False
) -> bool:
    """
    Download a single object from Git LFS.

    :param target_path:  path to save the downloaded object
    :param oid:          object ID
    :param size:         size of the object

    :return:             True if the object was downloaded successfully, False otherwise
    """
    lfs_metadata = {
        "operation": "download",
        "transfer": ["basic"],
        "objects": [{"oid": oid, "size": size}],
    }
    lfs_metadata_json = json.dumps(lfs_metadata)

    try:
        github_url = f"https://github.com/cBioPortal/datahub.git/info/lfs/objects/batch"

        curl_command = [
            "curl",
            "-X", "POST",
            "-H", "Accept: application/vnd.git-lfs+json",
            "-H", "Content-Type: application/json",
            "-d", lfs_metadata_json,
            github_url,
        ]

        result = subprocess.run(curl_command, capture_output=True, text=True)
        response_json = json.loads(result.stdout)

        href = response_json["objects"][0]["actions"]["download"]["href"]

        response = requests.get(href)

        if not response.ok:
            logger.error(f"Failed to download object {oid} to {target_path}")
            return False

        with open(target_path, "wb") as f:
            f.write(response.content)

        if verbose:
            logger.info(f"Downloaded object {oid} to {target_path}")

    except Exception as e:
        logger.error(f"Error downloading object {oid} to {target_path}: {e}")
        return False

    return True


class _LFSDownloadPlan:
    def __init__(self, verbose: bool = False):
        self.objects: list[tuple[str, tuple[str, int]]] = []
        """(target_path, (oid, size))"""

        self.verbose = verbose

    @property
    def total_size(self) -> int:
        return sum(size for _, (_, size) in self.objects)

    def add(self, target_path: str, oid: str, size: int):
        self.objects.append((target_path, (oid, size)))

    def download(self) -> bool:
        """
        Download all objects in the plan.

        :return:    True if all objects were downloaded successfully, False otherwise
        """
        success = True

        for target_path, (oid, size) in self.objects:
            success &= _download_file_from_git_lfs(target_path, oid, size)

        return success


def download_cbioportal_data(
    study_ids: list[str],
    verbose: bool = False,
    out_dir: str | None = None,
    confirm_download: bool = False,
) -> bool:
    """
    Download data from cBioPortal studies.

    Args:

    :param study_ids:           list of cBioPortal study IDs
    :param verbose:             whether to print out progress
    :param out_dir:             output directory to save the data, default `<current_directory>/gget_cbio_cache/`
    :param confirm_download:    whether to confirm the download before proceeding

    Return:

    :return:            True if successfully downloaded all needed data, False otherwise
    """

    actual_out_dir = os.path.abspath(out_dir or "gget_cbio_cache")

    os.makedirs(actual_out_dir, exist_ok=True)

    success = True

    if confirm_download:
        plan = _LFSDownloadPlan(verbose=verbose)
    else:
        plan = None

    for study_id in study_ids:
        file_types = ["mutations", "cna", "sv", "clinical_patient", "clinical_sample"]

        for file_type in file_types:
            try:
                data_folder = os.path.join(actual_out_dir, study_id)
                filename = os.path.join(data_folder, f"{file_type}.txt")

                if not os.path.exists(data_folder):
                    os.makedirs(data_folder, exist_ok=True)

                if os.path.exists(filename):
                    if verbose:
                        logger.info(f"File {filename} already exists, skipping.")
                    continue

                # URL to be downloaded
                url = f"https://raw.githubusercontent.com/cBioPortal/datahub/master/public/{study_id}/data_{file_type}.txt"

                response = requests.get(url)

                if not response.ok:
                    logger.error(
                        f"Failed to download {file_type} data for study {study_id}"
                    )
                    success = False
                    continue

                lines = response.content.decode().splitlines(keepends=True)

                fields: dict[str, str] = {}

                for line in lines:
                    k, v = line.split(" ", 1)
                    k = k.strip()
                    v = v.strip()
                    fields[k] = v

                assert (
                    fields["version"] == "https://git-lfs.github.com/spec/v1"
                ), f"Cannot handle git-lfs version {fields['version']}"

                oid: str = fields["oid"].split(":")[1].strip()
                size: int = int(fields["size"])

                if plan:
                    plan.add(filename, oid, size)
                else:
                    success &= _download_file_from_git_lfs(
                        filename, oid, size, verbose=verbose
                    )

            except Exception as e:
                logger.error(
                    f"Error downloading {file_type} data for study {study_id}: {e}"
                )
                success = False

        if verbose and not confirm_download:
            logger.info(f"Downloaded all data for study {study_id}")

    # If using a download plan AND there are actually objects to download, ask for confirmation
    if plan and plan.objects:
        do_download = (
            input(
                f"Do you want to download {_describe_bytes(plan.total_size)} to {actual_out_dir}? (y/n) "
            )
            .lower()
            .strip()
            == "y"
        )
        if do_download:
            success &= plan.download()
        else:
            print("Download aborted")
            return False
    elif verbose:
        logger.info("Downloaded all data for all studies")

    return success


def _extract_study_name(name: str) -> str:
    parenthesis_index = name.find("(")
    if parenthesis_index != -1:
        return name[:parenthesis_index].strip()
    return name


def find_study_ids_by_keywords(key_words: list[str] | set[str]) -> list[str]:
    """
    Find cBioPortal study IDs by keyword.

    Args:

    :param key_words:    list of keywords to search for

    Return:

    :return:            list of study IDs that match the keywords
    """

    try:
        from bravado.client import SwaggerClient
    except ImportError:
        logger.error(
            """
            Some third-party dependencies are missing. Please run the following command: 
            >>> gget.setup('cbio') or $ gget setup cbio

            Alternative: Install the bravado package using pip (https://pypi.org/project/bravado).
            """
        )
        return []

    api = SwaggerClient.from_url(
        "https://www.cbioportal.org/api/v2/api-docs",
        config={
            "validate_requests": False,
            "validate_responses": False,
            "validate_swagger_spec": False,
        },
    )

    for a in dir(api):
        setattr(api, a.replace(" ", "_").lower(), getattr(api, a))

    studies = api.Studies.getAllStudiesUsingGET().result()

    cancer_type_acronym_dict = {
        _extract_study_name(individual_study["name"]): individual_study["cancerTypeId"]
        for individual_study in studies
    }
    cancer_type_acronym_dict = OrderedDict(sorted(cancer_type_acronym_dict.items()))

    cancer_id_list = [
        cancer_type_acronym
        for cancer_type, cancer_type_acronym in cancer_type_acronym_dict.items()
        if any(
            key_word in cancer_type.lower() or key_word in cancer_type_acronym.lower()
            for key_word in key_words
        )
        and cancer_type_acronym.lower() != "mixed"
    ]

    matching_study_ids = [
        individual_study["studyId"]
        for individual_study in studies
        if individual_study["cancerTypeId"] in cancer_id_list
    ]

    return matching_study_ids


def _get_ensembl_gene_id(transcript_id: str, verbose: bool = False):
    try:
        url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
        response = requests.get(url, headers={"Content-Type": "application/json"})

        if not response.ok:
            response.raise_for_status()

        data = response.json()

        return data.get("Parent")
    except Exception as e:
        if verbose:
            print(f"Error for: {transcript_id}")
        return "Unknown"


def _get_valid_ensembl_gene_id(
    row, transcript_column: str = "seq_ID", gene_column: str = "gene_name"
):
    ensembl_gene_id = _get_ensembl_gene_id(row[transcript_column])
    if ensembl_gene_id == "Unknown":
        return row[gene_column]
    return ensembl_gene_id


def _nested_defaultdict() -> defaultdict[_K, _V | defaultdict[_K]]:
    return defaultdict(_nested_defaultdict)


_ENSEMBL = "Ensembl"
_SYMBOL = "Symbol"


class _GeneAnalysis:
    def __init__(
        self,
        study_ids: list[str],
        genes: list[str],
        merge_type: Literal["Symbol", "Ensembl"] = "Symbol",
        remove_non_ensembl_genes: bool = False,
        data_dir: str = "gget_cbio_cache",
        figure_output_dir: str = "gget_cbio_figures",
        verbose: bool = False
    ):
        self.study_ids = study_ids
        self.genes = genes
        self.merge_type = merge_type
        self.remove_non_ensembl_genes = remove_non_ensembl_genes
        self.data_dir = data_dir
        self.figure_output_dir = figure_output_dir
        self.verbose = verbose

        os.makedirs(self.figure_output_dir, exist_ok=True)

        assert self.merge_type in [
            _SYMBOL,
            _ENSEMBL,
        ], "merge_type must be either 'Symbol' or 'Ensembl'"

        self.columns_to_keep = [
            "Tumor_Sample_Barcode",
            "Hugo_Symbol",
            "Entrez_Gene_Id",
            "Consequence",
        ]
        self.column_for_merging: str = (
            "Hugo_Symbol" if self.merge_type == _SYMBOL else "Ensembl_Gene_ID"
        )

        self.df_collection: dict[str, dict[str, pd.DataFrame]] = {}
        self.big_combined_df = self._create_or_load_study_dataframes()

    def _create_single_study_dataframe(self, study_id: str) -> pd.DataFrame:
        data_folder = os.path.join(self.data_dir, study_id)

        mutation_df = pd.read_csv(os.path.join(data_folder, "mutations.txt"), sep="\t")
        sample_df = pd.read_csv(
            os.path.join(data_folder, "clinical_sample.txt"), sep="\t"
        )

        self.df_collection[study_id]["mutations"] = mutation_df
        self.df_collection[study_id]["samples"] = sample_df

        # Merge on gene name
        if self.merge_type == _ENSEMBL:
            if (
                "Gene" in mutation_df.columns
                and mutation_df["Gene"].str.startswith("ENSG").any()
            ):
                mutation_df.rename(columns={"Gene": "Ensembl_Gene_ID"}, inplace=True)
                if self.remove_non_ensembl_genes:
                    mutation_df = mutation_df[
                        mutation_df["Ensemble_Gene_ID"].str.startswith("ENSG")
                    ]
            elif (
                "Transcript_ID" in mutation_df.columns
                and mutation_df["Transcript_ID"].str.startswith("ENST").any()
            ):
                if self.verbose:
                    logger.info("Fetching gene IDs from Ensembl")
                mutation_df.progress_apply(
                    _get_valid_ensembl_gene_id,
                    axis=1,
                    transcript_column="Transcript_ID",
                    gene_column="Hugo_Symbol",
                )
                if self.remove_non_ensembl_genes:
                    mutation_df = mutation_df[
                        mutation_df["Ensemble_Gene_ID"].str.startswith("ENSG")
                    ]
            else:
                self.merge_type = _SYMBOL
                logger.warn(
                    "No Ensembl gene IDs found in the mutation data. Merging on gene symbol instead."
                )

        if self.merge_type == _ENSEMBL:
            self.column_for_merging = "Ensembl_Gene_ID"

            aggregated_df = (
                mutation_df.groupby(["Tumor_Sample_Barcode", self.column_for_merging])
                .agg(
                    {
                        "Hugo_Symbol": lambda x: ",".join(x.unique()),
                        "Entrez_Gene_Id": lambda x: ",".join(map(str, x.unique())),
                        "Consequence": lambda x: ",".join(x.unique()),
                    }
                )
                .reset_index()
            )
        elif self.merge_type == _SYMBOL:
            self.column_for_merging = "Hugo_Symbol"

            aggregated_df = (
                mutation_df.groupby(["Tumor_Sample_Barcode", self.column_for_merging])
                .agg(
                    {
                        "Entrez_Gene_Id": lambda x: ",".join(map(str, x.unique())),
                        "Consequence": lambda x: ",".join(x.unique()),
                    }
                )
                .reset_index()
            )
        else:
            raise AssertionError(f"Invalid merge type: {self.merge_type}")

        if self.column_for_merging not in self.columns_to_keep:
            self.columns_to_keep.append(self.column_for_merging)

        aggregated_df["Consequence"] = aggregated_df["Consequence"].apply(
            lambda x: "Multiple_consequences" if "," in x else x
        )

        occurrences_df = (
            aggregated_df.groupby([self.column_for_merging, "Tumor_Sample_Barcode"])
            .size()
            .reset_index(name="mutation_occurrences")
        )

        final_df = pd.merge(
            aggregated_df,
            occurrences_df,
            on=[self.column_for_merging, "Tumor_Sample_Barcode"],
        )

        # add CNA and SV info
        if os.path.exists(os.path.join(data_folder, "cna.txt")):
            cna_df = pd.read_csv(os.path.join(data_folder, "cna.txt"), sep="\t")
            self.df_collection[study_id]["cna"] = cna_df

            # Exclude 'Hugo_Symbol' column
            columns_to_transform = self.df_collection[study_id][
                "cna"
            ].columns.difference(["Hugo_Symbol"])

            # Apply binary transformation to the selected columns
            df_binary = self.df_collection[study_id]["cna"][columns_to_transform].map(
                lambda x: 1 if (pd.notna(x) and x != 0) else x
            )

            # Add 'Hugo_Symbol' column back to the DataFrame
            df_binary.insert(
                0, "Hugo_Symbol", self.df_collection[study_id]["cna"]["Hugo_Symbol"]
            )

            # Reassign the transformed DataFrame to the collection
            self.df_collection[study_id]["cna_binary"] = df_binary

            melted_cna = self.df_collection[study_id]["cna_binary"].melt(
                id_vars=["Hugo_Symbol"],
                var_name="Tumor_Sample_Barcode",
                value_name="cna_occurrences",
            )

            final_df = pd.merge(
                final_df,
                melted_cna[["Hugo_Symbol", "Tumor_Sample_Barcode", "cna_occurrences"]],
                on=["Hugo_Symbol", "Tumor_Sample_Barcode"],
                how="outer",
            )

        if os.path.exists(os.path.join(data_folder, "sv.txt")):
            sv_df = pd.read_csv(os.path.join(data_folder, "sv.txt"), sep="\t")
            self.df_collection[study_id]["sv"] = sv_df

            # Melt the DataFrame to combine Site1_Hugo_Symbol and Site2_Hugo_Symbol into a single column
            melted_sv = pd.melt(
                self.df_collection[study_id]["sv"],
                id_vars=["Sample_Id"],
                value_vars=["Site1_Hugo_Symbol", "Site2_Hugo_Symbol"],
                var_name="site",
                value_name="Hugo_Symbol",
            )

            # Drop duplicate rows to ensure each Hugo_Symbol is only counted once per Sample_Id
            melted_sv = melted_sv.drop_duplicates(subset=["Sample_Id", "Hugo_Symbol"])

            # Count the occurrences of each Hugo_Symbol in each Sample_Id
            sv_occurrences = (
                melted_sv.groupby(["Hugo_Symbol", "Sample_Id"])
                .size()
                .reset_index(name="sv_occurrences")
            )

            # Rename columns to match the desired output
            sv_occurrences = sv_occurrences.rename(
                columns={"Sample_Id": "Tumor_Sample_Barcode"}
            )

            final_df = pd.merge(
                final_df,
                sv_occurrences[
                    ["Hugo_Symbol", "Tumor_Sample_Barcode", "sv_occurrences"]
                ],
                on=["Hugo_Symbol", "Tumor_Sample_Barcode"],
                how="outer",
            )

        final_df["study_id"] = study_id

        if "Sample Identifier" in sample_df.columns:
            sample_identifier_column = "Sample Identifier"
        elif "#Sample Identifier" in sample_df.columns:
            sample_identifier_column = "#Sample Identifier"
        else:
            raise AssertionError(
                "Sample Identifier column not found in the sample dataframe"
            )

        final_df = pd.merge(
            final_df,
            sample_df[
                [sample_identifier_column, "Cancer Type", "Cancer Type Detailed"]
            ],
            left_on="Tumor_Sample_Barcode",
            right_on=sample_identifier_column,
            how="left",
        )

        final_df.rename(
            columns={
                "Cancer Type": "cancer_type",
                "Cancer Type Detailed": "cancer_type_detailed",
            },
            inplace=True,
        )

        final_df["tissue"] = (
            final_df["cancer_type"]
            .map(CBIO_CANCER_TYPE_TO_TISSUE_DICTIONARY)
            .fillna("unclassified")
        )

        # Drop the redundant SAMPLE_ID column
        final_df.drop(columns=[sample_identifier_column], inplace=True)

        return final_df

    def _create_study_dataframes(self) -> pd.DataFrame:
        self.df_collection = {}

        dataframes = []

        for study_id in self.study_ids:
            self.df_collection[study_id] = {}

            # clean up data just in case
            with open(f"{self.data_dir}/{study_id}/mutations.txt", "r") as file:
                lines = file.readlines()[
                    1:
                ]  # cut out a comment (needed at least for msk_impact_2017)

            if lines[0].split("\t")[0] == "Hugo_Symbol":
                # Write the remaining lines back to the file
                with open(f"{self.data_dir}/{study_id}/mutations.txt", "w") as file:
                    file.writelines(lines)

            final_df = self._create_single_study_dataframe(study_id=study_id)
            dataframes.append(final_df)

        return pd.concat(dataframes, ignore_index=True)

    def _create_or_load_study_dataframes(self) -> pd.DataFrame:
        cache_key = {
            "study_ids": self.study_ids,
            "merge_type": self.merge_type,
            "remove_non_ensembl_genes": self.remove_non_ensembl_genes,
        }
        cache_key = json.dumps(cache_key)
        cache_key = hashlib.sha256(cache_key.encode()).hexdigest()

        cache_dir = os.path.join(self.data_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)

        cache_file = os.path.join(cache_dir, f"{cache_key}.cache.gz")

        if os.path.exists(cache_file):
            if self.verbose:
                logger.info(f"Loading data from cache: {cache_file}")
            try:
                with gzip.open(cache_file, "rt", encoding="UTF-8") as file:
                    cache_data: dict[str, str | dict[str, dict[str, str]]] = json.load(file)

                    df = pd.read_json(cache_data["df"], orient="records")
                    self.df_collection = {
                        k: {
                            k1: pd.read_json(v1, orient="records")
                            for k1, v1 in v.items()
                        }
                        for k, v in cache_data["df_collection"].items()
                    }

                if self.verbose:
                    logger.info(f"Loaded data from cache: {cache_file}")

                return df
            except Exception as e:
                logger.error(f"Error loading cache: {e}, creating new dataframes")

        df = self._create_study_dataframes()

        # write to cache
        if self.verbose:
            logger.info(f"Writing data to cache: {cache_file}")

        try:
            cache_data = {
                "df": df.to_json(orient="records"),
                "df_collection": {
                    k: {
                        k1: v1.to_json(orient="records")
                        for k1, v1 in v.items()
                    }
                    for k, v in self.df_collection.items()
                }
            }
            with gzip.open(cache_file, "wt", encoding="UTF-8") as file:
                json.dump(cache_data, file)

            if self.verbose:
                logger.info(f"Data written to cache: {cache_file}")
        except Exception as e:
            logger.error(f"Error writing cache: {e}")

        return df

    def plot_heatmap(
        self,
        gene_list: list[str],
        stratification: Literal[
            "tissue", "cancer_type", "cancer_type_detailed", "study_id", "sample"
        ] = "tissue",
        filter_category: str | None = None,
        filter_value: str | None = None,
        variation_type: Literal[
            "mutation_occurrences",
            "cna_nonbinary",
            "sv_occurrences",
            "cna_occurrences",
            "Consequence",
        ] = "mutation_occurrences",
    ):
        if variation_type == "cna_nonbinary" or variation_type == "Consequence":
            assert (
                stratification == "sample"
            ), "Stratification must be 'sample' for cna_nonbinary and Consequence variations"

        if variation_type != "cna_nonbinary":
            simple_merge_by_stratification: dict[str, list[str]] = {
                "tissue": ["Tumor_Sample_Barcode", "study_id", "cancer_type", "tissue"],
                "cancer_type": ["Tumor_Sample_Barcode", "study_id", "cancer_type"],
                "cancer_type_detailed": [
                    "Tumor_Sample_Barcode",
                    "study_id",
                    "cancer_type_detailed",
                ],
                "study_id": ["Tumor_Sample_Barcode", "study_id"],
                "sample": ["Tumor_Sample_Barcode"],
            }
            # stratify by merging specific columns
            if stratification in simple_merge_by_stratification:
                merge_on: list[str] = simple_merge_by_stratification[stratification]

                if filter_category is None:  # no filtering
                    final_df = self.big_combined_df
                else:
                    final_df = self.big_combined_df[
                        self.big_combined_df[filter_category] == filter_value
                    ]

                columns_to_keep_copy = self.columns_to_keep.copy()

                unique_samples_info = final_df[
                    [
                        "Tumor_Sample_Barcode",
                        "cancer_type",
                        "cancer_type_detailed",
                        "study_id",
                        "tissue",
                    ]
                ].drop_duplicates()

                hugo_mask = final_df["Hugo_Symbol"].isin(
                    [gene for gene in gene_list if not gene.startswith("ENSG")]
                )

                if self.merge_type == _ENSEMBL:
                    ensg_mask = final_df["Ensembl_Gene_ID"].isin(
                        [gene for gene in gene_list if gene.startswith("ENSG")]
                    )
                    combined_mask = ensg_mask | hugo_mask
                elif self.merge_type == _SYMBOL:
                    combined_mask = hugo_mask
                else:
                    raise AssertionError(f"Invalid merge type: {self.merge_type}")

                filtered_genes_df: pd.DataFrame = final_df[combined_mask]

                if self.merge_type == _ENSEMBL:
                    existing_genes = set(filtered_genes_df["Ensembl_Gene_ID"]).union(
                        set(filtered_genes_df["Hugo_Symbol"])
                    )
                elif self.merge_type == _SYMBOL:
                    existing_genes = set(filtered_genes_df["Hugo_Symbol"])
                else:
                    raise AssertionError(f"Invalid merge type: {self.merge_type}")

                unexpressed_genes = [
                    gene for gene in gene_list if gene not in existing_genes
                ]

                # Get all unique Tumor_Sample_Barcode from the original DataFrame
                all_samples = final_df[merge_on].drop_duplicates()

                all_samples = pd.merge(
                    all_samples, unique_samples_info, on=merge_on, how="left"
                )

                if variation_type not in columns_to_keep_copy:
                    columns_to_keep_copy.append(variation_type)

                must_keep = [
                    "study_id",
                    "tissue",
                    "cancer_type",
                    "cancer_type_detailed",
                ]
                for column_name in must_keep:
                    if (
                        column_name in merge_on
                        and column_name not in columns_to_keep_copy
                    ):
                        columns_to_keep_copy.append(column_name)

                # Merge the filtered genes DataFrame with all samples to ensure all samples are included
                merged_df = pd.merge(
                    all_samples,
                    filtered_genes_df[columns_to_keep_copy],
                    on=merge_on,
                    how="left",
                    suffixes=("_sample", "_gene"),
                )

                if (
                    filter_category is None or stratification != "sample"
                ):  # no filtering
                    df_for_heatmap_very_final: pd.DataFrame = (
                        merged_df.groupby([self.column_for_merging, stratification])[
                            variation_type
                        ]
                        .sum()
                        .reset_index()
                    )
                else:
                    df_for_heatmap_very_final: pd.DataFrame = merged_df

                df = df_for_heatmap_very_final.copy()
            else:
                raise ValueError(
                    "Invalid stratification value. Please choose from 'tissue', 'cancer_type', 'cancer_type_detailed', 'study_id', 'sample'"
                )

            # plot heatmap code
            pivot_column = stratification

            # simple stratification (round 2)
            if stratification in {
                "tissue",
                "cancer_type",
                "cancer_type_detailed",
                "study_id",
            }:
                pass
            elif stratification == "sample":
                pivot_column = "Tumor_Sample_Barcode"
            else:
                raise ValueError(
                    "Invalid stratification value. Please choose from 'tissue', 'cancer_type', 'cancer_type_detailed', 'study_id', 'sample'"
                )

            pivot_df1 = df.pivot(
                index=self.column_for_merging,
                columns=pivot_column,
                values=variation_type,
            )

            pivot_df1 = pivot_df1.dropna(how="all")

            sorted_columns = pivot_df1.count().sort_values(ascending=False).index
            pivot_df1 = pivot_df1[sorted_columns]

            if unexpressed_genes:
                new_rows = pd.DataFrame(
                    np.nan, index=unexpressed_genes, columns=pivot_df1.columns
                )
                pivot_df1 = pd.concat([pivot_df1, new_rows])
            title = f"Heatmap of Gene mutations per gene across {stratification}"

            if filter_value:
                title = title + " - " + filter_value

            if self.column_for_merging == "Ensembl_Gene_ID":
                # todo - add ensembl to hugo conversion
                raise NotImplementedError(
                    "Ensembl to Hugo conversion not implemented yet"
                )
                # pivot_df1.rename(index=self.top_mutant_gene_list_ensembl_to_hugo, inplace=True)

        else:  # variation_type == "cna_nonbinary"
            assert (
                stratification == "sample"
            ), "stratification must be 'sample' for CNA data"
            assert (
                filter_category == "study_id"
            ), "filter_category must be 'study_id' for CNA data"
            pivot_df1 = self.df_collection[filter_value]["cna"].copy()
            pivot_df1.set_index("Hugo_Symbol", inplace=True)
            pivot_df1 = pivot_df1[pivot_df1.index.isin(gene_list)]

            pivot_df1 = pivot_df1.reset_index()
            if "Hugo_Symbol" not in pivot_df1.columns:
                pivot_df1["Hugo_Symbol"] = pivot_df1.index

            # Iterate over the top_mutant_gene_list and add missing entries
            for gene in gene_list:
                if gene not in pivot_df1["Hugo_Symbol"].values:
                    new_row = pd.Series(
                        {col: np.nan for col in pivot_df1.columns}, name=gene
                    )
                    new_row["Hugo_Symbol"] = gene
                    pivot_df1 = pivot_df1.append(new_row)

            # Set 'Hugo_Symbol' back as index if needed
            pivot_df1 = pivot_df1.set_index("Hugo_Symbol")

            title = "Heatmap of CNA data across samples"

        if pivot_df1.isna().all().all():
            print(f"No data to plot for {stratification}")
            return

        nas_present = True

        # ensure pivot_df1 is sorted by columns before plotting
        pivot_df1: pd.DataFrame
        pivot_df1.sort_index(axis="columns", inplace=True)

        # limit to first 500 columns
        render_divider_lines = True
        render_column_ids = pivot_df1.shape[1] < 100
        if (
            pivot_df1.shape[1] > 372
        ):  # 372 is fine, 373 is not. There's something wrong with pyplot...
            print("Warning: Too many columns to plot. Limiting to first 372 columns")
            pivot_df1 = pivot_df1.iloc[:, :372]
            render_divider_lines = False

        if variation_type == "cna_nonbinary":
            min_value = -3
            max_value = 2

            levels = list(range(min_value + 1, max_value + 1))
            pivot_df1 = pivot_df1.fillna(min_value)

            colors_list = plt.get_cmap("RdBu_r", max_value - min_value + 1)(
                range(max_value - min_value + 1)
            )
            colors_list = np.vstack(
                ([[0.5, 0.5, 0.5, 0.3]], colors_list[1:])
            )  # Grey color for -3
            cmap = ListedColormap(colors_list)

            # Define the norm with the diverging palette centered at 0
            norm = TwoSlopeNorm(vmin=min_value, vcenter=0, vmax=max_value)

        elif variation_type == "Consequence":
            consequences: list[str | float] = list(
                self.big_combined_df["Consequence"].unique()
            )

            colors_list = plt.get_cmap("tab20", len(consequences))(
                range(len(consequences))
            )

            # if consequences contains nan, ensure the nan value is at the beginning
            if np.nan in consequences:
                colors_list = np.vstack(([[1.0, 1.0, 1.0, 0.3]], colors_list[:-1]))
                consequences = [np.nan] + sorted(
                    v for v in consequences if not isinstance(v, float)
                )
            else:
                consequences.sort()
                nas_present = False

            cmap = ListedColormap(colors_list)
            min_value = 0
            max_value = len(consequences)
            norm = BoundaryNorm(
                boundaries=np.arange(min_value - 0.5, max_value + 0.5, 1),
                ncolors=cmap.N,
                clip=False,
            )
            levels = list(range(min_value, max_value))

            string_to_int: dict[str | float, int] = {
                consequence: i for i, consequence in enumerate(consequences)
            }

            pivot_df1 = pivot_df1.map(lambda x: string_to_int[x])

        else:
            min_value = int(pivot_df1.min().min())

            if pivot_df1.isna().sum().sum() != 0:
                min_value -= 1
            else:
                nas_present = False

            max_value = max(int(pivot_df1.max().max()), 1)

            levels = list(range(min_value, max_value + 1))
            pivot_df1 = pivot_df1.fillna(min_value)

            # Create a custom colormap
            colors_list = plt.get_cmap("Reds", len(levels))(range(len(levels)))
            if nas_present:
                colors_list = np.vstack(
                    ([[0.5, 0.5, 0.5, 0.3]], colors_list)
                )  # Grey color for -1
            cmap = ListedColormap(colors_list)

            # Define the norm, with vmin set to -1 and vmax to max_value
            norm = BoundaryNorm(
                boundaries=np.arange(min_value - 0.5, max_value + 1.5, 1),
                ncolors=cmap.N,
                clip=False,
            )

        # Create a figure
        plt.figure(figsize=(18, 4))

        # Display the heatmap with the filled array
        plt.imshow(pivot_df1, cmap=cmap, norm=norm, aspect="auto")

        colorbar_label = "Number of Mutations"
        if variation_type == "Consequence":
            colorbar_label = "Consequence Type"
        else:
            # Add colorbar with ticks for all levels, including the placeholder for -1 (NaN)
            levels = _ints_between(min_value, max_value, 25, 7)

        cbar = plt.colorbar(label=colorbar_label, ticks=levels)

        labels: list[str | int] = levels.copy()
        if nas_present:
            labels[0] = "NaN"

        if variation_type == "Consequence":
            # Guaranteed to be bound by code flow
            # noinspection PyUnboundLocalVariable
            # start at 0 if there are no NaN values, otherwise at 1
            for i, consequence in enumerate(consequences[int(nas_present) :]):
                labels[i + int(nas_present)] = consequence

        cbar.ax.set_yticklabels(labels)

        plt.grid(
            which="both",
            axis="both" if render_divider_lines else "y",
            color="black",
            linestyle="-",
            linewidth=0.5,
        )
        if render_column_ids:
            x_labels = pivot_df1.columns
        else:
            x_labels = [""] * len(pivot_df1.columns)
        plt.xticks(np.arange(len(pivot_df1.columns)), x_labels, rotation=90)
        plt.yticks(np.arange(len(pivot_df1.index)), pivot_df1.index)

        plt.gca().set_xticks(np.arange(-0.5, len(pivot_df1.columns), 1), minor=True)
        plt.gca().set_yticks(np.arange(-0.5, len(pivot_df1.index), 1), minor=True)

        plt.grid(which="major", color="white", linestyle="-", linewidth=0.5, alpha=0)
        plt.gca().tick_params(which="minor", size=0)

        plt.xlabel(stratification)
        plt.ylabel("Genes")
        plt.title(title)

        filename = f"Heatmap_{stratification}"
        if filter_value:
            filename = f"{filename}_{filter_value}"

        filepath = os.path.join(self.figure_output_dir, f"{filename}.png")

        plt.savefig(filepath, bbox_inches="tight")

        plt.show()

        plt.close()

        if nas_present:
            pivot_df1 = pivot_df1.replace(min_value, np.nan)

        # if filter_category is None:  # * no filtering
        #     self.pivot_df_dict[variation_type] = pivot_df1
        # else:
        #     self.pivot_df_dict[selected_group_category][selected_group][variation_type] = pivot_df1


def plot(
    study_ids: list[str],
    genes: list[str],
    figure_output_dir: str,

    stratification: Literal[
        "tissue", "cancer_type", "cancer_type_detailed", "study_id", "sample"
    ] = "tissue",
    variation_type: Literal[
        "mutation_occurrences",
        "cna_nonbinary",
        "sv_occurrences",
        "cna_occurrences",
        "Consequence",
    ] = "mutation_occurrences",
    filter_: tuple[str, str] | None = None,

    merge_type: Literal["Ensembl", "Symbol"] = "Symbol",
    remove_non_ensembl_genes: bool = False,

    data_dir: str = "gget_cbio_cache",
    verbose: bool = False,
    confirm_download: bool = False,
) -> bool:
    """
    Plot a heatmap of given genes in the given studies.

    Args:

    :param study_ids:                   list of cBioPortal study IDs
    :param genes:                       list of gene names
    :param figure_output_dir:           directory to save the figure in

    :param stratification:              column to group the data by, default 'tissue'
    :param variation_type:              column to use for the heatmap, default 'mutation_occurrences'
    :param filter_:                     filter to apply to the data in the format (column_name, value), default None

    :param merge_type:                  gene ID type to group by (e.g. 'IL13' or 'ENSG00000169194'), default 'Symbol'
    :param remove_non_ensembl_genes:    whether to remove genes without Ensembl IDs, default False

    :param data_dir:                    directory to store the downloaded data, default 'gget_cbio_cache'
    :param verbose:                     whether to print out progress, default False
    :param confirm_download:            whether to confirm the download before proceeding, default False

    Return:

    :return:                            True if the plot was successful, False otherwise
    """
    if verbose:
        logger.info("Downloading data")

    if not download_cbioportal_data(
        study_ids, verbose=verbose, out_dir=data_dir, confirm_download=confirm_download
    ):
        logger.error("Failed to download data")
        return False

    if verbose:
        logger.info("Loading data")
    gene_analyzer = _GeneAnalysis(
        study_ids,
        genes,
        data_dir=data_dir,
        figure_output_dir=figure_output_dir,
        merge_type=merge_type,
        remove_non_ensembl_genes=remove_non_ensembl_genes,
        verbose=verbose
    )
    if verbose:
        logger.info("Data loaded")
        logger.info("Plotting data")

    gene_analyzer.plot_heatmap(
        genes,
        stratification=stratification,
        filter_category=filter_[0] if filter_ else None,
        filter_value=filter_[1] if filter_ else None,
        variation_type=variation_type,
    )

    del gene_analyzer

    return True
