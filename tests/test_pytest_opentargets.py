# import pytest
# import gget
# import pandas as pd
# import json
# import hashlib

# def test_opentargets():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "diseases",
#         limit = 5
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [
#         [
#             "EFO_0000274",
#             "atopic eczema",
#             "A common chronic pruritic inflammatory skin disease with a strong genetic component. Onset typically occurs during the first 2 years of life.",
#             0.6898147017882402
#         ],
#         [
#             "MONDO_0004979",
#             "asthma",
#             "A bronchial disease that is characterized by chronic inflammation and narrowing of the airways, which is caused by a combination of environmental and genetic factors resulting in recurring periods of wheezing (a whistling sound while breathing), chest tightness, shortness of breath, mucus production and coughing. The symptoms appear due to a variety of triggers such as allergens, irritants, respiratory infections, weather changes, exercise, stress, reflux disease, medications, foods and emotional anxiety.",
#             0.6876614432123751
#         ],
#         [
#             "HP_0000964",
#             "Eczematoid dermatitis",
#             "Eczema is a form of dermatitis that is characterized by scaly, pruritic, erythematous lesions located on flexural surfaces.",
#             0.5960895778362398
#         ],
#         [
#             "EFO_0005854",
#             "allergic rhinitis",
#             "Inflammation of the nasal mucous membranes caused by an IgE-mediated response to external allergens. The inflammation may also involve the mucous membranes of the sinuses, eyes, middle ear, and pharynx. Symptoms include sneezing, nasal congestion, rhinorrhea, and itching. It may lead to fatigue, drowsiness, and malaise thus causing impairment of the quality of life.",
#             0.5352314952446948
#         ],
#         [
#             "EFO_0000676",
#             "psoriasis",
#             "An autoimmune condition characterized by red, well-delineated plaques with silvery scales that are usually on the extensor surfaces and scalp. They can occasionally present with these manifestations: pustules; erythema and scaling in intertriginous areas, and erythroderma, that are often distributed on extensor surfaces and scalp.",
#             0.5244620404842566
#         ]
#     ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"



# def test_opentargets_drugs_filter():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "drugs",
#         filters = {"disease_id": "EFO_0000274"},
#         limit = 2
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [[
#         "CHEMBL1743081",
#         "TRALOKINUMAB",
#         "Antibody",
#         "Interleukin-13 inhibitor",
#         "Antibody drug with a maximum clinical trial phase of IV (across all indications) that was first approved in 2021 and is indicated for atopic eczema and eczematoid dermatitis and has 5 investigational indications.",
#         [
#             "CAT-354",
#             "Tralokinumab"
#         ],
#         [
#             "Adbry",
#             "Adtralza"
#         ],
#         "EFO_0000274",
#         "atopic eczema",
#         4,
#         [],
#         True
#     ],
#     [
#         "CHEMBL1743035",
#         "LEBRIKIZUMAB",
#         "Antibody",
#         "Interleukin-13 inhibitor",
#         "Antibody drug with a maximum clinical trial phase of III (across all indications) and has 4 investigational indications.",
#         ['Lebrikizumab', 'MILR-1444A', 'MILR1444A', 'PRO-301444', 'PRO-301444 RG-3637', 'PRO301444', 'RG-3637'],
#         [],
#         "EFO_0000274",
#         "atopic eczema",
#         3,
#         ['NCT05916365', 'NCT05372419', 'NCT04392154', 'NCT05369403'],
#         False
#     ]]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"

# def test_opentargets_expression():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "expression",
#         limit = 3
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [
#         [
#             "UBERON_0000473",
#             "testis",
#             5,
#             1026.0,
#             "TPM",
#             3,
#             [
#                 "reproductive system"
#             ],
#             [
#                 "reproductive organ",
#                 "reproductive structure"
#             ]
#         ],
#         [
#             "CL_0000542",
#             "EBV-transformed lymphocyte",
#             1,
#             54.0,
#             "",
#             2,
#             [
#                 "hemolymphoid system",
#                 "immune system",
#                 "lymphoid system"
#             ],
#             [
#                 "immune organ"
#             ]
#         ],
#         [
#             "UBERON_0002371",
#             "bone marrow",
#             1,
#             45.5,
#             "",
#             2,
#             [
#                 "hemolymphoid system",
#                 "musculoskeletal system",
#                 "hematopoietic system",
#                 "immune system"
#             ],
#             [
#                 "skeletal element"
#             ]
#         ]
#     ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"

# def test_opentargets_expression_filter():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "expression",
#         filters = {
#             "tissue_id": "UBERON_0000473",
#             "organ": "reproductive organ"
#         },
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [
#         [
#             "UBERON_0000473",
#             "testis",
#             5,
#             1026,
#             "TPM",
#             3,
#             [
#                 "reproductive system"
#             ],
#             [
#                 "reproductive organ",
#                 "reproductive structure"
#             ]
#         ]
#     ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"

# def test_opentargets_interactions():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "interactions",
#         limit = 2
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [
#         [
#             0.999,
#             3,
#             "string",
#             "ENSP00000304915",
#             "ENSG00000169194",
#             "IL13",
#             "unspecified role",
#             9606,
#             "ENSP00000361004",
#             "ENSG00000123496",
#             "IL13RA2",
#             "unspecified role",
#             9606
#         ],
#         [
#             0.999,
#             3,
#             "string",
#             "ENSP00000304915",
#             "ENSG00000169194",
#             "IL13",
#             "unspecified role",
#             9606,
#             "ENSP00000360730",
#             "ENSG00000131724",
#             "IL13RA1",
#             "unspecified role",
#             9606
#         ]
#     ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"


# def test_opentargets_no_specified_resource():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         limit = 5
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()

#     expected_result = [
#             [
#                 "EFO_0000274",
#                 "atopic eczema",
#                 "A common chronic pruritic inflammatory skin disease with a strong genetic component. Onset typically occurs during the first 2 years of life.",
#                 0.6993140228020157
#             ],
#             [
#                 "MONDO_0004979",
#                 "asthma",
#                 "A bronchial disease that is characterized by chronic inflammation and narrowing of the airways, which is caused by a combination of environmental and genetic factors resulting in recurring periods of wheezing (a whistling sound while breathing), chest tightness, shortness of breath, mucus production and coughing. The symptoms appear due to a variety of triggers such as allergens, irritants, respiratory infections, weather changes, exercise, stress, reflux disease, medications, foods and emotional anxiety.",
#                 0.6780711161785388
#             ],
#             [
#                 "HP_0000964",
#                 "Eczematoid dermatitis",
#                 "Eczema is a form of dermatitis that is characterized by scaly, pruritic, erythematous lesions located on flexural surfaces.",
#                 0.6148607200185839
#             ],
#             [
#                 "EFO_0005854",
#                 "allergic rhinitis",
#                 "Inflammation of the nasal mucous membranes caused by an IgE-mediated response to external allergens. The inflammation may also involve the mucous membranes of the sinuses, eyes, middle ear, and pharynx. Symptoms include sneezing, nasal congestion, rhinorrhea, and itching. It may lead to fatigue, drowsiness, and malaise thus causing impairment of the quality of life.",
#                 0.5666922121972621
#             ],
#             [
#                 "EFO_0000676",
#                 "psoriasis",
#                 "An autoimmune condition characterized by red, well-delineated plaques with silvery scales that are usually on the extensor surfaces and scalp. They can occasionally present with these manifestations: pustules; erythema and scaling in intertriginous areas, and erythroderma, that are often distributed on extensor surfaces and scalp.",
#                 0.5492499302524128
#             ]
#         ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"

# def test_opentargets_pharmacogenetics():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "pharmacogenetics",
#         limit = 1
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()
    
#     # If any value of result is a 1-row dataframe, then convert to dict
#     for h, sublist in enumerate(result):
#         for i, item in enumerate(sublist):
#             if isinstance(item, pd.DataFrame) and len(item) == 1:
#                 result[h][i] = item.iloc[0].to_dict()

#     expected_result = [[
#         'rs1800925',
#         '5_132657117_C_T,T',
#         'TT',
#         'SO:0001631',
#         'upstream_gene_variant',
#         {
#             "id": "CHEMBL535",
#             "name": "sunitinib"
#         },
#         'decreased severity of drug-induced toxicity',
#         'Patients with renal cell carcinoma and the TT genotype may have a decreased severity of drug-induced toxicity when administered sunitinib as compared to patients with the CC or CT genotypes. Other clinical and genetic factors may also influence severity of drug-induced toxicity in patients with renal cell carcinoma who are administered sunitinib.',
#         'toxicity',
#         False,
#         '3',
#         'pharmgkb',
#         [
#             '26387812', '26387812'
#         ],
#     ]]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"


# def test_opentargets_pharmacogenetics_json():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "pharmacogenetics",
#         limit = 1,
#         json = True
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()
    
#     # If any value of result is a 1-row dataframe, then convert to dict
#     for h, sublist in enumerate(result):
#         for i, item in enumerate(sublist):
#             if isinstance(item, pd.DataFrame) and len(item) == 1:
#                 result[h][i] = item.iloc[0].to_dict()

#     expected_result = [{
#                 "rs_id": "rs1800925",
#                 "genotype_id": "5_132657117_C_T,T",
#                 "genotype": "TT",
#                 "variant_consequence_id": "SO:0001631",
#                 "variant_consequence_label": "upstream_gene_variant",
#                 "drugs": [
#                     {
#                         "id": "CHEMBL535",
#                         "name": "sunitinib"
#                     }
#                 ],
#                 "phenotype": "decreased severity of drug-induced toxicity",
#                 "genotype_annotation": "Patients with renal cell carcinoma and the TT genotype may have a decreased severity of drug-induced toxicity when administered sunitinib as compared to patients with the CC or CT genotypes. Other clinical and genetic factors may also influence severity of drug-induced toxicity in patients with renal cell carcinoma who are administered sunitinib.",
#                 "response_category": "toxicity",
#                 "direct_target": False,
#                 "evidence_level": "3",
#                 "source": "pharmgkb",
#                 "literature": [
#                     "26387812", "26387812"
#                 ]
#             }]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"


# def test_opentargets_interactions_complex_filter():
#     result = gget.opentargets(
#         ensembl_id = "ENSG00000169194",
#         resource = "interactions",
#         filters = {
#             "protein_a_id": "P35225",
#             "gene_b_id": [
#                 "ENSG00000077238",
#                 "ENSG00000111537"
#             ]
#         },
#         filter_mode = "or",
#         limit = 5
#     )

#     if isinstance(result, pd.DataFrame):
#         result = result.dropna(axis=1).values.tolist()
    
#     # If any value of result is a 1-row dataframe, then convert to dict
#     for h, sublist in enumerate(result):
#         for i, item in enumerate(sublist):
#             if isinstance(item, pd.DataFrame) and len(item) == 1:
#                 result[h][i] = item.iloc[0].to_dict()

#     expected_result = [
#             [
#                 0.999,
#                 3,
#                 "string",
#                 "ENSP00000304915",
#                 "ENSG00000169194",
#                 "IL13",
#                 "unspecified role",
#                 9606,
#                 "ENSP00000379111",
#                 "ENSG00000077238",
#                 "IL4R",
#                 "unspecified role",
#                 9606
#             ],
#             [
#                 0.961,
#                 2,
#                 "string",
#                 "ENSP00000304915",
#                 "ENSG00000169194",
#                 "IL13",
#                 "unspecified role",
#                 9606,
#                 "ENSP00000229135",
#                 "ENSG00000111537",
#                 "IFNG",
#                 "unspecified role",
#                 9606
#             ],
#             [
#                 0.8,
#                 7,
#                 "intact",
#                 "P35225",
#                 "ENSG00000169194",
#                 "IL13",
#                 "unspecified role",
#                 9606,
#                 "P78552",
#                 "ENSG00000131724",
#                 "IL13RA1",
#                 "unspecified role",
#                 9606
#             ],
#             [
#                 0.8,
#                 9,
#                 "intact",
#                 "P35225",
#                 "ENSG00000169194",
#                 "IL13",
#                 "unspecified role",
#                 9606,
#                 "Q14627",
#                 "ENSG00000123496",
#                 "IL13RA2",
#                 "unspecified role",
#                 9606
#             ],
#             [
#                 0.53,
#                 2,
#                 "intact",
#                 "P35225",
#                 "ENSG00000169194",
#                 "IL13",
#                 "unspecified role",
#                 9606,
#                 "Q86YS7",
#                 "ENSG00000111731",
#                 "C2CD5",
#                 "unspecified role",
#                 9606
#             ]
#         ]

#     assert result == expected_result, f"Expected {expected_result}, but got {result}"


# # def test_opentargets_drugs_no_limit():
# #     result = gget.opentargets(
# #         ensembl_id = "ENSG00000169194",
# #         resource = "drugs",
# #     )

# #     if isinstance(result, pd.DataFrame):
# #         result = result.dropna(axis=1).values.tolist()

# #     result = json.dumps(result)
# #     result = hashlib.md5(result.encode()).hexdigest()

# #     assert result == '872ca75cba98a2383e4fd116682ffe2b'
