import unittest
from gget.gget_info import info


class TestInfo(unittest.TestCase):
    maxDiff = None

    def test_info_WB_gene(self):
        df = info("WBGene00043981")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "WBGene00043981",
                "Q5WRS0",
                "3565421",
                "caenorhabditis_elegans",
                "WBcel235",
                "aaim-1",
                "T14E8.4",
                ["T14E8.4", "aaim-1"],
                "Protein aaim-1 (Antibacterial and aids invasion by microsporidia 1 protein)",
                "Uncharacterized protein [Source:NCBI gene (formerly Entrezgene);Acc:3565421]",
                "FUNCTION: Plays a role in promoting resistance to bacterial pathogens such as P.aeruginosa by inhibiting bacterial intestinal colonization. {ECO:0000269|PubMed:34994689}.; FUNCTION: (Microbial infection) Promotes infection by microsporidian pathogens such as N.parisii in the early larval stages of development (PubMed:34994689). Involved in ensuring the proper orientation and location of the spore proteins of N.parisii during intestinal cell invasion (PubMed:34994689). {ECO:0000269|PubMed:34994689}.",
                "Gene",
                "protein_coding",
                "T14E8.4.1.",
                "X",
                -1,
                6559466,
                6562428,
                ["T14E8.4.1"],
                ["protein_coding"],
                [unittest.mock.ANY],
            ]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_info_WB_transcript(self):
        df = info("T14E8.4.1")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "T14E8.4.1",
                "Q5WRS0",
                "caenorhabditis_elegans",
                "WBcel235",
                "aaim-1",
                ["T14E8.4", "aaim-1"],
                "WBGene00043981",
                "Protein aaim-1 (Antibacterial and aids invasion by microsporidia 1 protein)",
                "FUNCTION: Plays a role in promoting resistance to bacterial pathogens such as P.aeruginosa by inhibiting bacterial intestinal colonization. {ECO:0000269|PubMed:34994689}.; FUNCTION: (Microbial infection) Promotes infection by microsporidian pathogens such as N.parisii in the early larval stages of development (PubMed:34994689). Involved in ensuring the proper orientation and location of the spore proteins of N.parisii during intestinal cell invasion (PubMed:34994689). {ECO:0000269|PubMed:34994689}.",
                "Transcript",
                "protein_coding",
                "X",
                -1,
                6559466,
                6562428,
                [
                    "T14E8.4.1.e1",
                    "T14E8.4.1.e2",
                    "T14E8.4.1.e3",
                    "T14E8.4.1.e4",
                    "T14E8.4.1.e5",
                    "T14E8.4.1.e6",
                    "T14E8.4.1.e7",
                    "T14E8.4.1.e8",
                    "T14E8.4.1.e9",
                ],
                [
                    6562330,
                    6562225,
                    6562110,
                    6561817,
                    6561252,
                    6560727,
                    6560492,
                    6560197,
                    6559466,
                ],
                [
                    6562428,
                    6562286,
                    6562183,
                    6562059,
                    6561378,
                    6561206,
                    6560671,
                    6560443,
                    6559667,
                ],
                [
                    "T14E8.4.1.e1",
                    "T14E8.4.1.e2",
                    "T14E8.4.1.e3",
                    "T14E8.4.1.e4",
                    "T14E8.4.1.e5",
                    "T14E8.4.1.e6",
                    "T14E8.4.1.e7",
                    "T14E8.4.1.e8",
                    "T14E8.4.1.e9",
                ],
                [
                    6562330,
                    6562225,
                    6562110,
                    6561817,
                    6561252,
                    6560727,
                    6560492,
                    6560197,
                    6559466,
                ],
                [
                    6562428,
                    6562286,
                    6562183,
                    6562059,
                    6561378,
                    6561206,
                    6560671,
                    6560443,
                    6559667,
                ],
            ]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_info_FB_gene(self):
        df = info("FBgn0003656")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "FBgn0003656",
                "Q9U969",
                "31716",
                "drosophila_melanogaster",
                "BDGP6.32",
                "sws",
                "sws",
                ["CG2212", "Dmel\\CG2212", "PNPLA6", "SWS", "Sws", "olfE", "sws"],
                "Neuropathy target esterase sws (Swiss cheese) (DSWS) (EC 3.1.1.5)",
                "swiss cheese",
                "FUNCTION: Phospholipase B that deacylates intracellular phosphatidylcholine (PtdCho), generating glycerophosphocholine (GroPtdCho). This deacylation occurs at both sn-2 and sn-1 positions of PtdCho. Its specific chemical modification by certain organophosphorus (OP) compounds leads to distal axonopathy. Plays a role in the signaling mechanism between neurons and glia that regulates glia wrapping during development of the adult brain. Essential for membrane lipid homeostasis and cell survival in both neurons and glia of the adult brain. {ECO:0000269|PubMed:15772346, ECO:0000269|PubMed:18945896, ECO:0000269|PubMed:9295388}.",
                "Enables lysophospholipase activity and protein kinase A catalytic subunit binding activity. Involved in several processes, including negative regulation of cAMP-dependent protein kinase activity; photoreceptor cell maintenance; and sensory perception of smell. Located in endoplasmic reticulum membrane and plasma membrane. Is expressed in adult head and interface glial cell. Used to study blindness; cerebellar ataxia; hereditary spastic paraplegia; and neurodegenerative disease. Human ortholog(s) of this gene implicated in Boucher-Neuhauser syndrome; Laurence-Moon syndrome; Oliver-McFarlane syndrome; and hereditary spastic paraplegia 39. Orthologous to human PNPLA6 (patatin like phospholipase domain containing 6) and PNPLA7 (patatin like phospholipase domain containing 7). [provided by Alliance of Genome Resources, Apr 2022]",
                "Gene",
                "protein_coding",
                "FBtr0071125.",
                "X",
                -1,
                7956820,
                7968236,
                ["FBtr0301675", "FBtr0071125", "FBtr0071126"],
                ["protein_coding", "protein_coding", "protein_coding"],
                ["sws-RC", "sws-RA", "sws-RB"],
            ]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene(self):
        df = info("ENSMUSG00000000001")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUSG00000000001.5",
                "Q9DC51",
                "14679",
                "mus_musculus",
                "GRCm39",
                "Gnai3",
                "Gnai3",
                ["AI158965", "AW537698", "Galphai3", "Gnai-3", "Gnai3", "Hg1a"],
                "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)",
                "guanine nucleotide binding protein (G protein), alpha inhibiting 3 [Source:MGI Symbol;Acc:MGI:95773]",
                "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}.",
                "Predicted to enable several functions, including G-protein beta/gamma-subunit complex binding activity; GDP binding activity; and GTPase activating protein binding activity. Predicted to be involved in several processes, including positive regulation of NAD(P)H oxidase activity; positive regulation of superoxide anion generation; and positive regulation of vascular associated smooth muscle cell proliferation. Predicted to act upstream of or within G protein-coupled receptor signaling pathway. Located in Golgi apparatus. Is expressed in early conceptus; inner ear; and oocyte. Orthologous to human GNAI3 (G protein subunit alpha i3). [provided by Alliance of Genome Resources, Apr 2022]",
                "Gene",
                "protein_coding",
                "ENSMUST00000000001.5",
                "3",
                -1,
                108014596,
                108053462,
                ["ENSMUST00000000001"],
                ["protein_coding"],
                ["Gnai3-201"],
            ]
        ]
        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene_list_non_model(self):
        df = info(
            ["ENSMMUG00000054106.1", "ENSMMUG00000053116.1", "ENSMMUG00000021246.4"]
        )
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMMUG00000054106.1",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "lncRNA",
                "ENSMMUT00000080640.1",
                "8",
                1,
                64990191,
                65000159,
                ["ENSMMUT00000080640", "ENSMMUT00000100253"],
                ["lncRNA", "lncRNA"],
                [unittest.mock.ANY, unittest.mock.ANY],
            ],
            [
                "ENSMMUG00000053116.1",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000091015.1",
                "3",
                -1,
                111461994,
                111475279,
                ["ENSMMUT00000091015"],
                ["protein_coding"],
                [unittest.mock.ANY],
            ],
            [
                "ENSMMUG00000021246.4",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000029894.4",
                "2",
                -1,
                98646979,
                98755023,
                [
                    "ENSMMUT00000029893",
                    "ENSMMUT00000053619",
                    "ENSMMUT00000104418",
                    "ENSMMUT00000087615",
                    "ENSMMUT00000103912",
                    "ENSMMUT00000086824",
                    "ENSMMUT00000029894",
                    "ENSMMUT00000104481",
                    "ENSMMUT00000090481",
                    "ENSMMUT00000026408",
                ],
                [
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                ],
                [
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                    unittest.mock.ANY,
                ],
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene_list_non_model_json(self):
        result_to_test = info(
            ["ENSMMUG00000054106.1", "ENSMMUG00000053116.1", "ENSMMUT00000091015.1"],
            json=True,
        )
        expected_result = {
            "ENSMMUG00000054106": {
                "ensembl_id": "ENSMMUG00000054106.1",
                "uniprot_id": None,
                "ncbi_gene_id": None,
                "species": "macaca_mulatta",
                "assembly_name": "Mmul_10",
                "primary_gene_name": None,
                "ensembl_gene_name": None,
                "synonyms": None,
                "parent_gene": None,
                "protein_names": None,
                "ensembl_description": None,
                "uniprot_description": None,
                "ncbi_description": None,
                "object_type": "Gene",
                "biotype": "lncRNA",
                "canonical_transcript": "ENSMMUT00000080640.1",
                "seq_region_name": "8",
                "strand": 1,
                "start": 64990191,
                "end": 65000159,
                "all_transcripts": [
                    {
                        "transcript_id": "ENSMMUT00000080640",
                        "transcript_biotype": "lncRNA",
                        "transcript_name": None,
                    },
                    {
                        "transcript_id": "ENSMMUT00000100253",
                        "transcript_biotype": "lncRNA",
                        "transcript_name": None,
                    },
                ],
                "all_exons": [],
                "all_translations": [],
            },
            "ENSMMUG00000053116": {
                "ensembl_id": "ENSMMUG00000053116.1",
                "uniprot_id": "A0A5F8AEA0",
                "ncbi_gene_id": None,
                "species": "macaca_mulatta",
                "assembly_name": "Mmul_10",
                "primary_gene_name": None,
                "ensembl_gene_name": None,
                "synonyms": None,
                "parent_gene": None,
                "protein_names": "Uncharacterized protein",
                "ensembl_description": None,
                "uniprot_description": None,
                "ncbi_description": None,
                "object_type": "Gene",
                "biotype": "protein_coding",
                "canonical_transcript": "ENSMMUT00000091015.1",
                "seq_region_name": "3",
                "strand": -1,
                "start": 111461994,
                "end": 111475279,
                "all_transcripts": [
                    {
                        "transcript_id": "ENSMMUT00000091015",
                        "transcript_biotype": "protein_coding",
                        "transcript_name": None,
                    }
                ],
                "all_exons": [],
                "all_translations": [],
            },
            "ENSMMUT00000091015": {
                "ensembl_id": "ENSMMUT00000091015.1",
                "uniprot_id": "A0A5F8AEA0",
                "ncbi_gene_id": None,
                "species": "macaca_mulatta",
                "assembly_name": "Mmul_10",
                "primary_gene_name": None,
                "ensembl_gene_name": None,
                "synonyms": None,
                "parent_gene": "ENSMMUG00000053116",
                "protein_names": "Uncharacterized protein",
                "ensembl_description": None,
                "uniprot_description": None,
                "ncbi_description": None,
                "object_type": "Transcript",
                "biotype": "protein_coding",
                "canonical_transcript": None,
                "seq_region_name": "3",
                "strand": -1,
                "start": 111461994,
                "end": 111475279,
                "all_transcripts": [],
                "all_exons": [
                    {
                        "exon_id": "ENSMMUE00000479945",
                        "exon_start": 111475069,
                        "exon_end": 111475279,
                    },
                    {
                        "exon_id": "ENSMMUE00000514646",
                        "exon_start": 111461994,
                        "exon_end": 111468349,
                    },
                ],
                "all_translations": [
                    {
                        "translation_id": "ENSMMUE00000479945",
                        "translation_start": 111475069,
                        "translation_end": 111475279,
                    },
                    {
                        "translation_id": "ENSMMUE00000514646",
                        "translation_start": 111461994,
                        "translation_end": 111468349,
                    },
                ],
            },
        }

        self.assertEqual(result_to_test, expected_result)

    def test_info_transcript(self):
        df = info("ENSMUST00000000001.1")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUST00000000001.5",
                "Q9DC51",
                "14679",
                "mus_musculus",
                "GRCm39",
                "Gnai3",
                "Gnai3-201",
                ["AI158965", "AW537698", "Galphai3", "Gnai-3", "Gnai3", "Hg1a"],
                "ENSMUSG00000000001",
                "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)",
                "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}.",
                "Predicted to enable several functions, including G-protein beta/gamma-subunit complex binding activity; GDP binding activity; and GTPase activating protein binding activity. Predicted to be involved in several processes, including positive regulation of NAD(P)H oxidase activity; positive regulation of superoxide anion generation; and positive regulation of vascular associated smooth muscle cell proliferation. Predicted to act upstream of or within G protein-coupled receptor signaling pathway. Located in Golgi apparatus. Is expressed in early conceptus; inner ear; and oocyte. Orthologous to human GNAI3 (G protein subunit alpha i3). [provided by Alliance of Genome Resources, Apr 2022]",
                "Transcript",
                "protein_coding",
                "3",
                -1,
                108014596,
                108053462,
                [
                    "ENSMUSE00000334714",
                    "ENSMUSE00000276500",
                    "ENSMUSE00000276490",
                    "ENSMUSE00000276482",
                    "ENSMUSE00000565003",
                    "ENSMUSE00000565001",
                    "ENSMUSE00000565000",
                    "ENSMUSE00000404895",
                    "ENSMUSE00000363317",
                ],
                [
                    108053204,
                    108031111,
                    108030858,
                    108025617,
                    108023079,
                    108019789,
                    108019251,
                    108016719,
                    108014596,
                ],
                [
                    108053462,
                    108031153,
                    108030999,
                    108025774,
                    108023207,
                    108019918,
                    108019404,
                    108016928,
                    108016632,
                ],
                [
                    "ENSMUSE00000334714",
                    "ENSMUSE00000276500",
                    "ENSMUSE00000276490",
                    "ENSMUSE00000276482",
                    "ENSMUSE00000565003",
                    "ENSMUSE00000565001",
                    "ENSMUSE00000565000",
                    "ENSMUSE00000404895",
                    "ENSMUSE00000363317",
                ],
                [
                    108053204,
                    108031111,
                    108030858,
                    108025617,
                    108023079,
                    108019789,
                    108019251,
                    108016719,
                    108014596,
                ],
                [
                    108053462,
                    108031153,
                    108030999,
                    108025774,
                    108023207,
                    108019918,
                    108019404,
                    108016928,
                    108016632,
                ],
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_mix(self):
        df = info(["ENSTGUT00000027003.1", "ENSG00000169174"])
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSTGUT00000027003.1",
                "A0A674GVD2",
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "FUNDC1",
                "FUNDC1-202",
                ["FUNDC1"],
                "Uncharacterized protein",
                "Transcript",
                "protein_coding",
                "1",
                -1,
                107513786,
                107526965,
            ],
            [
                "ENSG00000169174.11",
                "Q8NBP7",
                "homo_sapiens",
                "GRCh38",
                "PCSK9",
                "PCSK9",
                [
                    "FH3",
                    "FHCL3",
                    "HCHOLA3",
                    "LDLCQ1",
                    "NARC-1",
                    "NARC1",
                    "PC9",
                    "PCSK9",
                    "PSEC0052",
                ],
                "Proprotein convertase subtilisin/kexin type 9 (EC 3.4.21.-) (Neural apoptosis-regulated convertase 1) (NARC-1) (Proprotein convertase 9) (PC9) (Subtilisin/kexin-like protease PC9)",
                "Gene",
                "protein_coding",
                "1",
                1,
                55039447,
                55064852,
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_exon_expand(self):
        df = info("ENSTGUEE00000179311")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSTGUEE00000179311.1",
                "taeniopygia_guttata",
                "bTaeGut1_v1.p",
                "Exon",
                "1",
                -1,
                107526792,
                107526965,
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_bad_id(self):
        result = info(["banana"])
        self.assertIsNone(result, "Invalid ID output is not None.")
