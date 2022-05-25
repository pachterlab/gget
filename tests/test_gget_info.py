import unittest
from gget.gget_info import info


class TestInfo(unittest.TestCase):
    maxDiff = None

    def test_info_gene(self):
        df = info("ENSMUSG00000000001")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUSG00000000001",
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
                "ENSMMUG00000054106",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "lncRNA",
                "ENSMMUT00000080640.1",
                "8",
                1,
                64990191,
                65000159,
            ],
            [
                "ENSMMUG00000053116",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000091015.1",
                "3",
                -1,
                111461994,
                111475279,
            ],
            [
                "ENSMMUG00000021246",
                "macaca_mulatta",
                "Mmul_10",
                "Gene",
                "protein_coding",
                "ENSMMUT00000029894.4",
                "2",
                -1,
                98646979,
                98755023,
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_gene_list_non_model_json(self):
        result_to_test = info(
            ["ENSMMUG00000054106.1", "ENSMMUG00000053116.1", "ENSMMUG00000021246.4"],
            json=True,
        )
        expected_result = {
            "ENSMMUG00000054106": {
                "ensembl_id": "ENSMMUG00000054106",
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
            },
            "ENSMMUG00000053116": {
                "ensembl_id": "ENSMMUG00000053116",
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
            },
            "ENSMMUG00000021246": {
                "ensembl_id": "ENSMMUG00000021246",
                "uniprot_id": [
                    "A0A1D5QWP5",
                    "A0A5F7ZUG9",
                    "A0A5F7ZY65",
                    "A0A5F7ZZI0",
                    "A0A5F8AMK9",
                    "F7D2F4",
                    "F7HRJ1",
                    "F7HRJ3",
                    "G7MIX6",
                ],
                "ncbi_gene_id": None,
                "species": "macaca_mulatta",
                "assembly_name": "Mmul_10",
                "primary_gene_name": "HHATL",
                "ensembl_gene_name": None,
                "synonyms": ["EGK_11753", "HHATL", "HIGD1A"],
                "parent_gene": None,
                "protein_names": [
                    "HIG1 domain family member 1A isoform a",
                    "HIG1 domain-containing protein",
                    "Uncharacterized protein",
                ],
                "ensembl_description": None,
                "uniprot_description": [
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
                "ncbi_description": None,
                "object_type": "Gene",
                "biotype": "protein_coding",
                "canonical_transcript": "ENSMMUT00000029894.4",
                "seq_region_name": "2",
                "strand": -1,
                "start": 98646979,
                "end": 98755023,
            },
        }

        self.assertEqual(result_to_test, expected_result)

    def test_info_gene_expand(self):
        df = info("ENSMUSG00000000001", expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUSG00000000001",
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

    def test_info_transcript(self):
        df = info("ENSMUST00000000001.1")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUST00000000001",
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
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_transcript_expand(self):
        df = info("ENSMUST00000000001.1", expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSMUST00000000001",
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
        df = info(["ENSTGUT00000027003.1", "ENSMUSG00000000001"])
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSTGUT00000027003",
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
                "ENSMUSG00000000001",
                "Q9DC51",
                "mus_musculus",
                "GRCm39",
                "Gnai3",
                "Gnai3",
                ["AI158965", "AW537698", "Galphai3", "Gnai-3", "Gnai3", "Hg1a"],
                "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)",
                "Gene",
                "protein_coding",
                "3",
                -1,
                108014596,
                108053462,
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_info_exon(self):
        df = info(["ENSTGUEE00000179311"])
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSTGUEE00000179311",
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

    def test_info_exon_expand(self):
        df = info(["ENSTGUEE00000179311"], expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                "ENSTGUEE00000179311",
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

    def test_info_bad_id_expand(self):
        result = info(["banana"], expand=True)
        self.assertIsNone(result, "Invalid ID output is not None.")
