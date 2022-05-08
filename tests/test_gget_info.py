import unittest
from gget.gget_info import info


class TestInfo(unittest.TestCase):
    def test_info_gene(self):
        df = info("ENSMUSG00000000001")
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                ["Q9DC51"],
                "14679",
                "mus_musculus",
                "GRCm39",
                ["Gnai3"],
                "Gnai3",
                ([["Gnai-3", "Hg1a", "Gnai3", "AW537698", "AI158965", "Galphai3"]],),
                [
                    "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)"
                ],
                "guanine nucleotide binding protein (G protein), alpha inhibiting 3 [Source:MGI Symbol;Acc:MGI:95773]",
                [
                    "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}."
                ],
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

    def test_info_gene_expand(self):
        df = info("ENSMUSG00000000001", expand=True)
        # Drop NaN columns, since np.nan != np.nan
        result_to_test = df.dropna(axis=1).values.tolist()
        expected_result = [
            [
                ["Q9DC51"],
                "14679",
                "mus_musculus",
                "GRCm39",
                ["Gnai3"],
                "Gnai3",
                ([["AW537698", "Gnai3", "AI158965", "Galphai3", "Hg1a", "Gnai-3"]],),
                [
                    "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)"
                ],
                "guanine nucleotide binding protein (G protein), alpha inhibiting 3 [Source:MGI Symbol;Acc:MGI:95773]",
                [
                    "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}."
                ],
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
                ["Q9DC51"],
                "14679",
                "mus_musculus",
                "GRCm39",
                ["Gnai3"],
                "Gnai3-201",
                ([["AW537698", "Gnai3", "AI158965", "Galphai3", "Hg1a", "Gnai-3"]],),
                "ENSMUSG00000000001",
                [
                    "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)"
                ],
                [
                    "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}."
                ],
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
                ["Q9DC51"],
                "14679",
                "mus_musculus",
                "GRCm39",
                ["Gnai3"],
                "Gnai3-201",
                ([["AW537698", "Gnai3", "AI158965", "Galphai3", "Hg1a", "Gnai-3"]],),
                "ENSMUSG00000000001",
                [
                    "Guanine nucleotide-binding protein G(i) subunit alpha-3 (G(i) alpha-3)"
                ],
                [
                    "FUNCTION: Heterotrimeric guanine nucleotide-binding proteins (G proteins) function as transducers downstream of G protein-coupled receptors (GPCRs) in numerous signaling cascades. The alpha chain contains the guanine nucleotide binding site and alternates between an active, GTP-bound state and an inactive, GDP-bound state. Signaling by an activated GPCR promotes GDP release and GTP binding. The alpha subunit has a low GTPase activity that converts bound GTP to GDP, thereby terminating the signal. Both GDP release and GTP hydrolysis are modulated by numerous regulatory proteins. Signaling is mediated via effector proteins, such as adenylate cyclase. Inhibits adenylate cyclase activity, leading to decreased intracellular cAMP levels. Stimulates the activity of receptor-regulated K(+) channels. The active GTP-bound form prevents the association of RGS14 with centrosomes and is required for the translocation of RGS14 from the cytoplasm to the plasma membrane. May play a role in cell division. {ECO:0000250|UniProtKB:P08754}."
                ],
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
                "Q9DC51",
                "mus_musculus",
                "GRCm39",
                "Gnai3",
                "Gnai3",
                ([["Gnai3", "Gnai-3", "AW537698", "AI158965", "Hg1a", "Galphai3"]],),
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
