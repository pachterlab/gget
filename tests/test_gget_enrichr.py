import unittest
import matplotlib
import matplotlib.pyplot as plt

# Prevent matplotlib from opening windows
matplotlib.use("Agg")
from gget.gget_enrichr import enrichr


class TestEnrichr(unittest.TestCase):
    def test_enrichr_pathway(self):
        df = enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="pathway",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            [
                1,
                "One carbon pool by folate",
                0.008965773496317424,
                131.39473684210526,
                619.439581424405,
                ["MTHFD2L"],
                0.01793154699263485,
                "KEGG_2021_Human",
            ],
            [
                2,
                "Mitophagy",
                0.030192894641305958,
                37.17164179104478,
                130.10627217002806,
                ["FUNDC1"],
                0.030192894641305958,
                "KEGG_2021_Human",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_none(self):
        df = enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                None,
                "ASF1A",
                None,
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="pathway",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            [
                1,
                "One carbon pool by folate",
                0.008965773496317424,
                131.39473684210526,
                619.439581424405,
                ["MTHFD2L"],
                0.01793154699263485,
                "KEGG_2021_Human",
            ],
            [
                2,
                "Mitophagy",
                0.030192894641305958,
                37.17164179104478,
                130.10627217002806,
                ["FUNDC1"],
                0.030192894641305958,
                "KEGG_2021_Human",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_transcription(self):
        df = enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="transcription",
        )
        result_to_test = df.values.tolist()[:5]
        expected_result = [
            [
                1,
                "PRDM5 23873026 ChIP-Seq MEFs Mouse",
                0.0007126827430407046,
                14.802731707317074,
                107.26761336583695,
                ["AIMP2", "FUNDC2", "MFHAS1", "ASF1A"],
                0.1860101959336239,
                "ChEA_2016",
            ],
            [
                2,
                "DMRT1 21621532 ChIP-ChIP FETAL Ovary",
                0.0015097917368846782,
                43.65054945054945,
                283.5445215172139,
                ["AIMP2", "FUNDC2"],
                0.1970278216634505,
                "ChEA_2016",
            ],
            [
                3,
                "FLI1 21571218 ChIP-Seq MEGAKARYOCYTES Human",
                0.003587828117602931,
                8.507636862879698,
                47.89976717448956,
                ["AIMP1", "TRMT112", "FUNDC2", "MFHAS1", "MTHFD2L", "BFAR", "ASF1A"],
                0.312141046231455,
                "ChEA_2016",
            ],
            [
                4,
                "CTNNB1 20460455 ChIP-Seq HCT116 Human",
                0.008067670326181822,
                9.64771573604061,
                46.50093363148564,
                ["AIMP2", "MFHAS1", "MTHFD2L"],
                0.5264154887833639,
                "ChEA_2016",
            ],
            [
                5,
                "RUNX2 22187159 ChIP-Seq PCA Human",
                0.05229842506931902,
                3.877625036560398,
                11.44205338839197,
                ["AIMP1", "MFHAS1", "FUNDC1", "MTHFD2L"],
                0.6933970446278469,
                "ChEA_2016",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_ontology(self):
        df = enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="ontology",
        )
        result_to_test = df.values.tolist()[:5]
        expected_result = [
            [
                1,
                "translation (GO:0006412)",
                9.675824801901589e-05,
                46.872037914691944,
                433.25207270512743,
                ["AIMP1", "AIMP2", "TRMT112"],
                0.00328427763325889,
                "GO_Biological_Process_2021",
            ],
            [
                2,
                "mitochondrion disassembly (GO:0061726)",
                0.00011250476197410829,
                167.7058823529412,
                1524.868252317394,
                ["FUNDC2", "FUNDC1"],
                0.00328427763325889,
                "GO_Biological_Process_2021",
            ],
            [
                3,
                "tRNA aminoacylation (GO:0043039)",
                0.00012548254326066126,
                158.37301587301587,
                1422.719267111456,
                ["AIMP1", "AIMP2"],
                0.00328427763325889,
                "GO_Biological_Process_2021",
            ],
            [
                4,
                "autophagy of mitochondrion (GO:0000422)",
                0.0001391618496477865,
                150.02255639097746,
                1332.181235322806,
                ["FUNDC2", "FUNDC1"],
                0.00328427763325889,
                "GO_Biological_Process_2021",
            ],
            [
                5,
                "tRNA aminoacylation for protein translation (GO:0006418)",
                0.00016099400163033776,
                139.02439024390245,
                1214.2589675263725,
                ["AIMP1", "AIMP2"],
                0.00328427763325889,
                "GO_Biological_Process_2021",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_diseases_drugs(self):
        df = enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="diseases_drugs",
        )
        result_to_test = df.values.tolist()[:5]
        expected_result = [
            [
                1,
                "Conduct disorder (maternal expressed emotions interaction)",
                0.003146171972763632,
                416.3541666666667,
                2398.8531813701998,
                ["MFHAS1"],
                0.029616604573025555,
                "GWAS_Catalog_2019",
            ],
            [
                2,
                "Skin aging (microtopography measurement)",
                0.00583590295074904,
                208.11458333333334,
                1070.4844511001115,
                ["MFHAS1"],
                0.029616604573025555,
                "GWAS_Catalog_2019",
            ],
            [
                3,
                "Immune reponse to smallpox (secreted IL-1beta)",
                0.006283564258093123,
                192.09615384615384,
                973.8925197668605,
                ["MTHFD2L"],
                0.029616604573025555,
                "GWAS_Catalog_2019",
            ],
            [
                4,
                "Mood instability",
                0.006283564258093123,
                192.09615384615384,
                973.8925197668605,
                ["MFHAS1"],
                0.029616604573025555,
                "GWAS_Catalog_2019",
            ],
            [
                5,
                "Personality dimensions",
                0.006731046493869444,
                178.36607142857142,
                892.013120036697,
                ["FUNDC1"],
                0.029616604573025555,
                "GWAS_Catalog_2019",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_celltypes(self):
        df = enrichr(["HAND1", "HAND2", "HEG1", "HRAT17"], database="celltypes")
        result_to_test = df.values.tolist()[:5]
        expected_result = [
            [
                1,
                "Cardiomyocytes",
                0.0004463780456458986,
                115.25581395348837,
                889.1230347572201,
                ["HAND2", "HAND1"],
                0.002231890228229493,
                "PanglaoDB_Augmented_2021",
            ],
            [
                2,
                "Airway Smooth Muscle Cells",
                0.019654783165825862,
                67.68027210884354,
                265.94519972543833,
                ["HAND2"],
                0.03218667151635187,
                "PanglaoDB_Augmented_2021",
            ],
            [
                3,
                "Vascular Smooth Muscle Cells",
                0.020442779355697382,
                65.01307189542484,
                252.90901184232303,
                ["HAND2"],
                0.03218667151635187,
                "PanglaoDB_Augmented_2021",
            ],
            [
                4,
                "Trophoblast Cells",
                0.025749337213081498,
                51.3359173126615,
                187.85590377250514,
                ["HAND1"],
                0.03218667151635187,
                "PanglaoDB_Augmented_2021",
            ],
            [
                5,
                "Enteric Neurons",
                0.03279117756136696,
                40.06262626262626,
                136.9178623130683,
                ["HAND2"],
                0.03279117756136696,
                "PanglaoDB_Augmented_2021",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_kinase_interactions(self):
        df = enrichr(
            ["HAND1", "HAND2", "HEG1", "HRAT17"], database="kinase_interactions"
        )
        result_to_test = df.values.tolist()[:5]
        expected_result = [
            [
                1,
                "PLK4",
                0.0005998942450539894,
                3332.3333333333335,
                24721.7718312461,
                ["HAND1"],
                0.0035993654703239363,
                "KEA_2015",
            ],
            [
                2,
                "PRKACG",
                0.034740452066662306,
                37.754285714285714,
                126.84875593802873,
                ["HAND1"],
                0.08269190935461691,
                "KEA_2015",
            ],
            [
                3,
                "AKT1",
                0.04134595467730846,
                31.558213716108455,
                100.53754798003185,
                ["HAND1"],
                0.08269190935461691,
                "KEA_2015",
            ],
            [
                4,
                "RPS6KA3",
                0.06476933555221323,
                19.803625377643506,
                54.20099786685237,
                ["HAND1"],
                0.08551838453962324,
                "KEA_2015",
            ],
            [
                5,
                "PRKACA",
                0.07631873969258186,
                16.670068027210885,
                42.88936390280122,
                ["HAND1"],
                0.08551838453962324,
                "KEA_2015",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_enrichr_plot(self):
        # Number of plots before running enrichr plot
        num_figures_before = plt.gcf().number
        enrichr(
            [
                "AIMP1",
                "MFHAS1",
                "BFAR",
                "FUNDC1",
                "AIMP2",
                "ASF1A",
                "FUNDC2",
                "TRMT112",
                "MTHFD2L",
            ],
            database="transcription",
            plot=True,
        )
        # Number of plots after running enrichr plot
        num_figures_after = plt.gcf().number

        self.assertGreater(
            num_figures_after,
            num_figures_before,
            "No matplotlib plt object was created.",
        )
