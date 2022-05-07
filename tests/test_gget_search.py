import unittest
from gget.gget_search import search


class TestSearch(unittest.TestCase):
    def test_search_gene_one_sw(self):
        searchwords = "swiss"
        species = "drosophila"
        df = search(
            searchwords,
            species,
            seqtype="gene",
            limit=None,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "FBgn0003656",
                "sws",
                "swiss cheese",
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBgn0003656",
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_or(self):
        searchwords = ["swiss", "cheese"]
        species = "drosophi"
        df = search(
            searchwords,
            species,
            seqtype="gene",
            andor="or",
            limit=None,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "FBgn0003656",
                "sws",
                "swiss cheese",
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBgn0003656",
            ],
            [
                "FBgn0043362",
                "bchs",
                "blue cheese",
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBgn0043362",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_and(self):
        searchwords = ["swiss", "cheese"]
        species = "droso"
        df = search(
            searchwords,
            species,
            seqtype="gene",
            andor="and",
            limit=None,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "FBgn0003656",
                "sws",
                "swiss cheese",
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBgn0003656",
            ]
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_one_sw_limit(self):
        searchwords = "fun"
        species = "human"
        limit = 3
        df = search(
            searchwords,
            species,
            seqtype="gene",
            limit=limit,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "ENSG00000069509",
                "FUNDC1",
                "FUN14 domain containing 1 [Source:HGNC Symbol;Acc:HGNC:28746]",
                "FUN14 domain containing 1",
                "protein_coding",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000069509",
            ],
            [
                "ENSG00000084754",
                "HADHA",
                "hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit alpha [Source:HGNC Symbol;Acc:HGNC:4801]",
                "hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit alpha",
                "protein_coding",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000084754",
            ],
            [
                "ENSG00000103429",
                "BFAR",
                "bifunctional apoptosis regulator [Source:HGNC Symbol;Acc:HGNC:17613]",
                "bifunctional apoptosis regulator",
                "protein_coding",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000103429",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_one_sw(self):
        searchwords = "nep3"
        species = "drosophila_melanogaster"
        df = search(
            searchwords,
            species,
            seqtype="transcript",
            limit=None,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "FBtr0070000",
                "Nep3-RA",
                None,
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBtr0070000",
            ],
            [
                "FBtr0307554",
                "Nep3-RB",
                None,
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBtr0307554",
            ],
            [
                "FBtr0307555",
                "Nep3-RC",
                None,
                None,
                "protein_coding",
                "https://uswest.ensembl.org/drosophila_melanogaster/Gene/Summary?g=FBtr0307555",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_two_sw_or(self):
        searchwords = ["SNORA71", "201"]
        species = "accipiter"
        df = search(
            searchwords,
            species,
            seqtype="transcript",
            andor="or",
            limit=5,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "ENSANIT00000001698",
                "SNORA71-201",
                None,
                "Small nucleolar RNA SNORA71",
                "snoRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001698",
            ],
            [
                "ENSANIT00000001705",
                "SNORA71-201",
                None,
                "Small nucleolar RNA SNORA71",
                "snoRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001705",
            ],
            [
                "ENSANIT00000000884",
                "RNaseP_nuc-201",
                None,
                "Nuclear RNase P",
                "ribozyme",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000000884",
            ],
            [
                "ENSANIT00000001802",
                "5_8S_rRNA-201",
                None,
                "5.8S ribosomal RNA",
                "rRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001802",
            ],
            [
                "ENSANIT00000001839",
                "U2-201",
                None,
                "U2 spliceosomal RNA",
                "snRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001839",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_transcript_two_sw_and(self):
        searchwords = ["SNORA71", "201"]
        species = "accipiter"
        df = search(
            searchwords,
            species,
            seqtype="transcript",
            andor="and",
            limit=None,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "ENSANIT00000001698",
                "SNORA71-201",
                None,
                "Small nucleolar RNA SNORA71",
                "snoRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001698",
            ],
            [
                "ENSANIT00000001705",
                "SNORA71-201",
                None,
                "Small nucleolar RNA SNORA71",
                "snoRNA",
                "https://uswest.ensembl.org/accipiter_nisus/Gene/Summary?g=ENSANIT00000001705",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_search_gene_two_sw_limit(self):
        searchwords = ["SNORA71", "201"]
        species = "homo"
        limit = 3
        df = search(
            searchwords,
            species,
            seqtype="transcript",
            andor="and",
            limit=limit,
        )

        result_to_test = df.values.tolist()
        expected_result = [
            [
                "ENST00000362582",
                "SNORA71.3-201",
                None,
                None,
                "snoRNA",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENST00000362582",
            ],
            [
                "ENST00000363484",
                "SNORA71D-201",
                None,
                "small nucleolar RNA, H/ACA box 71D",
                "snoRNA",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENST00000363484",
            ],
            [
                "ENST00000364523",
                "SNORA71.1-201",
                None,
                None,
                "snoRNA",
                "https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENST00000364523",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)


    ## Test bad input errors
    def test_search_gene_bad_species(self):
        searchwords = "fun"
        species = "banana"

        with self.assertRaises(ValueError):
            search(
                searchwords,
                species,
                seqtype="gene",
                andor="and",
                limit=None,
            )

    def test_search_transcript_bad_species(self):
        searchwords = "fun"
        species = "banana"

        with self.assertRaises(ValueError):
            search(
                searchwords,
                species,
                seqtype="transcript",
                andor="or",
                limit=None,
            )

    def test_search_gene_bad_andor(self):
        searchwords = "fun"
        species = "mouse"
        andor = "sneeze"

        with self.assertRaises(ValueError):
            search(
                searchwords,
                species,
                seqtype="gene",
                andor=andor,
                limit=None,
            )

    def test_search_transcript_bad_andor(self):
        searchwords = "fun"
        species = "mouse"
        andor = "sneeze"

        with self.assertRaises(ValueError):
            search(
                searchwords,
                species,
                seqtype="transcript",
                andor=andor,
                limit=None,
            )
