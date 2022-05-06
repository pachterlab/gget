import unittest
from gget.utils import (
    n_colors,
    aa_colors,
    get_uniprot_seqs,
    get_uniprot_info,
    rest_query,
    find_latest_ens_rel,
    gget_species_options,
    ref_species_options,
)

from gget.constants import UNIPROT_REST_API, ENSEMBL_REST_API

from .fixtures import LATEST_ENS_RELEASE, SPECIES_OPTIONS, REF_SPECIES_OPTIONS


class TestUtils(unittest.TestCase):
    def test_n_colors(self):
        result_to_test = n_colors("A")
        expected_result = "\x1b[38;5;15m\x1b[48;5;9mA\x1b[0;0m"

        self.assertEqual(result_to_test, expected_result)

    def test_aa_colors(self):
        result_to_test = aa_colors("L")
        expected_result = "\x1b[38;5;15m\x1b[48;5;12mL\x1b[0;0m"

        self.assertEqual(result_to_test, expected_result)

    def test_get_uniprot_seqs(self):
        df = get_uniprot_seqs(UNIPROT_REST_API, ["ENST00000392653", "ENST00000392657"])
        result_to_test = df.values
        expected_result = [
            [
                "P35326",
                "SPRR2A",
                "Homo sapiens (Human)",
                "MSYQQQQCKQPCQPPPVCPTPKCPEPCPPPKCPEPCPPPKCPQPCPPQQCQQKYPPVTPSPPCQSKYPPKSK",
                72,
                "ENST00000392653",
                list(["ENST00000392653"]),
            ],
            [
                "A7KAX9",
                "ARHGAP32 GRIT KIAA0712 RICS",
                "Homo sapiens (Human)",
                "METESESSTLGDDSVFWLESEVIIQVTDCEEEEREEKFRKMKSSVHSEEDDFVPELHRNVHPRERPDWEETLSAMARGADVPEIPGDLTLKTCGSTASMKVKHVKKLPFTKGHFPKMAECAHFHYENVEFGSIQLSLSEEQNEVMKNGCESKELVYLVQIACQGKSWIVKRSYEDFRVLDKHLHLCIYDRRFSQLSELPRSDTLKDSPESVTQMLMAYLSRLSAIAGNKINCGPALTWMEIDNKGNHLLVHEESSINTPAVGAAHVIKRYTARAPDELTLEVGDIVSVIDMPPKVLSTWWRGKHGFQVGLFPGHCVELINQKVPQSVTNSVPKPVSKKHGKLITFLRTFMKSRPTKQKLKQRGILKERVFGCDLGEHLLNSGFEVPQVLQSCTAFIERYGIVDGIYRLSGVASNIQRLRHEFDSEHVPDLTKEPYVQDIHSVGSLCKLYFRELPNPLLTYQLYEKFSDAVSAATDEERLIKIHDVIQQLPPPHYRTLEFLMRHLSLLADYCSITNMHAKNLAIVWAPNLLRSKQIESACFSGTAAFMEVRIQSVVVEFILNHVDVLFSGRISMAMQEGAASLSRPKSLLVSSPSTKLLTLEEAQARTQAQVNSPIVTENKYIEVGEGPAALQGKFHTIIEFPLERKRPQNKMKKSPVGSWRSFFNLGKSSSVSKRKLQRNESEPSEMKAMALKGGRAEGTLRSAKSEESLTSLHAVDGDSKLFRPRRPRSSSDALSASFNGEMLGNRCNSYDNLPHDNESEEEGGLLHIPALMSPHSAEDVDLSPPDIGVASLDFDPMSFQCSPPKAESECLESGASFLDSPGYSKDKPSANKKDAETGSSQCQTPGSTASSEPVSPLQEKLSPFFTLDLSPTEDKSSKPSSFTEKVVYAFSPKIGRKLSKSPSMSISEPISVTLPPRVSEVIGTVSNTTAQNASSSTWDKCVEERDATNRSPTQIVKMKTNETVAQEAYESEVQPLDQVAAEEVELPGKEDQSVSSSQSKAVASGQTQTGAVTHDPPQDSVPVSSVSLIPPPPPPKNVARMLALALAESAQQASTQSLKRPGTSQAGYTNYGDIAVATTEDNLSSSYSAVALDKAYFQTDRPAEQFHLQNNAPGNCDHPLPETTATGDPTHSNTTESGEQHHQVDLTGNQPHQAYLSGDPEKARITSVPLDSEKSDDHVSFPEDQSGKNSMPTVSFLDQDQSPPRFYSGDQPPSYLGASVDKLHHPLEFADKSPTPPNLPSDKIYPPSGSPEENTSTATMTYMTTTPATAQMSTKEASWDVAEQPTTADFAAATLQRTHRTNRPLPPPPSQRSAEQPPVVGQVQAATNIGLNNSHKVQGVVPVPERPPEPRAMDDPASAFISDSGAAAAQCPMATAVQPGLPEKVRDGARVPLLHLRAESVPAHPCGFPAPLPPTRMMESKMIAAIHSSSADATSSSNYHSFVTASSTSVDDALPLPLPVPQPKHASQKTVYSSFARPDVTTEPFGPDNCLHFNMTPNCQYRPQSVPPHHNKLEQHQVYGARSEPPASMGLRYNTYVAPGRNASGHHSKPCSRVEYVSSLSSSVRNTCYPEDIPPYPTIRRVQSLHAPPSSMIRSVPISRTEVPPDDEPAYCPRPLYQYKPYQSSQARSDYHVTQLQPYFENGRVHYRYSPYSSSSSSYYSPDGALCDVDAYGTVQLRPLHRLPNRDFAFYNPRLQGKSLYSYAGLAPRPRANVTGYFSPNDHNVVSMPPAADVKHTYTSWDLEDMEKYRMQSIRRESRARQKVKGPVMSQYDNMTPAVQDDLGGIYVIHLRSKSDPGKTGLLSVAEGKESRHAAKAISPEGEDRFYRRHPEAEMDRAHHHGGHGSTQPEKPSLPQKQSSLRSRKLPDMGCSLPEHRAHQEASHRQFCESKNGPPYPQGAGQLDYGSKGIPDTSEPVSYHNSGVKYAASGQESLRLNHKEVRLSKEMERPWVRQPSAPEKHSRDCYKEEEHLTQSIVPPPKPERSHSLKLHHTQNVERDPSVLYQYQPHGKRQSSVTVVSQYDNLEDYHSLPQHQRGVFGGGGMGTYVPPGFPHPQSRTYATALGQGAFLPAELSLQHPETQIHAE",
                2087,
                "ENST00000392657",
                list(["ENST00000392657"]),
            ],
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_get_uniprot_seqs_bad_type(self):
        result_to_test = get_uniprot_seqs(UNIPROT_REST_API, "banana")
        expected_result = None

        self.assertEqual(result_to_test, expected_result)

    def test_get_uniprot_info_gene(self):
        df = get_uniprot_info(
            UNIPROT_REST_API,
            ["ENSPFOG00000009362", "ENSGMOG00000033866"],
            id_type="Gene",
        )
        result_to_test = df.values
        expected_result = [
            ["A0A087XU86", nan, nan, "Fish-egg lectin-like", nan, "ENSPFOG00000009362"],
            [
                "A0A8C5B536",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
            [
                "A0A8C5B611",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
            [
                "A0A8C5BQR4",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
            [
                "A0A8C5C537",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
            [
                "A0A8C5CE04",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
            [
                "A0A8C5FPH9",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOG00000033866",
            ],
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_get_uniprot_info_transcript(self):
        df = get_uniprot_info(
            UNIPROT_REST_API,
            ["ENSPFOT00000009352", "ENSGMOT00000059188"],
            id_type="Transcript",
        ).values
        result_to_test = df.values
        expected_result = [
            ["A0A087XU86", nan, nan, "Fish-egg lectin-like", nan, "ENSPFOT00000009352"],
            [
                "A0A8C5B536",
                "LOC115543029",
                list(["LOC115543029"]),
                "Fish-egg lectin-like",
                nan,
                "ENSGMOT00000059188",
            ],
        ]

        self.assertEqual(result_to_test, expected_result)

    def test_get_uniprot_info_bad_type(self):
        result_to_test = get_uniprot_info(UNIPROT_REST_API, "banana")
        expected_result = None

        self.assertEqual(result_to_test, expected_result)

    def test_rest_query(self):
        server = ENSEMBL_REST_API
        content_type = "application/json"
        ensembl_ID = "ENSPFOG00000009362"
        result_to_test = rest_query(server, query, content_type)
        expected_result = {
            "source": "ensembl",
            "id": "ENSPFOG00000009362",
            "db_type": "core",
            "species": "poecilia_formosa",
            "canonical_transcript": "ENSPFOT00000009352.1",
            "end": 504269,
            "logic_name": "ensembl",
            "version": 1,
            "strand": 1,
            "biotype": "protein_coding",
            "assembly_name": "PoeFor_5.1.2",
            "object_type": "Gene",
            "description": "fish-egg lectin-like [Source:NCBI gene;Acc:103152571]",
            "start": 502511,
            "seq_region_name": "KI519928.1",
        }

        self.assertEqual(result_to_test, expected_result)

    def test_rest_query_bad_type(self):
        server = ENSEMBL_REST_API
        content_type = "application/json"
        ensembl_ID = "banana"
        with self.assertRaises(RuntimeError):
            result = rest_query(server, query, content_type)

    def test_find_latest_ens_rel(self):
        result_to_test = find_latest_ens_rel()
        expected_result = LATEST_ENS_RELEASE

        self.assertEqual(result_to_test, expected_result)

    def test_gget_species_options(self):
        result_to_test = gget_species_options(release=106)
        expected_result = SPECIES_OPTIONS

        self.assertEqual(result_to_test, expected_result)

    def test_gget_species_options_bad_type(self):
        with self.assertRaises(ValueError):
            result = gget_species_options(release=2000)

    def test_ref_species_options(self):
        result_to_test = ref_species_options("gtf", 105)
        expected_result = REF_SPECIES_OPTIONS

        self.assertEqual(result_to_test, expected_result)

    def test_ref_species_options_bad_type(self):
        with self.assertRaises(ValueError):
            result = ref_species_options("gtf", release=2000)
