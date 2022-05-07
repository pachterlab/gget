import unittest
import numpy as np
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
        result_to_test = df.values.tolist()
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

        self.assertListEqual(result_to_test, expected_result)

    def test_get_uniprot_seqs_bad_type(self):
        result_to_test = get_uniprot_seqs(UNIPROT_REST_API, "banana")
        # Expect an empty data frame
        self.assertTrue(result_to_test.empty)

    def test_get_uniprot_info_gene(self):
        df = get_uniprot_info(
            UNIPROT_REST_API,
            ["ENSG00000187140", "ENSG00000181449"],
            id_type="Gene",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            [
                "Q9UJU5",
                "FOXD3",
                ["FOXD3", "HFH2"],
                "Forkhead box protein D3 (HNF3/FH transcription factor genesis)",
                "FUNCTION: Binds to the consensus sequence 5'-A[AT]T[AG]TTTGTTT-3' and acts as a transcriptional repressor. Also acts as a transcriptional activator. Promotes development of neural crest cells from neural tube progenitors. Restricts neural progenitor cells to the neural crest lineage while suppressing interneuron differentiation. Required for maintenance of pluripotent cells in the pre-implantation and peri-implantation stages of embryogenesis. {ECO:0000269|PubMed:11891324}.",
                "ENSG00000187140",
            ],
            [
                "P48431",
                "SOX2",
                ["SOX2"],
                "Transcription factor SOX-2",
                "FUNCTION: Transcription factor that forms a trimeric complex with OCT4 on DNA and controls the expression of a number of genes involved in embryonic development such as YES1, FGF4, UTF1 and ZFP206 (By similarity). Binds to the proximal enhancer region of NANOG (By similarity). Critical for early embryogenesis and for embryonic stem cell pluripotency (PubMed:18035408). Downstream SRRT target that mediates the promotion of neural stem cell self-renewal (By similarity). Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation (By similarity). May function as a switch in neuronal development (By similarity). {ECO:0000250|UniProtKB:P48430, ECO:0000250|UniProtKB:P48432, ECO:0000269|PubMed:18035408}.",
                "ENSG00000181449",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_get_uniprot_info_transcript(self):
        df = get_uniprot_info(
            UNIPROT_REST_API,
            ["ENST00000325404", "ENST00000251020"],
            id_type="Transcript",
        )
        result_to_test = df.values.tolist()
        expected_result = [
            [
                "P48431",
                "SOX2",
                ["SOX2"],
                "Transcription factor SOX-2",
                "FUNCTION: Transcription factor that forms a trimeric complex with OCT4 on DNA and controls the expression of a number of genes involved in embryonic development such as YES1, FGF4, UTF1 and ZFP206 (By similarity). Binds to the proximal enhancer region of NANOG (By similarity). Critical for early embryogenesis and for embryonic stem cell pluripotency (PubMed:18035408). Downstream SRRT target that mediates the promotion of neural stem cell self-renewal (By similarity). Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation (By similarity). May function as a switch in neuronal development (By similarity). {ECO:0000250|UniProtKB:P48430, ECO:0000250|UniProtKB:P48432, ECO:0000269|PubMed:18035408}.",
                "ENST00000325404",
            ],
            [
                "Q9NSC2",
                "SALL1",
                ["SALL1", "SAL1", "ZNF794"],
                "Sal-like protein 1 (Spalt-like transcription factor 1) (Zinc finger protein 794) (Zinc finger protein SALL1) (Zinc finger protein Spalt-1) (HSal1) (Sal-1)",
                "FUNCTION: Transcriptional repressor involved in organogenesis. Plays an essential role in ureteric bud invasion during kidney development. {ECO:0000250|UniProtKB:Q9ER74}.",
                "ENST00000251020",
            ],
        ]

        self.assertListEqual(result_to_test, expected_result)

    def test_get_uniprot_info_bad_type(self):
        result_to_test = get_uniprot_info(
            UNIPROT_REST_API, "banana", id_type="Transcript"
        )
        # expected_result = None
        self.assertIsNone(result_to_test)

    def test_rest_query(self):
        server = ENSEMBL_REST_API
        content_type = "application/json"
        ensembl_ID = "ENSPFOG00000009362"
        query = "lookup/id/" + ensembl_ID + "?"
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
        query = "lookup/id/" + ensembl_ID + "?"
        with self.assertRaises(RuntimeError):
            rest_query(server, query, content_type)

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
