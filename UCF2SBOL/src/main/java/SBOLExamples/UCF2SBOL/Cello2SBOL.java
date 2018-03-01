package SBOLExamples.UCF2SBOL;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TimeZone;

import javax.xml.namespace.QName;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.Activity;
import org.sbolstandard.core2.Component;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.DirectionType;
import org.sbolstandard.core2.FunctionalComponent;
import org.sbolstandard.core2.Interaction;
import org.sbolstandard.core2.ModuleDefinition;
import org.sbolstandard.core2.OrientationType;
import org.sbolstandard.core2.SBOLConversionException;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.SBOLValidate;
import org.sbolstandard.core2.SBOLValidationException;
import org.sbolstandard.core2.Sequence;
import org.sbolstandard.core2.SequenceAnnotation;
import org.sbolstandard.core2.SequenceOntology;
import org.sbolstandard.core2.SystemsBiologyOntology;
import org.synbiohub.frontend.SynBioHubException;
import org.synbiohub.frontend.SynBioHubFrontend;

public class Cello2SBOL {
	
	static String uriPrefix = "http://cellocad.org/"; 
	static String version = "1";
	static URI derivedFrom = URI.create("https://github.com/CIDARLAB/cello/blob/master/resources/UCF/Eco1C1G1T0.UCF.json");
	static String so = "http://identifiers.org/so/";
	static String provNS = "http://www.w3.org/ns/prov#";
	static String dcNS = "http://purl.org/dc/elements/1.1/";
	static String dcTermsNS = "http://purl.org/dc/terms/";
	static String celloNS = "http://cellocad.org/Terms/cello#";

	static URI activityURI;
	static String createdDate;
	
	
	
    static final String obj_ver = "1.0";
    
    private static void createSensorsReporters(SBOLDocument doc) throws URISyntaxException,
    	SBOLValidationException, SBOLConversionException, IOException 
	{ 
    	// TODO: add wasDerivedFrom/wasGeneratedBy
		//creating ModuleDefinitions 
		ModuleDefinition LacI_protein_production = doc.createModuleDefinition("LacI_protein_production", obj_ver);
		ModuleDefinition LacI_pTac_repression = doc.createModuleDefinition("LacI_pTac_repression", obj_ver);
		ModuleDefinition LacI_protein_degradation = doc.createModuleDefinition("LacI_protein_degradation", obj_ver);

		ModuleDefinition TetR_protein_production = doc.createModuleDefinition("TetR_protein_production", obj_ver);
		ModuleDefinition TetR_pTet_repression = doc.createModuleDefinition("TetR_pTet_repression", obj_ver);
		ModuleDefinition TetR_protein_degradation = doc.createModuleDefinition("TetR_protein_degradation", obj_ver);
		
		ModuleDefinition AraC_protein_production = doc.createModuleDefinition("AraC_protein_production", obj_ver);
		ModuleDefinition AraC_pBAD_repression = doc.createModuleDefinition("AraC_pBAD_repression", obj_ver);

		ModuleDefinition LuxR_protein_production = doc.createModuleDefinition("LuxR_protein_production", obj_ver);
		ModuleDefinition LuxI_protein_production = doc.createModuleDefinition("LuxI_protein_production", obj_ver);

		ModuleDefinition LuxRLuxI_complex = doc.createModuleDefinition("LuxRLuxI_complexFormation", obj_ver);
		ModuleDefinition LuxRLuxI_pLuxStar_repression = doc.createModuleDefinition("LuxRLuxI_pLuxStar_repression", obj_ver);
	
		ModuleDefinition IPTG_LacI_binding = doc.createModuleDefinition("IPTG_LacI_complexFormation", obj_ver);
		ModuleDefinition aTc_TetR_binding = doc.createModuleDefinition("aTc_TetR_complexFormation", obj_ver);
		ModuleDefinition IPTG_LacI_complex_degradation = doc.createModuleDefinition("IPTG_LacI_complex_degradation", obj_ver);
		ModuleDefinition aTc_TetR_complex_degradation = doc.createModuleDefinition("aTc_TetR_complex_degradation", obj_ver);

		ModuleDefinition GFP_protein_production = doc.createModuleDefinition("GFP_protein_production", obj_ver);
		ModuleDefinition RFP_protein_production = doc.createModuleDefinition("RFP_protein_production", obj_ver);
		ModuleDefinition BFP_protein_production = doc.createModuleDefinition("BFP_protein_production", obj_ver);
		ModuleDefinition YFP_protein_production = doc.createModuleDefinition("YFP_protein_production", obj_ver);

		ModuleDefinition GFP_protein_degradation = doc.createModuleDefinition("GFP_protein_degradation", obj_ver);
		ModuleDefinition RFP_protein_degradation = doc.createModuleDefinition("RFP_protein_degradation", obj_ver);
		ModuleDefinition BFP_protein_degradation = doc.createModuleDefinition("BFP_protein_degradation", obj_ver);
		ModuleDefinition YFP_protein_degradation = doc.createModuleDefinition("YFP_protein_degradation", obj_ver);
		
		ModuleDefinition BetI_protein_degradation = doc.createModuleDefinition("BetI_protein_degradation", obj_ver);
		ModuleDefinition AmtR_protein_degradation = doc.createModuleDefinition("AmtR_protein_degradation", obj_ver);
		
		
		//Creating output ComponentDefinitions 
		ComponentDefinition YFP_protein = createComponentDefinition(doc, "YFP_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition RFP_protein = createComponentDefinition(doc, "RFP_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition BFP_protein = createComponentDefinition(doc, "BFP_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition GFP_protein = createComponentDefinition(doc, "GFP_protein", ComponentDefinition.PROTEIN, null);

		//Creating protein ComponentDefinitions
		ComponentDefinition LacI_protein = createComponentDefinition(doc, "LacI_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition TetR_protein = createComponentDefinition(doc, "TetR_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition AraC_protein = createComponentDefinition(doc, "AraC_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition LuxI_protein = createComponentDefinition(doc, "LuxI_protein", ComponentDefinition.PROTEIN, null);
		ComponentDefinition LuxR_protein = createComponentDefinition(doc, "LuxR_protein", ComponentDefinition.PROTEIN, null);
		
		//Creating small molecules ComponentDefinitions
		ComponentDefinition aTc_smallMolecule = createComponentDefinition(doc, "aTc_small_molecule", ComponentDefinition.SMALL_MOLECULE, null);
		ComponentDefinition IPTG_smallMolecule = createComponentDefinition(doc, "IPTG_small_molecule", ComponentDefinition.SMALL_MOLECULE, null);
		ComponentDefinition aTc_TetR_complex = createComponentDefinition(doc, "aTc_TetR_complex", ComponentDefinition.COMPLEX, null);
		ComponentDefinition IPTG_LacI_complex = createComponentDefinition(doc, "IPTG_LacI_complex", ComponentDefinition.COMPLEX, null);
		
		//Creating CDS ComponentDefinitions
		ComponentDefinition lacI_cds = createCDS(doc, "lacI");
		ComponentDefinition tetR_cds = createCDS(doc, "tetR");
		ComponentDefinition araC_cds = createCDS(doc, "araC");
		ComponentDefinition luxI_cds = createCDS(doc, "luxI");
		ComponentDefinition luxR_cds = createCDS(doc, "luxR");
		
		ComponentDefinition gfp_cds = createCDS(doc, "gfp");
		ComponentDefinition rfp_cds = createCDS(doc, "rfp");
		ComponentDefinition bfp_cds = createCDS(doc, "bfp");
		ComponentDefinition yfp_cds = createCDS(doc, "yfp");

		//Creating promoters ComponentDefinitions
		ComponentDefinition pTac_prom = createPromoter(doc, "pTac");
		ComponentDefinition pTet_prom = createPromoter(doc, "pTet");
		ComponentDefinition pBAD_prom = createPromoter(doc, "pBAD");
		ComponentDefinition pLuxStar_prom = createPromoter(doc, "pLuxStar");
		ComponentDefinition pCONST_prom = createPromoter(doc, "pCONST_prom");

		ComponentDefinition LuxRLuxI_comp = doc.createComponentDefinition("LuxRLuxI_complex", obj_ver, ComponentDefinition.COMPLEX);

		//creating FunctionalComponents to store in appropriate ModuleDefinitions 
		FunctionalComponent lacI_cds_func = createFunctionalComponent(LacI_protein_production, "lacI_functionalComponent", lacI_cds);
		FunctionalComponent LacI_prot_func = createFunctionalComponent(LacI_protein_production, "LacI_protein_functionalComponent", LacI_protein);
		FunctionalComponent LacI_prot2_func = createFunctionalComponent(LacI_pTac_repression,"LacI_protein_functionalComponent", LacI_protein);
		FunctionalComponent LacI_prot4_func = createFunctionalComponent(LacI_protein_degradation,"LacI_protein_functionalComponent", LacI_protein);
		
		FunctionalComponent pTac_func = createFunctionalComponent(LacI_pTac_repression, "pTac_functionalComponent", pTac_prom);

		FunctionalComponent tetR_cds_func = createFunctionalComponent(TetR_protein_production, "tetR_functionalComponent", tetR_cds);
		FunctionalComponent TetR_prot_func = createFunctionalComponent(TetR_protein_production, "TetR_protein_functionalComponent", TetR_protein);
		FunctionalComponent TetR_prot2_func = createFunctionalComponent(TetR_pTet_repression, "TetR_protein_functionalComponent", TetR_protein);
		FunctionalComponent TetR_prot4_func = createFunctionalComponent(TetR_protein_degradation, "TetR_protein_functionalComponent", TetR_protein);

		FunctionalComponent pTet_func = createFunctionalComponent(TetR_pTet_repression, "pTet_functionalComponent", pTet_prom);

		FunctionalComponent araC_cds_func = createFunctionalComponent(AraC_protein_production, "araC_functionalComponent", araC_cds);
		FunctionalComponent AraC_prot_func = createFunctionalComponent(AraC_protein_production, "AraC_protein_functionalComponent", AraC_protein);

		FunctionalComponent AraC_prot2_func = createFunctionalComponent(AraC_pBAD_repression, "AraC_protein_functionalComponent", AraC_protein);
		FunctionalComponent pBAD_func = createFunctionalComponent(AraC_pBAD_repression, "pBAD_functionalComponent", pBAD_prom);

		FunctionalComponent luxR_cds_func = createFunctionalComponent(LuxR_protein_production, "luxR_functionalComponent", luxR_cds);
		FunctionalComponent LuxR_prot_func = createFunctionalComponent(LuxR_protein_production, "LuxR_protein_functionalComponent", LuxR_protein);

		FunctionalComponent luxI_cds_func = createFunctionalComponent(LuxI_protein_production, "luxI_functionalComponent", luxI_cds);
		FunctionalComponent LuxI_prot_func = createFunctionalComponent(LuxI_protein_production, "LuxI_protein_functionalComponent", LuxI_protein);

		FunctionalComponent LuxILuxR_func = createFunctionalComponent(LuxRLuxI_pLuxStar_repression, "LuxRLuxI_complex_functionalComponent", LuxRLuxI_comp);
		FunctionalComponent pluxStar_func = createFunctionalComponent(LuxRLuxI_pLuxStar_repression, "pluxStar_functionalComponent", pLuxStar_prom);

		FunctionalComponent LuxI_prot2_func = createFunctionalComponent(LuxRLuxI_complex, "LuxI_protein_functionalComponent", LuxI_protein);
		FunctionalComponent LuxR_prot2_func = createFunctionalComponent(LuxRLuxI_complex, "LuxR_protein_functionalComponent", LuxR_protein);
		FunctionalComponent LuxILuxR2_func = createFunctionalComponent(LuxRLuxI_complex, "LuxRLuxI_complex_functionalComponent", LuxRLuxI_comp);

		FunctionalComponent IPTG_smallMolecule_func = createFunctionalComponent(IPTG_LacI_binding, "IPTG_smallMolecule_functionalComponent", IPTG_smallMolecule);
		FunctionalComponent LacI_prot3_func = createFunctionalComponent(IPTG_LacI_binding, "LacI_protein_functionalComponent", LacI_protein);
		FunctionalComponent IPTG_LacI_comp_func = createFunctionalComponent(IPTG_LacI_binding, "IPTG_LacI_complex_functionalComponent", IPTG_LacI_complex);
		FunctionalComponent IPTG_LacI_complex_deg_func = createFunctionalComponent(IPTG_LacI_complex_degradation, "IPTG_LacI_complex_deg_functionalComponent", IPTG_LacI_complex);
		
		FunctionalComponent aTc_smallMolecule_func = createFunctionalComponent(aTc_TetR_binding, "aTc_smallMolecule_functionalComponent", aTc_smallMolecule);
		FunctionalComponent TetR_prot3_func = createFunctionalComponent(aTc_TetR_binding, "TetR_protein_functionalComponent", TetR_protein);
		FunctionalComponent aTc_TetR_comp_func = createFunctionalComponent(aTc_TetR_binding, "aTc_TetR_complex_functionalComponent", aTc_TetR_complex);
		FunctionalComponent aTc_TetR_complex_deg_func = createFunctionalComponent(aTc_TetR_complex_degradation, "aTc_TetR_complex_deg_functionalComponent", aTc_TetR_complex);
		
		FunctionalComponent gfp_cds_func = createFunctionalComponent(GFP_protein_production, "gfp_functionalComponent", gfp_cds);
		FunctionalComponent GFP_prot_func = createFunctionalComponent(GFP_protein_production, "GFP_protein_functionalComponent", GFP_protein);
		FunctionalComponent GFP_prot2_func = createFunctionalComponent(GFP_protein_degradation, "GFP_protein_functionalComponent", GFP_protein);
		

		FunctionalComponent rfp_cds_func = createFunctionalComponent(RFP_protein_production, "rfp_functionalComponent", rfp_cds);
		FunctionalComponent RFP_prot_func = createFunctionalComponent(RFP_protein_production, "RFP_protein_functionalComponent", RFP_protein);
		FunctionalComponent RFP_prot2_func = createFunctionalComponent(RFP_protein_degradation, "RFP_protein_functionalComponent", RFP_protein);
		
		FunctionalComponent bfp_cds_func = createFunctionalComponent(BFP_protein_production, "bfp_functionalComponent", bfp_cds);
		FunctionalComponent BFP_prot_func = createFunctionalComponent(BFP_protein_production, "BFP_protein_functionalComponent", BFP_protein);
		FunctionalComponent BFP_prot2_func = createFunctionalComponent(BFP_protein_degradation, "BFP_protein_functionalComponent", BFP_protein);
		
		FunctionalComponent yfp_cds_func = createFunctionalComponent(YFP_protein_production, "yfp_functionalComponent", yfp_cds);
		FunctionalComponent YFP_prot_func = createFunctionalComponent(YFP_protein_production, "YFP_protein_functionalComponent", YFP_protein);
		FunctionalComponent YFP_prot2_func = createFunctionalComponent(YFP_protein_degradation, "YFP_protein_functionalComponent", YFP_protein);

		ComponentDefinition AmtR_prot = doc.getComponentDefinition("AmtR_protein", "1"); 
		ComponentDefinition BetI_prot = doc.getComponentDefinition("BetI_protein", "1"); 
		FunctionalComponent AmtR_prot_func = createFunctionalComponent(AmtR_protein_degradation, "AmtR_protein_functionalComponent", AmtR_prot);
		FunctionalComponent BetI_prot_func = createFunctionalComponent(BetI_protein_degradation, "BetI_protein_functionalComponent", BetI_prot);

		
		//creating LacI Interactions
		createProductionInteraction(LacI_protein_production, "LacI", lacI_cds_func, LacI_prot_func);
		createInhibitionInteraction(LacI_pTac_repression, "LacI_pTac", LacI_prot2_func, pTac_func);
		createDegradationInteraction(LacI_protein_degradation, LacI_prot4_func);

		//creating TetR Interactions
		createProductionInteraction(TetR_protein_production, "TetR", tetR_cds_func, TetR_prot_func);
		createInhibitionInteraction(TetR_pTet_repression, "TetR_pTet", TetR_prot2_func, pTet_func);
		createDegradationInteraction(TetR_protein_degradation, TetR_prot4_func);

		//creating AraC Interactions
		createProductionInteraction(AraC_protein_production, "AraC", araC_cds_func, AraC_prot_func);
		createInhibitionInteraction(AraC_pBAD_repression, "AraC_pBad", AraC_prot2_func, pBAD_func);
		
		//creating LuxR Interactions
		createProductionInteraction(LuxR_protein_production, "LuxR", luxR_cds_func, LuxR_prot_func); 
		
		//creating LuxI Interactions
		createProductionInteraction(LuxI_protein_production, "LuxI", luxI_cds_func, LuxI_prot_func); 
		
		//creating luxR-luxI complex Interactions
		createComplexInteraction(LuxRLuxI_complex, LuxI_prot2_func, LuxR_prot2_func, LuxILuxR2_func);
		createInhibitionInteraction(LuxRLuxI_pLuxStar_repression, "LuxRLuxI_pLuxStar", LuxILuxR_func, pluxStar_func);

		//creating IPTG Interactions
		createComplexInteraction(IPTG_LacI_binding, IPTG_smallMolecule_func, LacI_prot3_func, IPTG_LacI_comp_func);
		createDegradationInteraction(IPTG_LacI_complex_degradation, IPTG_LacI_complex_deg_func);
		
		//creating aTc-TetR Interactions
		createComplexInteraction(aTc_TetR_binding, aTc_smallMolecule_func, TetR_prot3_func, aTc_TetR_comp_func);
		createDegradationInteraction(aTc_TetR_complex_degradation, aTc_TetR_complex_deg_func);
		
		//creating output Interactions
		createProductionInteraction(GFP_protein_production, "GFP", gfp_cds_func, GFP_prot_func);
		createProductionInteraction(RFP_protein_production, "RFP", rfp_cds_func, RFP_prot_func);
		createProductionInteraction(YFP_protein_production, "YFP", yfp_cds_func, YFP_prot_func);
		createProductionInteraction(BFP_protein_production, "BFP", bfp_cds_func, BFP_prot_func);
		
		createDegradationInteraction(GFP_protein_degradation, GFP_prot2_func);
		createDegradationInteraction(RFP_protein_degradation, RFP_prot2_func);
		createDegradationInteraction(YFP_protein_degradation, YFP_prot2_func);
		createDegradationInteraction(BFP_protein_degradation, BFP_prot2_func);
		
		//creating BetI and AmtR Interactions
		createDegradationInteraction(AmtR_protein_degradation, AmtR_prot_func);
		createDegradationInteraction(BetI_protein_degradation, BetI_prot_func);
		
		//creating promoter Sequences 
		
		String pCONST_prom_sequences = "GATAAGTCCCTAACTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC";
		Sequence pCONST_prom_seq= doc.createSequence("pCONST_sequence", obj_ver, pCONST_prom_sequences, Sequence.IUPAC_DNA);
		pCONST_prom.addSequence(pCONST_prom_seq);
		
		String pTac_prom_sequences = "AACGATCGTTGGCTGTGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGCTCACAATT";
		Sequence pTac_prom_seq= doc.createSequence("pTac_sequence", obj_ver, pTac_prom_sequences, Sequence.IUPAC_DNA);
		pTac_prom.addSequence(pTac_prom_seq);

		String pTet_prom_sequences = "TACTCCACCGTTGGCTTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCAC";
		Sequence pTet_prom_seq= doc.createSequence("pTet_sequence", obj_ver, pTet_prom_sequences, Sequence.IUPAC_DNA);
		pTet_prom.addSequence(pTet_prom_seq);

		String pBAD_prom_sequences = "ACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGC";
		Sequence pBAD_prom_seq= doc.createSequence("pBAD_sequence", obj_ver, pBAD_prom_sequences, Sequence.IUPAC_DNA);
		pBAD_prom.addSequence(pBAD_prom_seq);

		String pLuxStar_prom_sequences = "ATAGCTTCTTACCGGACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTACTTTCGAATAAA";
		Sequence pLuxStar_prom_seq= doc.createSequence("pLuxStar_sequence", obj_ver, pLuxStar_prom_sequences, Sequence.IUPAC_DNA);
		pLuxStar_prom.addSequence(pLuxStar_prom_seq);
		
		//creating input sensor sequences
		String araC_cds_sequences = "atggctgaagcgcaaaatgatcccctgctgccgggatactcgtttaatgcccatctggtggcgggtttaacgccgattgaggccaacggttatctcgatttttttatcgaccgaccgctgggaatgaaaggttatattctcaatctcaccattcgcggtcagggggtggtgaaaaatcagggacgagaatttgtttgccgaccgggtgatattttgctgttcccgccaggagagattcatcactacggtcgtcatccggaggctcgcgaatggtatcaccagtgggtttactttcgtccgcgcgcctactggcatgaatggcttaactggccgtcaatatttgccaatacggggttctttcgcccggatgaagcgcaccagccgcatttcagcgacctgtttgggcaaatcattaacgccgggcaaggggaagggcgctattcggagctgctggcgataaatctgcttgagcaattgttactgcggcgcatggaagcgattaacgagtcgctccatccaccgatggataatcgggtacgcgaggcttgtcagtacatcagcgatcacctggcagacagcaattttgatatcgccagcgtcgcacagcatgtttgcttgtcgccgtcgcgtctgtcacatcttttccgccagcagttagggattagcgtcttaagctggcgcgaggaccaacgtatcagccaggcgaagctgcttttgagcaccacccggatgcctatcgccaccgtcggtcgcaatgttggttttgacgatcaactctatttctcgcgggtatttaaaaaatgcaccggggccagcccgagcgagttccgtgccggttgtgaagaaaaagtgaatgatgtagccgtcaagttgtcataa"; 
		Sequence araC_cds_seq = doc.createSequence("araC_sequence", obj_ver, araC_cds_sequences, Sequence.IUPAC_DNA);
		araC_cds.addSequence(araC_cds_seq);

		String tetR_cds_sequences = "atgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaatgaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgccttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctttactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggtttttcactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacctactactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggattagaaaaacaacttaaatgtgaaagtgggtcctaa";
		Sequence tetR_cds_seq = doc.createSequence("tetR_sequence", obj_ver, tetR_cds_sequences, Sequence.IUPAC_DNA);
		tetR_cds.addSequence(tetR_cds_seq);

		String lacI_cds_sequences = "atgaaaccagtaacgttatacgatgtcgcagagtatgccggtgtctcttatcagaccgtttcccgcgt"+
				"ggtgaaccaggccagccacgtttctgcgaaaacgcgggaaaaagtggaagcggcgatggcggagctgaattacattcccaaccgcgtg"+
				"gcacaacaactggcgggcaaacagtcgttgctgattggcgttgccacctccagtctggccctgcacgcgccgtcgcaaattgtcgcgg"+
				"cgattaaatctcgcgccgatcaactgggtgccagcgtggtggtgtcgatggtagaacgaagcggcgtcgaagcctgtaaagcggcggt"+
				"gcacaatcttctcgcgcaacgcgtcagtgggctgatcattaactatccgctggatgaccaggatgccattgctgtggaagctgcctgc"+
				"actaatgttccggcgttatttcttgatgtctctgaccagacacccatcaacagtattattttctcccatgaggacggtacgcgactgg"+
				"gcgtggagcatctggtcgcattgggtcaccagcaaatcgcgctgttagcgggcccattaagttctgtctcggcgcgtctgcgtctggc"+
				"tggctggcataaatatctcactcgcaatcaaattcagccgatagcggaacgggaaggcgactggagtgccatgtccggttttcaacaa"+
				"accatgcaaatgctgaatgagggcatcgttcccactgcgatgctggttgccaacgatcagatggcgctgggcgcaatgcgcgccatta"+
				"ccgagtccgggctgcgcgttggtgcggatatctcggtagtgggatacgacgataccgaagatagctcatgttatatcccgccgttaac"+
				"caccatcaaacaggattttcgcctgctggggcaaaccagcgtggaccgcttgctgcaactctctcagggccaggcggtgaagggcaat"+
				"cagctgttgccagtctcactggtgaaaagaaaaaccaccctggcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcat"+
				"taatgcagctggcacgacaggtttcccgactggaaagcgggcagtgataa";
		Sequence lacI_cds_seq = doc.createSequence("lacI_sequence", obj_ver, lacI_cds_sequences, Sequence.IUPAC_DNA);
		lacI_cds.addSequence(lacI_cds_seq);
		
		//creating output Sequences 
		String yfp_sequence = "CTGAAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAATACTAGAGAAAGAGGGGAAATACTAGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACAGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCAATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTTAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC";
		Sequence yfp_seq= doc.createSequence("YFP_protein_sequence", obj_ver, yfp_sequence, Sequence.IUPAC_DNA);
		yfp_cds.addSequence(yfp_seq);

		String rfp_sequence = "CTGAAGTGGTCGTGATCTGAAACTCGATCACCTGATGAGCTCAAGGCAGAGCGAAACCACCTCTACAAATAATTTTGTTTAATACTAGAGTCACACAGGAAAGTACTAGATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAAGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAATAACAGATAAAAAAAATCCTTAGCTTTCGCTAAGGATGATTTCT";
		Sequence rfp_seq= doc.createSequence("RFP_protein_sequence", obj_ver, rfp_sequence, Sequence.IUPAC_DNA);
		rfp_cds.addSequence(rfp_seq);

		String bfp_sequence = "CTGAAGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence bfp_seq= doc.createSequence("BFP_protein_sequence", obj_ver, bfp_sequence, Sequence.IUPAC_DNA);
		bfp_cds.addSequence(bfp_seq);

		String gfp_sequence = "CTGAAGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence gfp_seq= doc.createSequence("GFP_protein_sequence", obj_ver, gfp_sequence, Sequence.IUPAC_DNA);
		gfp_cds.addSequence(gfp_seq);
	}
	
	private static FunctionalComponent createFunctionalComponent(ModuleDefinition md, String fcId, ComponentDefinition cd) throws SBOLValidationException
	{
		return md.createFunctionalComponent(fcId, AccessType.PUBLIC, cd.getIdentity(), DirectionType.INOUT);
	}
	
	private static Interaction createDegradationInteraction(ModuleDefinition md, FunctionalComponent r1) throws SBOLValidationException
	{
		String interactionId = r1.getDisplayId() + "_degradation_interaction";
		String r1_id = r1.getDisplayId() + "_react_part";
		Interaction interaction = md.createInteraction(interactionId, SystemsBiologyOntology.DEGRADATION);
		interaction.createParticipation(r1_id, r1.getIdentity(),  SystemsBiologyOntology.REACTANT);
		return interaction;
	}
	
	private static Interaction createComplexInteraction(ModuleDefinition md, FunctionalComponent r1, FunctionalComponent r2, FunctionalComponent complex) throws SBOLValidationException
	{
		String interactionId = r1.getDisplayId() + "_" + r2.getDisplayId()  + "_complex_interaction";
		String r1_id = r1.getDisplayId() + "_" + r2.getDisplayId() + "_react_part";
		String r2_id = r2.getDisplayId() + "_" + r1.getDisplayId() + "_react_part";

		String complexId = r1.getDisplayId() + "_" + r2.getDisplayId()  + "_comp_part";
		Interaction interaction = md.createInteraction(interactionId, SystemsBiologyOntology.NON_COVALENT_BINDING);
		interaction.createParticipation(r1_id, r1.getIdentity(), SystemsBiologyOntology.REACTANT);
		interaction.createParticipation(r2_id, r2.getIdentity(), SystemsBiologyOntology.REACTANT);
		interaction.createParticipation(complexId, complex.getIdentity(), SystemsBiologyOntology.PRODUCT);
		return interaction;
	}

	private static Interaction createInhibitionInteraction(ModuleDefinition md, String regulationId, FunctionalComponent inhibitor, FunctionalComponent inhibited) throws SBOLValidationException
	{
		String interactionId = regulationId + "_interaction";
		String templateId = regulationId + "_prot_rep_part";
		String productId = regulationId + "_prom_part";
		Interaction interaction = md.createInteraction(interactionId, SystemsBiologyOntology.INHIBITION);
		interaction.createParticipation(templateId, inhibitor.getIdentity(), SystemsBiologyOntology.INHIBITOR);
		interaction.createParticipation(productId, inhibited.getIdentity(), SystemsBiologyOntology.INHIBITED);

		return interaction;
	}
	
	private static Interaction createProductionInteraction(ModuleDefinition md, String proteinId, FunctionalComponent template, FunctionalComponent product) throws SBOLValidationException
	{
		String interactionId = proteinId + "_protein_interaction";
		String templateId = proteinId + "_cds_part";
		String productId = proteinId + "_prod_part";
		Interaction interaction = md.createInteraction(interactionId, SystemsBiologyOntology.GENETIC_PRODUCTION);
		interaction.createParticipation(templateId, template.getIdentity(), SystemsBiologyOntology.TEMPLATE);
		interaction.createParticipation(productId, product.getIdentity(), SystemsBiologyOntology.PRODUCT);

		return interaction;
	}
	
	private static ComponentDefinition createCDS(SBOLDocument doc, String display) throws SBOLValidationException
	{
		return createComponentDefinition(doc, display, ComponentDefinition.DNA, SequenceOntology.CDS);
	}
	
	private static ComponentDefinition createPromoter(SBOLDocument doc, String display) throws SBOLValidationException
	{
		return createComponentDefinition(doc, display, ComponentDefinition.DNA, SequenceOntology.PROMOTER);
	}
	
	private static ComponentDefinition createComponentDefinition(SBOLDocument doc, String display, URI type, URI role) throws SBOLValidationException
	{
		ComponentDefinition cds = doc.createComponentDefinition(display, obj_ver, type);
		
		if(role == null)
		{
			//cds.addRole(SequenceOntology.CDS);
		}
		else
		{
			cds.addRole(role);
		}
		return cds;
	}
		
	private static URI getRole(String type) {
		if (type.equals("ribozyme")) {
	        return URI.create(so + "SO:0000374");
	    }
	    else if (type.equals("scar")) {
	        return URI.create(so + "SO:0001953");
	    }
	    else if (type.equals("cds")) {
	        return URI.create(so + "SO:0000316");
	    }
	    else if (type.equals("promoter")) {
	        return URI.create(so + "SO:0000167");
	    }
	    else if (type.equals("rbs")) {
	        return URI.create(so + "SO:0000139");
	    }
	    else if (type.equals("terminator")) {
	        return URI.create(so + "SO:0000141");
	    } else {
	        System.err.println("Part Type not found");
	        return null;
	    }
	}
	
	private static void convertPartsToSBOL(SBOLDocument document,HashMap<String,JSONObject> partsMap) throws SBOLValidationException {
		for (JSONObject part : partsMap.values()) {
			String name = (String)part.get("name");
			String dnasequence = (String)part.get("dnasequence");
			Sequence sequence = document.createSequence(name + "_sequence", version, dnasequence, Sequence.IUPAC_DNA);
			sequence.setName(name+"_sequence");
			sequence.addWasDerivedFrom(derivedFrom);
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(name, version, ComponentDefinition.DNA);
			componentDefinition.setName(name);
			componentDefinition.addWasDerivedFrom(derivedFrom);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			String partType = (String)part.get("type");
			componentDefinition.addRole(getRole(partType));
			componentDefinition.addSequence(sequence);
			
			if (partType.equals("cds")) {
				ComponentDefinition proteinComponentDefinition =
						document.createComponentDefinition(name+"_protein", version, ComponentDefinition.PROTEIN);
				proteinComponentDefinition.setName(name+"_protein");
				proteinComponentDefinition.addWasDerivedFrom(derivedFrom);
				proteinComponentDefinition.addWasGeneratedBy(activityURI);
				proteinComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

				ModuleDefinition moduleDefinition = 
						document.createModuleDefinition(name+"_protein_production", version);
				moduleDefinition.setName(name+"_protein_production");
				moduleDefinition.addWasDerivedFrom(derivedFrom);
				moduleDefinition.addWasGeneratedBy(activityURI);
				moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
				moduleDefinition.createFunctionalComponent(name, AccessType.PUBLIC, 
						componentDefinition.getIdentity(), DirectionType.NONE);
				moduleDefinition.createFunctionalComponent(name+"_protein", AccessType.PUBLIC, 
						proteinComponentDefinition.getIdentity(), DirectionType.NONE);
				Interaction interaction = moduleDefinition.createInteraction(name+"_protein_interaction", 
						SystemsBiologyOntology.GENETIC_PRODUCTION);
				interaction.createParticipation(name, name, SystemsBiologyOntology.TEMPLATE);
				interaction.createParticipation(name+"_protein", name+"_protein", SystemsBiologyOntology.PRODUCT);
			}
		}
	}

	public static void convertGatePartsToSBOL(SBOLDocument document,HashSet<JSONObject> gate_partsArr,
			HashMap<String,JSONObject> gatesMap,HashMap<String,JSONObject> responseMap) throws SBOLValidationException {
		for (JSONObject gate : gate_partsArr) {
			String gate_name = (String)gate.get("gate_name");
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(gate_name, version, ComponentDefinition.DNA);
			componentDefinition.setName(gate_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasDerivedFrom(derivedFrom);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
	        componentDefinition.createAnnotation(new QName(celloNS,"family","cello"), 
	        		(String)gatesMap.get(gate_name).get("system"));
	        //componentDefinition.addUriAnnotation(regulatorSO, gatesMap[gpartName].regulator);
	        componentDefinition.createAnnotation(new QName(celloNS,"gate_type","cello"), 
	        		(String)gatesMap.get(gate_name).get("gate_type"));
	        componentDefinition.createAnnotation(new QName(celloNS,"group_name","cello"), 
	        		(String)gatesMap.get(gate_name).get("group_name"));
	        componentDefinition.createAnnotation(new QName(celloNS,"color_hexcode","cello"), 
	        		(String)gatesMap.get(gate_name).get("color_hexcode"));
	        componentDefinition.createAnnotation(new QName(celloNS,"response_function","cello"), 
	        		(String)responseMap.get(gate_name).get("equation"));
	        JSONArray parameters = (JSONArray)responseMap.get(gate_name).get("parameters");
	        for (Object obj : parameters) {
	        	String name = (String)((JSONObject)obj).get("name");
	        	componentDefinition.createAnnotation(new QName(celloNS,name,"cello"), 
	        			(Double)((JSONObject)obj).get("value"));
	        }
	        JSONArray variables = (JSONArray)responseMap.get(gate_name).get("variables");
	        for (Object obj : variables) {
	        	String name = (String)((JSONObject)obj).get("name");
	        	componentDefinition.createAnnotation(new QName(celloNS,name+"_off_threshold","cello"), 
	        			(Double)((JSONObject)obj).get("off_threshold"));
	        	componentDefinition.createAnnotation(new QName(celloNS,name+"_on_threshold","cello"), 
	        			(Double)((JSONObject)obj).get("on_threshold"));
	        }
			
			JSONArray expression_cassettes = (JSONArray) gate.get("expression_cassettes");
			String seq = "";
			for (Object obj : expression_cassettes) {
				int annotationCount = 0;
				int start = 1;
				//int constraintCount = 0;
				//Component previousComponent = null;
		        Component currentComponent = null;
		        
				JSONObject expression_cassette = (JSONObject) obj;
				JSONArray cassette_parts = (JSONArray)expression_cassette.get("cassette_parts");
				for (Object obj2 : cassette_parts) {
					String partId = (String)obj2;
					ComponentDefinition partComponentDefinition = document.getComponentDefinition(partId, version);
					String cass_seq = document.getSequence(partId+"_sequence",version).getElements();
					seq += cass_seq;
					currentComponent = componentDefinition.createComponent(partId, AccessType.PUBLIC, partId, version);
//					if (previousComponent != null) {
//						componentDefinition.createSequenceConstraint("constraint"+constraintCount, RestrictionType.PRECEDES,
//								previousComponent.getIdentity(), currentComponent.getIdentity());
//						constraintCount++;
//					}
//					previousComponent = currentComponent;
					SequenceAnnotation sa = componentDefinition.createSequenceAnnotation("annotation"+annotationCount, 
							"range", start, start + cass_seq.length() - 1, OrientationType.INLINE);
					sa.setComponent(partId);
					start += cass_seq.length();
					annotationCount++;
					
					if (partComponentDefinition.getRoles().contains(SequenceOntology.CDS)) {
						String promoter = (String)gate.get("promoter");
						if (document.getModuleDefinition(partId+"_"+promoter+"_repression", version)==null) {
							ModuleDefinition moduleDefinition = 
									document.createModuleDefinition(partId+"_"+promoter+"_repression", version);
							moduleDefinition.addWasDerivedFrom(derivedFrom);
							moduleDefinition.addWasGeneratedBy(activityURI);
							moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
							Interaction interaction = moduleDefinition.createInteraction(partId+"_"+promoter+"_repression", 
									SystemsBiologyOntology.INHIBITION);
							interaction.createParticipation(partId+"_protein_participation", partId+"_protein", 
									SystemsBiologyOntology.INHIBITOR);
							interaction.createParticipation(promoter+"_promoter_participation", promoter, 
									SystemsBiologyOntology.INHIBITED);
						}
					}
				}
				
			}
			
			Sequence sequence = document.createSequence(gate_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(gate_name+"_sequence");
			sequence.addWasDerivedFrom(derivedFrom);
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}
	
	// args[0] - login email
	// args[1] - password
	// args[2] - login user
	// args[3] - temporary directory
	// args[4] - databasePrefix
	// args[5] - path to UCF file
	public static void main( String[] args ) throws SBOLValidationException, SBOLConversionException, SynBioHubException, FileNotFoundException, IOException, ParseException, URISyntaxException
    {
		if (args.length < 6) {
			System.err.println("Usage:");
			System.err.println(" login email");
			System.err.println(" password");
			System.err.println(" login user");
			System.err.println(" temporary directory");
			System.err.println(" database prefix");
			System.err.println(" path to UCF file");
			return;
		}
		// Create an SBOLDocument
		String databasePrefix = args[4];
		
		SBOLDocument document = new SBOLDocument(); 
		document.setDefaultURIprefix(uriPrefix); 
		document.setComplete(true); 
		document.setCreateDefaults(true);
		
		// Create an Activity
		Activity genericTopLevel = document.createActivity("cello2sbol", version);
		genericTopLevel.setName("Cello UCF to SBOL conversion");
		genericTopLevel.setDescription("Conversion of the Cello UCF parts and metadata to SBOL2");
		genericTopLevel.createAnnotation(new QName(dcNS,"creator","dc"), "Prashant Vaidyanathan");
		genericTopLevel.createAnnotation(new QName(dcNS,"creator","dc"), "Chris J. Myers");
		TimeZone tz = TimeZone.getTimeZone("UTC");
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
		df.setTimeZone(tz);
		createdDate = df.format(new Date());
		//System.out.println(createdDate);
		//DateTimeFormatter fmt = ISODateTimeFormat.dateTime();
		//DateTime endedAtTime = fmt.parseDateTime(createdDate);
		//System.out.println(endedAtTime);
		genericTopLevel.createAnnotation(new QName(provNS,"endedAtTime","prov"), createdDate);
		activityURI = genericTopLevel.getIdentity();
		//document.write(System.out);
		
		// Parse JSON
		HashMap<String,JSONObject> partsMap = new HashMap<String,JSONObject>();
		HashSet<JSONObject> gate_partsArr = new HashSet<JSONObject>();
		HashMap<String,JSONObject> gatesMap = new HashMap<String,JSONObject>();
		HashMap<String,JSONObject> responseMap = new HashMap<String,JSONObject>();

		JSONParser parser = new JSONParser();
		JSONArray a = (JSONArray) parser.parse(new FileReader(args[5]));

		for (Object o : a)
		{
			JSONObject ucf = (JSONObject) o;

			String collection = (String) ucf.get("collection");

			if (collection.equals("parts")) {
				partsMap.put((String)ucf.get("name"),ucf);
			}
			else if (collection.equals("gate_parts")) {
				gate_partsArr.add(ucf);
			}
			else if (collection.equals("gates")) {
				gatesMap.put((String)ucf.get("gate_name"),ucf);
			}
			else if (collection.equals("response_functions")) {
				responseMap.put((String)ucf.get("gate_name"),ucf);
			}
		}

		convertPartsToSBOL(document,partsMap);
        convertGatePartsToSBOL(document,gate_partsArr,gatesMap,responseMap);
        createSensorsReporters(document);
        
//        Collection cdsCollection = document.createCollection("cdsCollection", version);
//        for (ComponentDefinition cd : document.getComponentDefinitions()) {
//        	if (cd.containsRole(SequenceOntology.CDS)) {
//        		cdsCollection.addMember(cd.getIdentity());
//        	}
//        }
        
        // Validate
        SBOLValidate.validateSBOL(document,true,true,true);
        if (SBOLValidate.getNumErrors()>0) {
        	for (String error : SBOLValidate.getErrors()) {
        		System.out.println(error);
        	}
        } else {   
        	// Upload to SynBioHub
        	SynBioHubFrontend sbh = new SynBioHubFrontend(databasePrefix);
        	sbh.login(args[0], args[1]);
        	sbh.createCollection("CelloParts", "1", "Cello Parts", "These are the Cello parts", "27034378", true, document);
        	JSONArray motif_library = new JSONArray();
    		for (Object o : a)
    		{
    			JSONObject ucf = (JSONObject) o;

    			String collection = (String) ucf.get("collection");

    			if (collection.equals("gate_toxicity")) {
    				String gateName = (String)ucf.get("gate_name");
    				File file = new File(args[3] + gateName+"_gate_toxicity.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/"+gateName+"/1"), 
    						args[3] + gateName+"_gate_toxicity.json");
    			} else if (collection.equals("gate_cytometry")) {
    				String gateName = (String)ucf.get("gate_name");
    				File file = new File("/Users/myers/Downloads/"+gateName+"_gate_cytometry.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/myers/CelloParts/"+gateName+"/1"), 
    						args[3] + gateName+"_gate_cytometry.json");
    			} else if (collection.equals("header")) {
    				File file = new File(args[3] + "header.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
    						args[3] + "header.json");
    			} else if (collection.equals("measurement_std")) {
    				File file = new File(args[3] + "measurement_std.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
    						args[3] + "measurement_std.json");
    			} else if (collection.equals("logic_constraints")) {
    				File file = new File(args[3] + "logic_constraints.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
    						args[3] + "logic_constraints.json");
    			} else if (collection.equals("eugene_rules")) {
    				File file = new File(args[3] + "eugene_rules.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
    						args[3] + "eugene_rules.json");
    			} else if (collection.equals("genetic_locations")) {
    				File file = new File(args[3] + "genetic_locations.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
    						args[3] + "genetic_locations.json");
    			} else if (collection.equals("motif_library")) {
    				motif_library.add(ucf);
    			} else if (collection.equals("gates")) {
    			} else if (collection.equals("gate_parts")) {
    			} else if (collection.equals("parts")) {
    			} else if (collection.equals("response_functions")) {
    			} else {
        			System.out.println(collection);
    			}
    		}
			File file = new File(args[3] + "motif_library.json");
			FileOutputStream stream = new FileOutputStream(file);
			BufferedOutputStream buffer = new BufferedOutputStream(stream);
			stream.write(motif_library.toJSONString().getBytes());
			stream.close();
			buffer.close();
			sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/CelloParts_collection/1"), 
					args[3] + "motif_library.json");
        	System.out.println("Conversion, validation, and upload successful");
        }

        //document.write(System.out);
    }
}
