package UCF2SBOL.UCF2SBOL;

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

import org.joda.time.DateTime;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.Activity;
import org.sbolstandard.core2.Agent;
import org.sbolstandard.core2.Attachment;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.DirectionType;
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
//	static URI derivedFrom = URI.create("https://github.com/CIDARLAB/cello/blob/master/resources/UCF/Eco1C1G1T0.UCF.json");
	static String so = "http://identifiers.org/so/";
	static String provNS = "http://www.w3.org/ns/prov#";
	static String dcNS = "http://purl.org/dc/elements/1.1/";
	static String dcTermsNS = "http://purl.org/dc/terms/";
	static String celloNS = "http://cellocad.org/Terms/cello#";

	static URI activityURI;
	static String createdDate;
	
	private static void createSensor(SBOLDocument doc,String id,ComponentDefinition prom,
			ComponentDefinition riboJ,ComponentDefinition rbs,ComponentDefinition cds,ComponentDefinition term) 
					throws SBOLValidationException {
		ComponentDefinition Sensor = doc.createComponentDefinition(id, version, ComponentDefinition.DNA_REGION);
		Sensor.addRole(SequenceOntology.ENGINEERED_REGION);

		Sensor.createComponent(prom.getDisplayId(), AccessType.PRIVATE, prom.getIdentity());
		int start = 1;
		int end = prom.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length();
		Sensor.createSequenceAnnotation("annot1", "range", start, end, OrientationType.INLINE).setComponent(prom.getDisplayId());

		Sensor.createComponent(riboJ.getDisplayId(), AccessType.PRIVATE, riboJ.getIdentity());
		start = end+1;
		end = start+riboJ.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		Sensor.createSequenceAnnotation("annot2", "range", start, end, OrientationType.INLINE).setComponent(riboJ.getDisplayId());

		Sensor.createComponent(rbs.getDisplayId(), AccessType.PRIVATE, rbs.getIdentity());
		start = end+1;
		end = start+rbs.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		Sensor.createSequenceAnnotation("annot3", "range", start, end, OrientationType.INLINE).setComponent(rbs.getDisplayId());

		Sensor.createComponent(cds.getDisplayId(), AccessType.PRIVATE, cds.getIdentity());
		start = end+1;
		end = start+cds.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		Sensor.createSequenceAnnotation("annot4", "range", start, end, OrientationType.INLINE).setComponent(cds.getDisplayId());

		Sensor.createComponent(term.getDisplayId(), AccessType.PRIVATE, term.getIdentity());
		start = end+1;
		end = start+term.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		Sensor.createSequenceAnnotation("annot5", "range", start, end, OrientationType.INLINE).setComponent(term.getDisplayId());

		Sequence Sensor_seq = doc.createSequence(id+"_seq", version, 
				Sensor.getImpliedNucleicAcidSequence(),Sequence.IUPAC_DNA);
		Sensor.addSequence(Sensor_seq);
	}
    
    private static void createSensorsReporters(SBOLDocument doc) throws URISyntaxException,
    	SBOLValidationException, SBOLConversionException, IOException 
	{ 
		
		// YFP Reporter
//		ComponentDefinition ribozyme_scar = doc.createComponentDefinition("Ribozyme_scar", version, ComponentDefinition.DNA_REGION);
//		ribozyme_scar.addRole(getRole("scar"));
//		String ribozyme_scar_sequence = "CTGA";
//		Sequence ribozyme_scar_seq = doc.createSequence("Ribozyme_scar_sequence", ribozyme_scar_sequence, Sequence.IUPAC_DNA);
//		ribozyme_scar.addSequence(ribozyme_scar_seq);
		
		ComponentDefinition riboJ = doc.createComponentDefinition("RiboJ", version, ComponentDefinition.DNA_REGION);
		riboJ.addRole(getRole("ribozyme"));
		String riboJ_sequence = "AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA";
		Sequence riboJ_seq = doc.createSequence("RiboJ_sequence", riboJ_sequence, Sequence.IUPAC_DNA);
		riboJ.addSequence(riboJ_seq);
		
//		ComponentDefinition RBS_upstream_scar = doc.createComponentDefinition("RBS_upstream_scar", version, ComponentDefinition.DNA_REGION);
//		RBS_upstream_scar.addRole(getRole("scar"));
//		String RBS_upstream_scar_sequence = "TACTAGAG";
//		Sequence RBS_upstream_scar_seq = doc.createSequence("RBS_upstream_scar_sequence", RBS_upstream_scar_sequence, Sequence.IUPAC_DNA);
//		RBS_upstream_scar.addSequence(RBS_upstream_scar_seq);

		ComponentDefinition BBa_B0064_rbs = doc.createComponentDefinition("BBa_B0064_rbs", version, ComponentDefinition.DNA_REGION);
		BBa_B0064_rbs.addRole(getRole("rbs"));
		String BBa_B0064_rbs_sequence = "AAAGAGGGGAAA";
		Sequence BBa_B0064_rbs_seq = doc.createSequence("BBa_B0064_rbs_sequence", BBa_B0064_rbs_sequence, Sequence.IUPAC_DNA);
		BBa_B0064_rbs.addSequence(BBa_B0064_rbs_seq);
		
//		ComponentDefinition RBS_CDS_scar = doc.createComponentDefinition("RBS_CDS_scar", version, ComponentDefinition.DNA_REGION);
//		RBS_CDS_scar.addRole(getRole("scar"));
//		String RBS_CDS_scar_sequence = "TACTAG";
//		Sequence RBS_CDS_scar_seq = doc.createSequence("RBS_CDS_scar_sequence", RBS_CDS_scar_sequence, Sequence.IUPAC_DNA);
//		RBS_CDS_scar.addSequence(RBS_CDS_scar_seq);
		
		ComponentDefinition yfp_cds = createCDS(doc, "YFP");
		createProtein(doc,"YFP",yfp_cds);
		String yfp_sequence = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACAGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCAATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTTAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA";
		Sequence yfp_seq= doc.createSequence("YFP_protein_sequence", version, yfp_sequence, Sequence.IUPAC_DNA);
		yfp_cds.addSequence(yfp_seq);
		
		ComponentDefinition L3S2P21 = doc.createComponentDefinition("L3S2P21_terminator", version, ComponentDefinition.DNA_REGION);
		L3S2P21.addRole(getRole("terminator"));
		String L3S2P21_sequence = "CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC";
		Sequence L3S2P21_seq = doc.createSequence("L3S2P21_terminator_sequence", L3S2P21_sequence, Sequence.IUPAC_DNA);
		L3S2P21.addSequence(L3S2P21_seq);

		ComponentDefinition YFP_Reporter = doc.createComponentDefinition("YFP_reporter", version, ComponentDefinition.DNA_REGION);
		YFP_Reporter.addRole(SequenceOntology.ENGINEERED_REGION);

//		YFP_Reporter.createComponent("Ribozyme_scar", AccessType.PRIVATE, ribozyme_scar.getIdentity());
//		int start = 1;
//		int end = ribozyme_scar.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length();
//		YFP_Reporter.createSequenceAnnotation("annot1", "range", start, end, OrientationType.INLINE).setComponent("Ribozyme_scar");

		YFP_Reporter.createComponent("RiboJ", AccessType.PRIVATE, riboJ.getIdentity());
//		start = end+1;
//		end = start+riboJ.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		int start = 1;
		int end = riboJ.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length();
		YFP_Reporter.createSequenceAnnotation("annot2", "range", start, end, OrientationType.INLINE).setComponent("RiboJ");

//		YFP_Reporter.createComponent("RBS_upstream_scar", AccessType.PRIVATE, RBS_upstream_scar.getIdentity());
//		start = end+1;
//		end = start+RBS_upstream_scar.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
//		YFP_Reporter.createSequenceAnnotation("annot3", "range", start, end, OrientationType.INLINE).setComponent("RBS_upstream_scar");

		YFP_Reporter.createComponent("BBa_B0064_rbs", AccessType.PRIVATE, BBa_B0064_rbs.getIdentity());
		start = end+1;
		end = start+BBa_B0064_rbs.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		YFP_Reporter.createSequenceAnnotation("annot4", "range", start, end, OrientationType.INLINE).setComponent("BBa_B0064_rbs");

//		YFP_Reporter.createComponent("RBS_CDS_scar", AccessType.PRIVATE, RBS_CDS_scar.getIdentity());
//		start = end+1;
//		end = start+RBS_CDS_scar.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
//		YFP_Reporter.createSequenceAnnotation("annot5", "range", start, end, OrientationType.INLINE).setComponent("RBS_CDS_scar");

		YFP_Reporter.createComponent("YFP", AccessType.PRIVATE, yfp_cds.getIdentity());
		start = end+1;
		end = start+yfp_cds.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		YFP_Reporter.createSequenceAnnotation("annot6", "range", start, end, OrientationType.INLINE).setComponent("YFP");

		YFP_Reporter.createComponent("L3S2P21", AccessType.PRIVATE, L3S2P21.getIdentity());
		start = end+1;
		end = start+L3S2P21.getSequenceByEncoding(Sequence.IUPAC_DNA).getElements().length()-1;
		YFP_Reporter.createSequenceAnnotation("annot7", "range", start, end, OrientationType.INLINE).setComponent("L3S2P21");

		Sequence YRP_Reporter_seq = doc.createSequence("YFP_reporter_seq", version, 
				YFP_Reporter.getImpliedNucleicAcidSequence(),Sequence.IUPAC_DNA);
		YFP_Reporter.addSequence(YRP_Reporter_seq);
		
    	// Constitutive Promoter
		ComponentDefinition pCONST_prom = createPromoter(doc, "pCONST");
		String pCONST_prom_sequences = "GATAAGTCCCTAACTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC";
		Sequence pCONST_prom_seq= doc.createSequence("pCONST_sequence", version, pCONST_prom_sequences, Sequence.IUPAC_DNA);
		pCONST_prom.addSequence(pCONST_prom_seq);

		// LacI Sensor
		ComponentDefinition pTac_prom = createPromoter(doc, "pTac");
		ComponentDefinition lacI_cds = createCDS(doc, "LacI");
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
		Sequence lacI_cds_seq = doc.createSequence("LacI_sequence", version, lacI_cds_sequences, Sequence.IUPAC_DNA);
		lacI_cds.addSequence(lacI_cds_seq);
		createProtein(doc,"LacI",lacI_cds);
		doc.createComponentDefinition("IPTG", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"IPTG","LacI_protein");
		createInhibition(doc,"LacI_protein","pTac",0.0034,2.8,null,null);
		String pTac_prom_sequences = "AACGATCGTTGGCTGTGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGCTCACAATT";
		Sequence pTac_prom_seq= doc.createSequence("pTac_sequence", version, pTac_prom_sequences, Sequence.IUPAC_DNA);
		pTac_prom.addSequence(pTac_prom_seq);
		createSensor(doc,"LacI_sensor",pCONST_prom,riboJ,BBa_B0064_rbs,lacI_cds,L3S2P21);

		// TetR Sensor
		ComponentDefinition pTet_prom = createPromoter(doc, "pTet");
		ComponentDefinition tetR_cds = createCDS(doc, "TetR");
		String tetR_cds_sequences = "atgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaatgaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgccttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctttactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggtttttcactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacctactactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggattagaaaaacaacttaaatgtgaaagtgggtcctaa";
		Sequence tetR_cds_seq = doc.createSequence("TetR_sequence", version, tetR_cds_sequences, Sequence.IUPAC_DNA);
		tetR_cds.addSequence(tetR_cds_seq);
		createProtein(doc,"TetR",tetR_cds);
		doc.createComponentDefinition("aTc", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"aTc","TetR_protein");
		createInhibition(doc,"TetR_protein","pTet",0.0013,4.4,null,null);
		String pTet_prom_sequences = "TACTCCACCGTTGGCTTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCAC";
		Sequence pTet_prom_seq= doc.createSequence("pTet_sequence", version, pTet_prom_sequences, Sequence.IUPAC_DNA);
		pTet_prom.addSequence(pTet_prom_seq);
		createSensor(doc,"TetR_sensor",pCONST_prom,riboJ,BBa_B0064_rbs,tetR_cds,L3S2P21);
		
		// AraC Sensor
		ComponentDefinition pBAD_prom = createPromoter(doc, "pBAD");
		ComponentDefinition araC_cds = createCDS(doc, "AraC");
		String araC_cds_sequences = "atggctgaagcgcaaaatgatcccctgctgccgggatactcgtttaatgcccatctggtggcgggtttaacgccgattgaggccaacggttatctcgatttttttatcgaccgaccgctgggaatgaaaggttatattctcaatctcaccattcgcggtcagggggtggtgaaaaatcagggacgagaatttgtttgccgaccgggtgatattttgctgttcccgccaggagagattcatcactacggtcgtcatccggaggctcgcgaatggtatcaccagtgggtttactttcgtccgcgcgcctactggcatgaatggcttaactggccgtcaatatttgccaatacggggttctttcgcccggatgaagcgcaccagccgcatttcagcgacctgtttgggcaaatcattaacgccgggcaaggggaagggcgctattcggagctgctggcgataaatctgcttgagcaattgttactgcggcgcatggaagcgattaacgagtcgctccatccaccgatggataatcgggtacgcgaggcttgtcagtacatcagcgatcacctggcagacagcaattttgatatcgccagcgtcgcacagcatgtttgcttgtcgccgtcgcgtctgtcacatcttttccgccagcagttagggattagcgtcttaagctggcgcgaggaccaacgtatcagccaggcgaagctgcttttgagcaccacccggatgcctatcgccaccgtcggtcgcaatgttggttttgacgatcaactctatttctcgcgggtatttaaaaaatgcaccggggccagcccgagcgagttccgtgccggttgtgaagaaaaagtgaatgatgtagccgtcaagttgtcataa"; 
		Sequence araC_cds_seq = doc.createSequence("AraC_sequence", version, araC_cds_sequences, Sequence.IUPAC_DNA);
		araC_cds.addSequence(araC_cds_seq);
		createProtein(doc,"AraC",araC_cds);
		doc.createComponentDefinition("Ara", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"Ara","AraC_protein");
		createActivation(doc,"Ara_AraC_protein","pBAD",0.0082,2.5,null,null);
		String pBAD_prom_sequences = "ACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGC";
		Sequence pBAD_prom_seq= doc.createSequence("pBAD_sequence", version, pBAD_prom_sequences, Sequence.IUPAC_DNA);
		pBAD_prom.addSequence(pBAD_prom_seq);
		createSensor(doc,"AraC_sensor",pCONST_prom,riboJ,BBa_B0064_rbs,araC_cds,L3S2P21);

		// LuxR Sensor
		ComponentDefinition pLuxStar_prom = createPromoter(doc, "pLuxStar");
		ComponentDefinition luxR_cds = createCDS(doc, "LuxR");
		String luxR_cds_sequences = "ATGAACATTAAAAATATAAATGCTAATGAGAAGATAATTGATAAAATTAAAACTTGTAATAATAATAAAG"+
				"ATATTAATCAATGTTTATCTGAAATAGCAAAGATAATACATTGTGAATATTACCTATTCGCTATTATCTA"+
				"TCCTCACTCAATAATTAAACCTGATGTTTCAATTATAGATAATTACCCTGAAAAATGGCGTAAATATTAT"+
				"GATGATGCCGGACTACTAGAATATGACCCTGTAGTCGATTACTCTAAGTCCCATCATTCACCAATTAATT"+
				"GGAACGTATTCGAAAAAAAAACAATAAAAAAAGAGTCTCCGAATGTAATAAAAGAAGCACAGGAATCGGG"+
				"ACTCATTACTGGATTTAGCTTTCCAATTCATACTGCAAGTAATGGTTTTGGAATGCTCAGTTTTGCTCAT"+
				"TCAGATAAAGATATTTATACTGACAGTTTATTTTTACACGCTAGTACAAATGTACCATTAATGCTTCCTT"+
				"CTTTAGTCGATAATTATCAAAAAATAAATACGACACGTAAAAAGTCAGATTCTATTTTAACAAAAAGAGA"+
				"AAAAGAATGCTTAGCGTGGGCGAGTGAAGGAAAAAGTACATGGGATATTTCAAAAATACTTGGCTGCAGT"+
				"GAGCGTACTGTCACTTTTCATTTAACCAATACTCAAATGAAACTCAATACAACTAACCGCTGCCAAAGTA"+
				"TTTCTAAAGCAATTTTAACTGGCGCCATTAATTGTCCATACCTTAAAAATTAA";
		Sequence luxR_cds_seq = doc.createSequence("LuxR_sequence", version, luxR_cds_sequences, Sequence.IUPAC_DNA);
		luxR_cds.addSequence(luxR_cds_seq);
		createProtein(doc,"LuxR",luxR_cds);
		doc.createComponentDefinition("HSL", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"HSL","LuxR_protein");
		createActivation(doc,"HSL_LuxR_protein","pLuxStar",0.025,0.31,null,null);
		String pLuxStar_prom_sequences = "ATAGCTTCTTACCGGACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTACTTTCGAATAAA";
		Sequence pLuxStar_prom_seq= doc.createSequence("pLuxStar_sequence", version, pLuxStar_prom_sequences, Sequence.IUPAC_DNA);
		pLuxStar_prom.addSequence(pLuxStar_prom_seq);
		createSensor(doc,"LuxR_sensor",pCONST_prom,riboJ,BBa_B0064_rbs,luxR_cds,L3S2P21);
		
		// TODO: likely not correct
		/*
		ComponentDefinition gfp_cds = createCDS(doc, "GFP");
		createProtein(doc,"GFP",gfp_cds);
		String gfp_sequence = "AGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence gfp_seq= doc.createSequence("GFP_protein_sequence", version, gfp_sequence, Sequence.IUPAC_DNA);
		gfp_cds.addSequence(gfp_seq);

		ComponentDefinition rfp_cds = createCDS(doc, "RFP");
		createProtein(doc,"RFP",rfp_cds);
		String rfp_sequence = "AGTGGTCGTGATCTGAAACTCGATCACCTGATGAGCTCAAGGCAGAGCGAAACCACCTCTACAAATAATTTTGTTTAATACTAGAGTCACACAGGAAAGTACTAGATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAAGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAATAACAGATAAAAAAAATCCTTAGCTTTCGCTAAGGATGATTTCT";
		Sequence rfp_seq= doc.createSequence("RFP_protein_sequence", version, rfp_sequence, Sequence.IUPAC_DNA);
		rfp_cds.addSequence(rfp_seq);
		
		ComponentDefinition bfp_cds = createCDS(doc, "BFP");
		createProtein(doc,"BFP",bfp_cds);
		String bfp_sequence = "AGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence bfp_seq= doc.createSequence("BFP_protein_sequence", version, bfp_sequence, Sequence.IUPAC_DNA);
		bfp_cds.addSequence(bfp_seq);
		
		ComponentDefinition sigmaK1FR_cds = createCDS(doc, "SigmaK1FR");
		createProtein(doc,"SigmaK1FR",sigmaK1FR_cds);
		String sigmaK1FR_sequence = "ATGAGCATCGCGGCGACCCTGGAGAACGATCTGGCGCGTCTGGAAAACGAAAACGCTCGTCTCGAAAAAGACATCGCGAACCTGGAACGTGACCTGGCGAAACTGGAGCGTGAAGAAGCGTACTTCGGAGGTTCAGGTGGTAAGAACACTGGTGAAATCTCTGAGAAAGTCAAGCTGGGCACTAAGGCACTGGCTGGTCAATGGCTGGCTTACGGTGTTACTCGCAGTGTGACTAAGCGTTCAGTCATGACGCTGGCTTACGGGTCCAAAGAGTTCGGCTTCCGTCAACAAGTGCTGGAAGATACCATTCAGCCAGCTATTGATTCCGGCAAGGGTCTGATGTTCACTCAGCCGAATCAGGCTGCTGGATACATGGCTAAGCTGATTTGGGAATCTGTGAGCGTGACGGTGGTAGCTGCGGTTGAAGCAATGAACTGGCTTAAGTCTGCTGCTAAGCTGCTGGCTGCTGAGGTCAAAGATAAGAAGACTGGAGAGATTCTTCGCAAGCGTTGCGCTGTGCATTGGGTAACTCCTGATGGTTTCCCTGTGTGGCAGGAATACAAGAAGCCTATTCAGACGCGCTTGAACCTGAGGTTCCTCGGTTCGTTCAACCTCCAGCCGACCGTCAACACCAACAAAGATAGCGAGATTGATGCACACAAACAGGAGTCTGGTATCGCTCCTAACTTTGTACACAGCCAAGACGGTAGCCACCTTCGTAAGACTGTAGTGTGGGCACACGAGAAGTACGGAATCGAATCTTTTGCACTGATTCACGACTCCTTCGGTACGATTCCGGCTGACGCTGCGAACCTGTTCAAAGCAGTGCGCGAAACTATGGTTGACACATATGAGTCTTGTGATGTACTGGCTGATTTCTACGACCAGTTCGCTGACCAGTTGCACGAGTCTCAATTGGACAAAATGCCAGCACTTCCGGCTAAAGGTAACTTGAACCTCCGTGACATCTTAGAGTCGGACTTCGCGTTCGCG";
		Sequence sigmaK1FR_seq= doc.createSequence("SigmaK1FR_protein_sequence", version, sigmaK1FR_sequence, Sequence.IUPAC_DNA);
		sigmaK1FR_cds.addSequence(sigmaK1FR_seq);
		
		ComponentDefinition sigmaT3_cds = createCDS(doc, "SigmaT3");
		createProtein(doc,"SigmaT3",sigmaT3_cds);
		String sigmaT3_sequence = "ATGAGCATCGCGGCGACCCTGGAGAACGATCTGGCGCGTCTGGAAAACGAAAACGCTCGTCTCGAAAAAGACATCGCGAACCTGGAACGTGACCTGGCGAAACTGGAGCGTGAAGAAGCGTACTTCGGAGGTTCAGGTGGTAAGAACACTGGTGAAATCTCTGAGAAAGTCAAGCTGGGCACTAAGGCACTGGCTGGTCAATGGCTGGCTTACGGTGTTACTCGCAGTGTGACTAAGCGTTCAGTCATGACGCTGGCTTACGGGTCCAAAGAGTTCGGCTTCCGTCAACAAGTGCTGGAAGATACCATTCAGCCAGCTATTGATTCCGGCAAGGGTCTGATGTTCACTCAGCCGAATCAGGCTGCTGGATACATGGCTAAGCTGATTTGGGAATCTGTGAGCGTGACGGTGGTAGCTGCGGTTGAAGCAATGAACTGGCTTAAGTCTGCTGCTAAGCTGCTGGCTGCTGAGGTCAAAGATAAGAAGACTGGAGAGATTCTTCGCAAGCGTTGCGCTGTGCATTGGGTAACTCCTGATGGTTTCCCTGTGTGGCAGGAATACAAGAAGCCTATTCAGAAGCGCCTGGACATGATTTTCTTGGGTCAATTTCGCTTGCAACCTACCATTAACACCAACAAAGATAGCGAGATTGATGCACACAAACAGGAGTCTGGTATCGCTCCTAACTTTGTACACAGCCAAGACGGTAGCCACCTTCGTAAGACTGTAGTGTGGGCACACGAGAAGTACGGAATCGAATCTTTTGCACTGATTCACGACTCCTTCGGTACGATTCCGGCTGACGCTGCGAACCTGTTCAAAGCAGTGCGCGAAACTATGGTTGACACATATGAGTCTTGTGATGTACTGGCTGATTTCTACGACCAGTTCGCTGACCAGTTGCACGAGTCTCAATTGGACAAAATGCCAGCACTTCCGGCTAAAGGTAACTTGAACCTCCGTGACATCTTAGAGTCGGACTTCGCGTTCGCG";
		Sequence sigmaT3_seq= doc.createSequence("SigmaT3_protein_sequence", version, sigmaT3_sequence, Sequence.IUPAC_DNA);
		sigmaT3_cds.addSequence(sigmaT3_seq);
		
		ComponentDefinition sigmaT7_cds = createCDS(doc, "SigmaT7");
		createProtein(doc,"SigmaT7",sigmaT7_cds);
		String sigmaT7_sequence = "ATGAGCATCGCGGCGACCCTGGAGAACGATCTGGCGCGTCTGGAAAACGAAAACGCTCGTCTCGAAAAAGACATCGCGAACCTGGAACGTGACCTGGCGAAACTGGAGCGTGAAGAAGCGTACTTCGGAGGTTCAGGTGGTAAGAACACTGGTGAAATCTCTGAGAAAGTCAAGCTGGGCACTAAGGCACTGGCTGGTCAATGGCTGGCTTACGGTGTTACTCGCAGTGTGACTAAGAGTTCAGTCATGACGCTGGCTTACGGGTCCAAAGAGTTCGGCTTCCGTCAACAAGTGCTGGAAGATACCATTCAGCCAGCTATTGATTCCGGCAAGGGTCTGATGTTCACTCAGCCGAATCAGGCTGCTGGATACATGGCTAAGCTGATTTGGGAATCTGTGAGCGTGACGGTGGTAGCTGCGGTTGAAGCAATGAACTGGCTTAAGTCTGCTGCTAAGCTGCTGGCTGCTGAGGTCAAAGATAAGAAGACTGGAGAGATTCTTCGCAAGCGTTGCGCTGTGCATTGGGTAACTCCTGATGGTTTCCCTGTGTGGCAGGAATACAAGAAGCCTATTCAGACGCGCTTGAACCTGATGTTCCTCGGTCAGTTCCGCTTACAGCCTACCATTAACACCAACAAAGATAGCGAGATTGATGCACACAAACAGGAGTCTGGTATCGCTCCTAACTTTGTACACAGCCAAGACGGTAGCCACCTTCGTAAGACTGTAGTGTGGGCACACGAGAAGTACGGAATCGAATCTTTTGCACTGATTCACGACTCCTTCGGTACGATTCCGGCTGACGCTGCGAACCTGTTCAAAGCAGTGCGCGAAACTATGGTTGACACATATGAGTCTTGTGATGTACTGGCTGATTTCTACGACCAGTTCGCTGACCAGTTGCACGAGTCTCAATTGGACAAAATGCCAGCACTTCCGGCTAAAGGTAACTTGAACCTCCGTGACATCTTAGAGTCGGACTTCGCGTTCGCG";
		Sequence sigmaT7_seq= doc.createSequence("SigmaT7_protein_sequence", version, sigmaT7_sequence, Sequence.IUPAC_DNA);
		sigmaT7_cds.addSequence(sigmaT7_seq);

		ComponentDefinition sigmaCGG_cds = createCDS(doc, "SigmaCGG");
		createProtein(doc,"SigmaCGG",sigmaT7_cds);
		String sigmaCGG_sequence = "ATGAGCATCGCGGCGACCCTGGAGAACGATCTGGCGCGTCTGGAAAACGAAAACGCTCGTCTCGAAAAAGACATCGCGAACCTGGAACGTGACCTGGCGAAACTGGAGCGTGAAGAAGCGTACTTCGGAGGTTCAGGTGGTAAGAACACTGGTGAAATCTCTGAGAAAGTCAAGCTGGGCACTAAGGCACTGGCTGGTCAATGGCTGGCTTACGGTGTTACTCGCAGTGTGACTAAGCGTTCAGTCATGACGCTGGCTTACGGGTCCAAAGAGTTCGGCTTCCGTCAACAAGTGCTGGAAGATACCATTCAGCCAGCTATTGATTCCGGCAAGGGTCTGATGTTCACTCAGCCGAATCAGGCTGCTGGATACATGGCTAAGCTGATTTGGGAATCTGTGAGCGTGACGGTGGTAGCTGCGGTTGAAGCAATGAACTGGCTTAAGTCTGCTGCTAAGCTGCTGGCTGCTGAGGTCAAAGATAAGAAGACTGGAGAGATTCTTCGCAAGCGTTGCGCTGTGCATTGGGTAACTCCTGATGGTTTCCCTGTGTGGCAGGAATACAAGAAGCCTATTAAAACGCGCGTGCATATTATGTTCCTCGGTCAGTTCGAAATGCAGCCTACCATTAACACCAACAAAGATAGCGAGATTGATGCACGCAAACAGGAGTCTGGTATCGCTCCTAACTTTGTACACAGCCAAGACGGTAGCCACCTTCGTAAGACTGTAGTGTGGGCACACGAGAAGTACGGAATCGAATCTTTTGCACTGATTCACGACTCCTTCGGTACGATTCCGGCTGACGCTGCGAACCTGTTCAAAGCAGTGCGCGAAACTATGGTTGACACATATGAGTCTTGTGATGTACTGGCTGATTTCTACGACCAGTTCGCTGACCAGTTGCACGAGTCTCAATTGGACAAAATGCCAGCACTTCCGGCTAAAGGTAACTTGAACCTCCGTGACATCTTAGAGTCGGACTTCGCGTTCGCG";
		Sequence sigmaCGG_seq= doc.createSequence("SigmaCGG_protein_sequence", version, sigmaCGG_sequence, Sequence.IUPAC_DNA);
		sigmaCGG_cds.addSequence(sigmaCGG_seq);
		*/
	}
	
	private static ComponentDefinition createCDS(SBOLDocument doc, String display) throws SBOLValidationException
	{
		return createComponentDefinition(doc, display, ComponentDefinition.DNA_REGION, SequenceOntology.CDS);
	}
	
	private static ComponentDefinition createPromoter(SBOLDocument doc, String display) throws SBOLValidationException
	{
		return createComponentDefinition(doc, display, ComponentDefinition.DNA_REGION, SequenceOntology.PROMOTER);
	}
	
	private static ComponentDefinition createComponentDefinition(SBOLDocument doc, String display, URI type, URI role) throws SBOLValidationException
	{
		ComponentDefinition cds = doc.createComponentDefinition(display, version, type);
		cds.setName(display);
		cds.addWasGeneratedBy(activityURI);
		cds.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		
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
	        return URI.create(so + "SO:0001977");
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
	    } else if (type.equals("grna")) {
	        return URI.create(so + "SO:0001264");
	    } else if (type.equals("UAS")) {
	        return URI.create(so + "SO:0001678");
	    } else if (type.equals("LinkingSequence")) {
	        return URI.create(so + "SO:0001678");
	    } else if (type.equals("TataBox")) {
	        return URI.create(so + "SO:0000174");
	    } else if (type.equals("tss")) {
	        return URI.create(so + "SO:0000315");
	    } else if (type.equals("Kozak")) {
	        return URI.create(so + "SO:0001647");
	    } else if (type.equals("Codon")) {
	        return URI.create(so + "SO:0000360");
	    } else if (type.equals("operator")) {
	        return URI.create(so + "SO:0000057");
	    } else if (type.equals("backbone")) {
	        return URI.create(so + "SO:0000755");
	    } else if (type.equals("spacer")) {
	    	return URI.create(so + "SO:0002223");
	    } else {
	        System.err.println("Part Type " + type + " not found");
	        return null;
	    }
	}

	private static void convertPartsToSBOL(SBOLDocument document,HashMap<String,JSONObject> partsMap) throws SBOLValidationException {
		for (JSONObject part : partsMap.values()) {
			String name = (String)part.get("name");
			name = name.replace("-", "_");
			System.out.println(name);
			String dnasequence = (String)part.get("dnasequence");
			Sequence sequence = document.createSequence(name + "_sequence", version, dnasequence, Sequence.IUPAC_DNA);
			sequence.setName(name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(name);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			String partType = (String)part.get("type");
			componentDefinition.addRole(getRole(partType));
			componentDefinition.addSequence(sequence);
			
			if (partType.equals("cds")) {
				createProtein(document,name,componentDefinition);
			}
			if (partType.equals("grna")) {
				createRNA(document,name,componentDefinition);
			}
		}
	}
	
	private static void createProtein(SBOLDocument document,String cdsId,ComponentDefinition cds) throws SBOLValidationException 
	{
		ComponentDefinition proteinComponentDefinition =
				document.createComponentDefinition(cdsId+"_protein", version, ComponentDefinition.PROTEIN);
		proteinComponentDefinition.setName(cdsId+"_protein");
		proteinComponentDefinition.addWasGeneratedBy(activityURI);
		proteinComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(cdsId+"_protein_production", version);
		moduleDefinition.setName(cdsId+"_protein_production");
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(cdsId, AccessType.PUBLIC, 
				cds.getIdentity(), DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(cdsId+"_protein", AccessType.PUBLIC, 
				proteinComponentDefinition.getIdentity(), DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(cdsId+"_protein_interaction", 
				SystemsBiologyOntology.GENETIC_PRODUCTION);
		interaction.createParticipation(cdsId, cdsId, SystemsBiologyOntology.TEMPLATE);
		interaction.createParticipation(cdsId+"_protein", cdsId+"_protein", SystemsBiologyOntology.PRODUCT);
		
		moduleDefinition = 
				document.createModuleDefinition(cdsId+"_protein_degradation", version);
		moduleDefinition.setName(cdsId+"_protein_degradation");
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(cdsId+"_protein", AccessType.PUBLIC, 
				proteinComponentDefinition.getIdentity(), DirectionType.NONE);
		String interactionId = cdsId + "_degradation_interaction";
		interaction = moduleDefinition.createInteraction(interactionId, SystemsBiologyOntology.DEGRADATION);
		interaction.createParticipation(cdsId+"_protein", cdsId+"_protein",  SystemsBiologyOntology.REACTANT);
	}
	
	
	private static void createRNA(SBOLDocument document,String rnaId,ComponentDefinition rna) throws SBOLValidationException 
	{
		ComponentDefinition rnaComponentDefinition =
				document.createComponentDefinition(rnaId+"_rna", version, ComponentDefinition.RNA_MOLECULE);
		rnaComponentDefinition.addRole(URI.create(so + "SO:0001998"));
		rnaComponentDefinition.setName(rnaId+"_rna");
		rnaComponentDefinition.addWasGeneratedBy(activityURI);
		rnaComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(rnaId+"_rna_production", version);
		moduleDefinition.setName(rnaId+"_rna_production");
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(rnaId, AccessType.PUBLIC, 
				rna.getIdentity(), DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(rnaId+"_protein", AccessType.PUBLIC, 
				rnaComponentDefinition.getIdentity(), DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(rnaId+"_rna_interaction", 
				SystemsBiologyOntology.GENETIC_PRODUCTION);
		interaction.createParticipation(rnaId, rnaId, SystemsBiologyOntology.TEMPLATE);
		interaction.createParticipation(rnaId+"_rna", rnaId+"_rna", SystemsBiologyOntology.PRODUCT);
		
		moduleDefinition = 
				document.createModuleDefinition(rnaId+"_rna_degradation", version);
		moduleDefinition.setName(rnaId+"_rna_degradation");
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(rnaId+"_rna", AccessType.PUBLIC, 
				rnaComponentDefinition.getIdentity(), DirectionType.NONE);
		String interactionId = rnaId + "_degradation_interaction";
		interaction = moduleDefinition.createInteraction(interactionId, SystemsBiologyOntology.DEGRADATION);
		interaction.createParticipation(rnaId+"_rna", rnaId+"_rna",  SystemsBiologyOntology.REACTANT);
	}
	
	private static void createInhibition(SBOLDocument document,String inhibitor,String inhibited,
			Double ymin,Double ymax,Double alpha,Double beta) throws SBOLValidationException 
	{
		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(inhibitor+"_"+inhibited+"_repression", version);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(inhibitor, AccessType.PUBLIC, 
				inhibitor, version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(inhibited, AccessType.PUBLIC, 
				inhibited, version, DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(inhibitor+"_"+inhibited+"_repression", 
				SystemsBiologyOntology.INHIBITION);
		interaction.createParticipation(inhibitor+"_participation", inhibitor, 
				SystemsBiologyOntology.INHIBITOR);
		interaction.createParticipation(inhibited+"_promoter_participation", inhibited, 
				SystemsBiologyOntology.INHIBITED);
		if (ymin != null) {
			interaction.createAnnotation(new QName(celloNS,"ymin","cello"),ymin); 
		}
		if (ymax != null) {
			interaction.createAnnotation(new QName(celloNS,"ymax","cello"),ymax); 
		}
		if (alpha != null) {
			interaction.createAnnotation(new QName(celloNS,"alpha","cello"),alpha); 
		}
		if (beta != null) {
			interaction.createAnnotation(new QName(celloNS,"beta","cello"),beta); 
		}
	}
	
	private static void createActivation(SBOLDocument document,String activator,String promoter,
			Double ymin,Double ymax,Double alpha,Double beta) throws SBOLValidationException 
	{
		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(activator+"_"+promoter+"_activation", version);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(activator, AccessType.PUBLIC, 
				activator, version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(promoter, AccessType.PUBLIC, 
				promoter, version, DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(activator+"_"+promoter+"_activation", 
				SystemsBiologyOntology.STIMULATION);
		interaction.createParticipation(activator+"_protein_participation", activator, 
				SystemsBiologyOntology.STIMULATOR);
		interaction.createParticipation(promoter+"_promoter_participation", promoter, 
				SystemsBiologyOntology.STIMULATED);
		if (ymin != null) {
			interaction.createAnnotation(new QName(celloNS,"ymin","cello"),ymin); 
		}
		if (ymax != null) {
			interaction.createAnnotation(new QName(celloNS,"ymax","cello"),ymax); 
		}
		if (alpha != null) {
			interaction.createAnnotation(new QName(celloNS,"alpha","cello"),alpha); 
		}
		if (beta != null) {
			interaction.createAnnotation(new QName(celloNS,"beta","cello"),beta); 
		}
	}
	
	private static void createComplex(SBOLDocument document,String reactant1,String reactant2) throws SBOLValidationException 
	{
		String complex = reactant1 + "_" + reactant2;
		ComponentDefinition complexComponentDefinition = 
				document.createComponentDefinition(complex, version, ComponentDefinition.COMPLEX);
		complexComponentDefinition.setName(complex);
		complexComponentDefinition.addWasGeneratedBy(activityURI);
		complexComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(complex+"_complex_formation", version);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(reactant1, AccessType.PUBLIC, 
				reactant1, version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(reactant2, AccessType.PUBLIC, 
				reactant2, version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(complex, AccessType.PUBLIC, 
				complex, version, DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(complex+"_complex_formation", 
				SystemsBiologyOntology.NON_COVALENT_BINDING);
		interaction.createParticipation(reactant1, reactant1, SystemsBiologyOntology.REACTANT);
		interaction.createParticipation(reactant2, reactant2, SystemsBiologyOntology.REACTANT);
		interaction.createParticipation(complex, complex, SystemsBiologyOntology.PRODUCT);
		
		moduleDefinition = 
				document.createModuleDefinition(complex+"_degradation", version);
		moduleDefinition.setName(complex+"_degradation");
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(complex, AccessType.PUBLIC, 
				complexComponentDefinition.getIdentity(), DirectionType.NONE);
		String interactionId = complex + "_degradation_interaction";
		interaction = moduleDefinition.createInteraction(interactionId, SystemsBiologyOntology.DEGRADATION);
		interaction.createParticipation(complex, complex,  SystemsBiologyOntology.REACTANT);
	}

	private static void convertGatePartsToSBOL(SBOLDocument document,HashSet<JSONObject> gate_partsArr,
			HashMap<String,JSONObject> gatesMap,HashMap<String,JSONObject> responseMap) throws SBOLValidationException {
		for (JSONObject gate : gate_partsArr) {
			String gate_name = (String)gate.get("gate_name");
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(gate_name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(gate_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
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
	        componentDefinition.createAnnotation(new QName(celloNS,"tandem_efficiency_factor","cello"), 
	        		(String)responseMap.get(gate_name).get("tandem_efficiency_factor"));
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
//		        Component currentComponent = null;
		        
				JSONObject expression_cassette = (JSONObject) obj;
				JSONArray cassette_parts = (JSONArray)expression_cassette.get("cassette_parts");
				for (Object obj2 : cassette_parts) {
					String partId = (String)obj2;
					ComponentDefinition partComponentDefinition = document.getComponentDefinition(partId, version);
					String cass_seq = document.getSequence(partId+"_sequence",version).getElements();
					seq += cass_seq;
					//currentComponent = 
					componentDefinition.createComponent(partId, AccessType.PUBLIC, partId, version);
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
						if (document.getModuleDefinition(partId+"_protein_"+promoter+"_repression", version)==null) {
							createInhibition(document,partId+"_protein",promoter,null,null,null,null);
						}
					}
					if (partComponentDefinition.getRoles().contains(URI.create(so + "SO:0001264"))) {
						String promoter = (String)gate.get("promoter");
						createComplex(document,partId+"_rna","dCAS9_Mxi1_protein");
						if (document.getModuleDefinition(partId+"_rna_dCAS9_Mxi1_protein_"+promoter+"_repression", version)==null) {
							createInhibition(document,partId+"_rna_dCAS9_Mxi1_protein",promoter,null,null,null,null);
						}
					}
				}
				
			}
			
			Sequence sequence = document.createSequence(gate_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(gate_name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}

	private static void convertInputSensorsToSBOL(SBOLDocument document,HashSet<JSONObject> input_sensorsArr) throws SBOLValidationException {
		for (JSONObject sensor : input_sensorsArr) {
			String sensor_name = (String)sensor.get("name");
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(sensor_name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(sensor_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
	        componentDefinition.createAnnotation(new QName(celloNS,"gateType","cello"), "input_sensor");
					        
			JSONArray parts = (JSONArray)sensor.get("parts");
			String seq = "";
			int annotationCount = 0;
			int start = 1;
			for (Object obj2 : parts) {
				String partId = (String)obj2;
				ComponentDefinition partComponentDefinition = document.getComponentDefinition(partId, version);
				String cass_seq = document.getSequence(partId+"_sequence",version).getElements();
				seq += cass_seq;
				//currentComponent = 
				componentDefinition.createComponent(partId, AccessType.PUBLIC, partId, version);
//						if (previousComponent != null) {
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
					String promoter = (String)sensor.get("promoter");
					String input_molecule = (String)sensor.get("input_molecule");
					Double signal_low = (Double)sensor.get("signal_low");
					Double signal_high = (Double)sensor.get("signal_high");
					Double alpha = null;
					Double beta = null;
			        JSONArray parameters = (JSONArray)sensor.get("parameters");
			        for (Object obj : parameters) {
			        	String name = (String)((JSONObject)obj).get("name");
			        	if (name.equals("signal_low")) {
			        		signal_low = (Double)((JSONObject)obj).get("value");
			        	} else if (name.equals("signal_high")) {
			        		signal_high = (Double)((JSONObject)obj).get("value");
			        	} else if (name.equals("alpha")) {
			        		alpha = (Double)((JSONObject)obj).get("value");
			        	} else if (name.equals("beta")) {
			        		beta = (Double)((JSONObject)obj).get("value");
			        	} 
			        }


					document.createComponentDefinition(input_molecule, version, ComponentDefinition.SMALL_MOLECULE);
					createComplex(document,input_molecule,partId+"_protein");
					if (((String)sensor.get("type")).equals("complex_stimulator")) {
						createActivation(document,input_molecule+"_"+partId+"_protein",promoter,signal_low,signal_high,alpha,beta);
					} else if (((String)sensor.get("type")).equals("sequester_inhibitor")) {
						createInhibition(document,partId+"_protein",promoter,signal_low,signal_high,alpha,beta);
					}
				}
				
			}
			
			Sequence sequence = document.createSequence(sensor_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(sensor_name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}

	private static void convertOutputReportersToSBOL(SBOLDocument document,HashSet<JSONObject> output_reportersArr) throws SBOLValidationException {
		for (JSONObject sensor : output_reportersArr) {
			String reporter_name = (String)sensor.get("name");
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(reporter_name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(reporter_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
	        componentDefinition.createAnnotation(new QName(celloNS,"gateType","cello"), "output_reporter");
					        
			JSONArray parts = (JSONArray)sensor.get("parts");
			String seq = "";
			int annotationCount = 0;
			int start = 1;
			for (Object obj2 : parts) {
				String partId = (String)obj2;
				document.getComponentDefinition(partId, version);
				String cass_seq = document.getSequence(partId+"_sequence",version).getElements();
				seq += cass_seq;
				//currentComponent = 
				componentDefinition.createComponent(partId, AccessType.PUBLIC, partId, version);
//						if (previousComponent != null) {
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
				
			}
			
			Sequence sequence = document.createSequence(reporter_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(reporter_name+"_sequence");
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
	// args[6] - databaseURL
	// args[7] - collection id
	// args[8] - collection version
	// args[9] - collection name
	// args[10] - collection description
	// args[11] - collection pubMedId
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
			System.err.println(" database URL");
			System.err.println(" collection id");
			System.err.println(" collection version");
			System.err.println(" collection name");
			System.err.println(" collection description");
			System.err.println(" collection pubMedId");
			return;
		}
		// Create an SBOLDocument
		String loginEmail = args[0];
		String password = args[1];
		String loginUser = args[2];
		String tmpDir = args[3];
		String databasePrefix = args[4];
		String pathToUCFFile = args[5];
		String databaseURL = databasePrefix;
		String collectionId = "Cello_Parts";
		String collectionVersion = "1";
		String collectionName = "Cello Parts";
		String collectionDescription = "These are the Cello parts";
		String collectionPubMedId = "27034378";
		if (args.length > 6) {
			databaseURL = args[6];
		}
		if (args.length > 7) {
			collectionId = args[7];
		}
		if (args.length > 8) {
			collectionVersion = args[8];
		}
		if (args.length > 9) {
			collectionName = args[9];
		}
		if (args.length > 10) {
			collectionDescription = args[10];
		}
		if (args.length > 11) {
			collectionPubMedId = args[11];
		}
		
		SBOLDocument document = new SBOLDocument(); 
		document.setDefaultURIprefix(uriPrefix); 
		document.setComplete(true); 
		document.setCreateDefaults(true);
		
		TimeZone tz = TimeZone.getTimeZone("UTC");
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
		df.setTimeZone(tz);
		createdDate = df.format(new Date());
		
		// Create an Activity
		Activity activity = document.createActivity("CelloUCF2sbol_Activity", version);
		activity.setName("Cello UCF to SBOL conversion");
		activity.setDescription("Conversion of the Cello UCF parts and metadata to SBOL 2 documents.");
		activity.setEndedAtTime(DateTime.now());
		activityURI = activity.getIdentity();
		Agent agent = document.createAgent("CelloUCF2SBOL", version);
		agent.setName("Cello UCF to SBOL");
		agent.setDescription("A script to convert Cello UCF parts and metadata to SBOL 2 documents.");
		agent.createAnnotation(new QName(dcNS,"source","dc"), URI.create("https://github.com/MyersResearchGroup/UCF2SBOL"));
		agent.createAnnotation(new QName(dcNS,"creator","dc"), "Prashant Vaidyanathan");
		agent.createAnnotation(new QName(dcNS,"creator","dc"), "Chris J. Myers");
		activity.createAssociation("association", agent.getIdentity());
		
//		Attachment attachmentUCF = document.createAttachment("Eco1C1G1T1_UCF", version, 
//				URI.create("https://github.com/MyersResearchGroup/UCF2SBOL/blob/master/UCF2SBOL/src/main/resources/Eco1C1G1T1.UCF.json"));
//		activity.createUsage("UCF_File", attachmentUCF.getIdentity());
//		document.write(System.out);
//		if (true) {
//			return;
//		}
		
		// Parse JSON
		HashMap<String,JSONObject> partsMap = new HashMap<String,JSONObject>();
		HashSet<JSONObject> gate_partsArr = new HashSet<JSONObject>();
		HashSet<JSONObject> input_sensorsArr = new HashSet<JSONObject>();
		HashSet<JSONObject> output_reportersArr = new HashSet<JSONObject>();
		HashMap<String,JSONObject> gatesMap = new HashMap<String,JSONObject>();
		HashMap<String,JSONObject> responseMap = new HashMap<String,JSONObject>();

		JSONParser parser = new JSONParser();
		JSONArray a = (JSONArray) parser.parse(new FileReader(pathToUCFFile));

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
			else if (collection.equals("input_sensors")) {
				input_sensorsArr.add(ucf);
			}
			else if (collection.equals("output_reporters")) {
				output_reportersArr.add(ucf);
			}
			else if (collection.equals("gates")) {
				gatesMap.put((String)ucf.get("gate_name"),ucf);
			}
			else if (collection.equals("response_functions")) {
				responseMap.put((String)ucf.get("gate_name"),ucf);
			}
		}
        
//		// dCAS9
//        ComponentDefinition dCas9 = createCDS(document,"dCAS9_Mxi1");
//        createProtein(document,"dCAS9_Mxi1",dCas9);
//        
//        // yeGFP
//		ComponentDefinition yegfp_cds = createCDS(document, "yeGFP");
//		createProtein(document,"yeGFP",yegfp_cds);
//		String yegfp_sequence = "ggattctagaactagtggatctacaaaatgtctaaaggtgaagaattattcactggtgttgtcccaattttggttgaattagatggtgatgttaatggtcacaaattttctgtctccggtgaaggtgaaggtgatgctacttacggtaaattgaccttaaaatttatttgtactactggtaaattgccagttccatggccaaccttagtcactactttcggttatggtgttcaatgttttgctagatacccagatcatatgaaacaacatgactttttcaagtctgccatgccagaaggttatgttcaagaaagaactatttttttcaaagatgacggtaactacaagaccagagctgaagtcaagtttgaaggtgataccttagttaatagaatcgaattaaaaggtattgattttaaagaagatggtaacattttaggtcacaaattggaatacaactataactctcacaatgtttacatcatggctgacaaacaaaagaatggtatcaaagttaacttcaaaattagacacaacattgaagatggttctgttcaattagctgaccattatcaacaaaatactccaattggtgatggtccagtcttgttaccagacaaccattacttatccactcaatctgccttatccaaagatccaaacgaaaagagagaccacatggtcttgttagaatttgttactgctgctggtattacccatggtatggatgaattgtacaaataatgataccgtcgacctcgagtc";
//		Sequence yegfp_seq= document.createSequence("yeGFP_protein_sequence", version, yegfp_sequence, Sequence.IUPAC_DNA);
//		yegfp_cds.addSequence(yegfp_seq);

		convertPartsToSBOL(document,partsMap);
        convertGatePartsToSBOL(document,gate_partsArr,gatesMap,responseMap);
        convertInputSensorsToSBOL(document,input_sensorsArr);
        convertOutputReportersToSBOL(document,output_reportersArr);
        
        //createSensorsReporters(document);
        
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
        	System.out.println("Uploading to SynBioHub " + databaseURL);
        	System.out.println("Database Prefix: " + databasePrefix + " Login " + loginUser + " CollectionId " + collectionId + " Version " + collectionVersion);
        	SynBioHubFrontend sbh = new SynBioHubFrontend(databaseURL,databasePrefix);
        	sbh.login(loginEmail, password);
        	sbh.createCollection(collectionId, collectionVersion, collectionName, collectionDescription,
        			collectionPubMedId, true);
        	sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), pathToUCFFile);
        	SBOLDocument doc = sbh.getSBOL(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion));
        	for (Attachment attachment : doc.getAttachments()) {
        		activity.createUsage("UCF_file", attachment.getIdentity());
        		break;
        	}
        	document.write("/Users/myers/"+collectionId + ".xml");
        	System.out.println(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion);
        	sbh.addToCollection(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), false, document);
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
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + gateName + "/" + collectionVersion), 
    						tmpDir + gateName+"_gate_toxicity.json");
    			} else if (collection.equals("gate_cytometry")) {
    				String gateName = (String)ucf.get("gate_name");
    				File file = new File(args[3] + gateName+"_gate_cytometry.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + gateName + "/" + collectionVersion), 
    						tmpDir + gateName+"_gate_cytometry.json");
    			} else if (collection.equals("header")) {
    				File file = new File(args[3] + "header.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "header.json");
    			} else if (collection.equals("measurement_std")) {
    				File file = new File(args[3] + "measurement_std.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "measurement_std.json");
    			} else if (collection.equals("logic_constraints")) {
    				File file = new File(args[3] + "logic_constraints.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "logic_constraints.json");
    			} else if (collection.equals("eugene_rules")) {
    				File file = new File(args[3] + "eugene_rules.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "eugene_rules.json");
    			} else if (collection.equals("genetic_locations")) {
    				File file = new File(args[3] + "genetic_locations.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "genetic_locations.json");
    			} else if (collection.equals("PartitionProfile")) {
    				File file = new File(args[3] + "PartitionProfile.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "PartitionProfile.json");
    			} else if (collection.equals("containers")) {
    				File file = new File(args[3] + "containers.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "containers.json");
    			} else if (collection.equals("motif_library")) {
    				motif_library.add(ucf);
    			} else if (collection.equals("gates")) {
    			} else if (collection.equals("gate_parts")) {
    			} else if (collection.equals("input_sensors")) {
    			} else if (collection.equals("output_reporters")) {
    			} else if (collection.equals("parts")) {
    			} else if (collection.equals("response_functions")) {
    			} else {
        			System.out.println(collection);
    			}
    		}
			File file = new File(tmpDir + "motif_library.json");
			FileOutputStream stream = new FileOutputStream(file);
			BufferedOutputStream buffer = new BufferedOutputStream(stream);
			stream.write(motif_library.toJSONString().getBytes());
			stream.close();
			buffer.close();
			sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
					tmpDir + "motif_library.json");
        	System.out.println("Conversion, validation, and upload successful");
        }

        //document.write(System.out);
    }
}
