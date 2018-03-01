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
    
    private static void createSensorsReporters(SBOLDocument doc) throws URISyntaxException,
    	SBOLValidationException, SBOLConversionException, IOException 
	{ 
		ComponentDefinition pTac_prom = createPromoter(doc, "pTac");
		ComponentDefinition lacI_cds = createCDS(doc, "LacI");
		createProtein(doc,"LacI",lacI_cds);
		doc.createComponentDefinition("IPTG", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"IPTG","LacI_protein");
		createInhibition(doc,"LacI","pTac");

		ComponentDefinition pTet_prom = createPromoter(doc, "pTet");
		ComponentDefinition tetR_cds = createCDS(doc, "TetR");
		createProtein(doc,"TetR",tetR_cds);
		doc.createComponentDefinition("aTc", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"aTc","TetR_protein");
		createInhibition(doc,"TetR","pTet");
		
		ComponentDefinition pBAD_prom = createPromoter(doc, "pBAD");
		ComponentDefinition araC_cds = createCDS(doc, "AraC");
		createProtein(doc,"AraC",araC_cds);
		doc.createComponentDefinition("Ara", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"Ara","AraC_protein");
		createActivation(doc,"Ara_AraC_protein","pBAD");
		
		//ComponentDefinition luxI_cds = createCDS(doc, "LuxI");
		//createProtein(doc,"LuxI",luxI_cds);
		
		ComponentDefinition pLuxStar_prom = createPromoter(doc, "pLuxStar");
		ComponentDefinition luxR_cds = createCDS(doc, "LuxR");
		createProtein(doc,"LuxR",luxR_cds);
		doc.createComponentDefinition("HSL", version, ComponentDefinition.SMALL_MOLECULE);
		createComplex(doc,"HSL","LuxR_protein");
		createActivation(doc,"HSL_LuxR_protein","pLuxStar");
		
		ComponentDefinition gfp_cds = createCDS(doc, "GFP");
		createProtein(doc,"GFP",gfp_cds);
		ComponentDefinition rfp_cds = createCDS(doc, "RFP");
		createProtein(doc,"RFP",rfp_cds);
		ComponentDefinition bfp_cds = createCDS(doc, "BFP");
		createProtein(doc,"BFP",bfp_cds);
		ComponentDefinition yfp_cds = createCDS(doc, "YFP");
		createProtein(doc,"YFP",yfp_cds);

		ComponentDefinition pCONST_prom = createPromoter(doc, "pCONST");
		
		//creating promoter Sequences 
		String pCONST_prom_sequences = "GATAAGTCCCTAACTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC";
		Sequence pCONST_prom_seq= doc.createSequence("pCONST_sequence", version, pCONST_prom_sequences, Sequence.IUPAC_DNA);
		pCONST_prom.addSequence(pCONST_prom_seq);
		
		String pTac_prom_sequences = "AACGATCGTTGGCTGTGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGCTCACAATT";
		Sequence pTac_prom_seq= doc.createSequence("pTac_sequence", version, pTac_prom_sequences, Sequence.IUPAC_DNA);
		pTac_prom.addSequence(pTac_prom_seq);

		String pTet_prom_sequences = "TACTCCACCGTTGGCTTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCAC";
		Sequence pTet_prom_seq= doc.createSequence("pTet_sequence", version, pTet_prom_sequences, Sequence.IUPAC_DNA);
		pTet_prom.addSequence(pTet_prom_seq);

		String pBAD_prom_sequences = "ACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGGCTAGC";
		Sequence pBAD_prom_seq= doc.createSequence("pBAD_sequence", version, pBAD_prom_sequences, Sequence.IUPAC_DNA);
		pBAD_prom.addSequence(pBAD_prom_seq);

		String pLuxStar_prom_sequences = "ATAGCTTCTTACCGGACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTACTTTCGAATAAA";
		Sequence pLuxStar_prom_seq= doc.createSequence("pLuxStar_sequence", version, pLuxStar_prom_sequences, Sequence.IUPAC_DNA);
		pLuxStar_prom.addSequence(pLuxStar_prom_seq);
		
		//creating input sensor sequences
		String araC_cds_sequences = "atggctgaagcgcaaaatgatcccctgctgccgggatactcgtttaatgcccatctggtggcgggtttaacgccgattgaggccaacggttatctcgatttttttatcgaccgaccgctgggaatgaaaggttatattctcaatctcaccattcgcggtcagggggtggtgaaaaatcagggacgagaatttgtttgccgaccgggtgatattttgctgttcccgccaggagagattcatcactacggtcgtcatccggaggctcgcgaatggtatcaccagtgggtttactttcgtccgcgcgcctactggcatgaatggcttaactggccgtcaatatttgccaatacggggttctttcgcccggatgaagcgcaccagccgcatttcagcgacctgtttgggcaaatcattaacgccgggcaaggggaagggcgctattcggagctgctggcgataaatctgcttgagcaattgttactgcggcgcatggaagcgattaacgagtcgctccatccaccgatggataatcgggtacgcgaggcttgtcagtacatcagcgatcacctggcagacagcaattttgatatcgccagcgtcgcacagcatgtttgcttgtcgccgtcgcgtctgtcacatcttttccgccagcagttagggattagcgtcttaagctggcgcgaggaccaacgtatcagccaggcgaagctgcttttgagcaccacccggatgcctatcgccaccgtcggtcgcaatgttggttttgacgatcaactctatttctcgcgggtatttaaaaaatgcaccggggccagcccgagcgagttccgtgccggttgtgaagaaaaagtgaatgatgtagccgtcaagttgtcataa"; 
		Sequence araC_cds_seq = doc.createSequence("araC_sequence", version, araC_cds_sequences, Sequence.IUPAC_DNA);
		araC_cds.addSequence(araC_cds_seq);

		String tetR_cds_sequences = "atgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaatgaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgccttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctttactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggtttttcactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacctactactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggattagaaaaacaacttaaatgtgaaagtgggtcctaa";
		Sequence tetR_cds_seq = doc.createSequence("tetR_sequence", version, tetR_cds_sequences, Sequence.IUPAC_DNA);
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
		Sequence lacI_cds_seq = doc.createSequence("lacI_sequence", version, lacI_cds_sequences, Sequence.IUPAC_DNA);
		lacI_cds.addSequence(lacI_cds_seq);
		
		//creating output Sequences 
		String yfp_sequence = "CTGAAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAATACTAGAGAAAGAGGGGAAATACTAGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACAGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCAATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTTAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCC";
		Sequence yfp_seq= doc.createSequence("YFP_protein_sequence", version, yfp_sequence, Sequence.IUPAC_DNA);
		yfp_cds.addSequence(yfp_seq);

		String rfp_sequence = "CTGAAGTGGTCGTGATCTGAAACTCGATCACCTGATGAGCTCAAGGCAGAGCGAAACCACCTCTACAAATAATTTTGTTTAATACTAGAGTCACACAGGAAAGTACTAGATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAAGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAATAACAGATAAAAAAAATCCTTAGCTTTCGCTAAGGATGATTTCT";
		Sequence rfp_seq= doc.createSequence("RFP_protein_sequence", version, rfp_sequence, Sequence.IUPAC_DNA);
		rfp_cds.addSequence(rfp_seq);

		String bfp_sequence = "CTGAAGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence bfp_seq= doc.createSequence("BFP_protein_sequence", version, bfp_sequence, Sequence.IUPAC_DNA);
		bfp_cds.addSequence(bfp_seq);

		String gfp_sequence = "CTGAAGTTCCAGTCGAGACCTGAAGTGGGTTTCCTGATGAGGCTGTGGAGAGAGCGAAAGCTTTACTCCCGCACAAGCCGAAACTGGAACCTCTACAAATAATTTTGTTTAAGAGTCACACAGGAAAGTACTAGATGAGCGAGCTGATTAAGGAGAACATGCACATGAAGCTGTACATGGAGGGCACCGTGGACAACCATCACTTCAAGTGCACATCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACTAGCTTCCTCTACGGCAGCAAGACCTTCATCAACCACACCCAGGGCATCCCCGACTTCTTCAAGCAGTCCTTCCCTGAGGGCTTCACATGGGAGAGAGTCACCACATACGAAGATGGGGGCGTGCTGACCGCTACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCACATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCTTCACCGAGACGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAAACGACATGGCCCTGAAGCTCGTGGGCGGGAGCCATCTGATCGCAAACATCAAGACCACATATAGATCCAAGAAACCCGCTAAGAACCTCAAGATGCCTGGCGTCTACTATGTGGACTACAGACTGGAAAGAATCAAGGAGGCCAACAACGAGACCTACGTCGAGCAGCACGAGGTGGCAGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGGCACTAACCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA";
		Sequence gfp_seq= doc.createSequence("GFP_protein_sequence", version, gfp_sequence, Sequence.IUPAC_DNA);
		gfp_cds.addSequence(gfp_seq);
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
		ComponentDefinition cds = doc.createComponentDefinition(display, version, type);
		
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
				createProtein(document,name,componentDefinition);
			}
		}
	}
	
	private static void createProtein(SBOLDocument document,String cdsId,ComponentDefinition cds) throws SBOLValidationException 
	{
		ComponentDefinition proteinComponentDefinition =
				document.createComponentDefinition(cdsId+"_protein", version, ComponentDefinition.PROTEIN);
		proteinComponentDefinition.setName(cdsId+"_protein");
		proteinComponentDefinition.addWasDerivedFrom(derivedFrom);
		proteinComponentDefinition.addWasGeneratedBy(activityURI);
		proteinComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(cdsId+"_protein_production", version);
		moduleDefinition.setName(cdsId+"_protein_production");
		moduleDefinition.addWasDerivedFrom(derivedFrom);
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
		moduleDefinition.addWasDerivedFrom(derivedFrom);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(cdsId+"_protein", AccessType.PUBLIC, 
				proteinComponentDefinition.getIdentity(), DirectionType.NONE);
		String interactionId = cdsId + "_degradation_interaction";
		interaction = moduleDefinition.createInteraction(interactionId, SystemsBiologyOntology.DEGRADATION);
		interaction.createParticipation(cdsId+"_protein", cdsId+"_protein",  SystemsBiologyOntology.REACTANT);
	}
	
	private static void createInhibition(SBOLDocument document,String partId,String promoter) throws SBOLValidationException 
	{
		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(partId+"_"+promoter+"_repression", version);
		moduleDefinition.addWasDerivedFrom(derivedFrom);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(partId+"_protein", AccessType.PUBLIC, 
				partId+"_protein", version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(promoter, AccessType.PUBLIC, 
				promoter, version, DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(partId+"_"+promoter+"_repression", 
				SystemsBiologyOntology.INHIBITION);
		interaction.createParticipation(partId+"_protein_participation", partId+"_protein", 
				SystemsBiologyOntology.INHIBITOR);
		interaction.createParticipation(promoter+"_promoter_participation", promoter, 
				SystemsBiologyOntology.INHIBITED);
	}
	
	private static void createActivation(SBOLDocument document,String activator,String promoter) throws SBOLValidationException 
	{
		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(activator+"_"+promoter+"_repression", version);
		moduleDefinition.addWasDerivedFrom(derivedFrom);
		moduleDefinition.addWasGeneratedBy(activityURI);
		moduleDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
		moduleDefinition.createFunctionalComponent(activator, AccessType.PUBLIC, 
				activator, version, DirectionType.NONE);
		moduleDefinition.createFunctionalComponent(promoter, AccessType.PUBLIC, 
				promoter, version, DirectionType.NONE);
		Interaction interaction = moduleDefinition.createInteraction(activator+"_"+promoter+"_repression", 
				SystemsBiologyOntology.STIMULATION);
		interaction.createParticipation(activator+"_protein_participation", activator, 
				SystemsBiologyOntology.STIMULATOR);
		interaction.createParticipation(promoter+"_promoter_participation", promoter, 
				SystemsBiologyOntology.STIMULATED);
	}
	
	private static void createComplex(SBOLDocument document,String reactant1,String reactant2) throws SBOLValidationException 
	{
		String complex = reactant1 + "_" + reactant2;
		ComponentDefinition complexComponentDefinition = 
				document.createComponentDefinition(complex, version, ComponentDefinition.COMPLEX);

		ModuleDefinition moduleDefinition = 
				document.createModuleDefinition(complex+"_complex_formation", version);
		moduleDefinition.addWasDerivedFrom(derivedFrom);
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
		moduleDefinition.addWasDerivedFrom(derivedFrom);
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
							createInhibition(document,partId,promoter);
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
    				File file = new File(args[3] + gateName+"_gate_cytometry.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + args[2] + "/CelloParts/"+gateName+"/1"), 
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
