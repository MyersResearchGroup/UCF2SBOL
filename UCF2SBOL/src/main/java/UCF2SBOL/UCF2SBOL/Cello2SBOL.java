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

/**
 * The Cello2SBOL class provides methods for converting biological parts and gates in CELLO UCF format
 * to SBOL (Synthetic Biology Open Language) format.
 * It handles components such as parts, gates, input sensors, output reporters,
 * as well as interactions like protein synthesis, inhibition, activation,
 * and complex formation.
 */
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
	
	private static URI getRole(String type) {
		if (type.equals("ribozyme")) {
	        return URI.create(so + "SO:0001977");
	    } else if (type.equals("scar")) {
	        return URI.create(so + "SO:0001953");
	    } else if (type.equals("cds")) {
	        return URI.create(so + "SO:0000316");
	    } else if (type.equals("promoter")) {
	        return URI.create(so + "SO:0000167");
	    } else if (type.equals("rbs")) {
	        return URI.create(so + "SO:0000139");
	    } else if (type.equals("cassette")) {
	    	return URI.create(so + "SO:0005853");
	    } else if (type.equals("terminator")) {
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

	/**
	 * Converts parts described in a HashMap to SBOL (Synthetic Biology Open Language) format.
	 *
	 * @param document The SBOLDocument to which the parts will be added.
	 * @param partsMap A HashMap containing parts represented as JSONObjects,
	 *                 where the key is a String identifier and the value is the corresponding JSONObject.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 * @throws IOException If there is an input/output error during the conversion process.
	 */
	private static void convertPartsToSBOL(SBOLDocument document,HashMap<String,JSONObject> partsMap) throws SBOLValidationException, IOException {

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
	
	/**
	 * Creates a protein component in the given SBOL document and its associated
	 * interactions and degradation pathways.
	 *
	 * @param document The SBOLDocument to which the protein component will be added.
	 * @param cdsId    The identifier for the coding sequence (CDS) that will be used to create the protein.
	 * @param cds      The ComponentDefinition for the CDS that will be linked to the protein production.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 */
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
	
	
	/**
	 * Creates an RNA component in the given SBOL document and its associated interactions
	 * and degradation pathways.
	 *
	 * @param document The SBOLDocument to which the RNA component will be added.
	 * @param rnaId    The identifier for the RNA that will be created.
	 * @param rna      The ComponentDefinition for the RNA that will be linked to the RNA production.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 */
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
	
	/**
	 * Creates an inhibition interaction in the given SBOL document between the inhibitor and the inhibited components,
	 * with optional parameters for ymin, ymax, alpha, beta, tau_on, and tau_off.
	 *
	 * @param document The SBOLDocument to which the inhibition interaction will be added.
	 * @param inhibitor The identifier for the inhibitor component.
	 * @param inhibited The identifier for the inhibited component.
	 * @param ymin Optional parameter representing the minimum output value.
	 * @param ymax Optional parameter representing the maximum output value.
	 * @param alpha Optional parameter representing the alpha value.
	 * @param beta Optional parameter representing the beta value.
	 * @param tau_on Optional parameter representing the tau_on value.
	 * @param tau_off Optional parameter representing the tau_off value.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 */
	private static void createInhibition(SBOLDocument document,String inhibitor,String inhibited,
			Double ymin,Double ymax,Double alpha,Double beta,Double tau_on,Double tau_off) throws SBOLValidationException 
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
		if (tau_on != null) {
			interaction.createAnnotation(new QName(celloNS,"tau_on","cello"),beta); 
		}
		if (tau_off != null) {
			interaction.createAnnotation(new QName(celloNS,"tau_off","cello"),beta); 
		}
	}
	
	/**
	 * Creates an activation interaction in the given SBOL document between the activator and the promoter components,
	 * with optional parameters for ymin, ymax, alpha, beta, tau_on, and tau_off.
	 *
	 * @param document The SBOLDocument to which the activation interaction will be added.
	 * @param activator The identifier for the activator component.
	 * @param promoter The identifier for the promoter component.
	 * @param ymin Optional parameter representing the minimum output value.
	 * @param ymax Optional parameter representing the maximum output value.
	 * @param alpha Optional parameter representing the alpha value.
	 * @param beta Optional parameter representing the beta value.
	 * @param tau_on Optional parameter representing the tau_on value.
	 * @param tau_off Optional parameter representing the tau_off value.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 */
	private static void createActivation(SBOLDocument document,String activator,String promoter,
			Double ymin,Double ymax,Double alpha,Double beta,Double tau_on,Double tau_off) throws SBOLValidationException 
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
		if (tau_on != null) {
			interaction.createAnnotation(new QName(celloNS,"tau_on","cello"),tau_on); 
		}
		if (tau_off != null) {
			interaction.createAnnotation(new QName(celloNS,"tau_off","cello"),tau_off); 
		}
	}
	
	/**
	 * Creates a complex component from two reactant components within the given SBOL document,
	 * including the necessary interactions and degradation pathways.
	 *
	 * @param document The SBOLDocument to which the complex component will be added.
	 * @param reactant1 The identifier for the first reactant component.
	 * @param reactant2 The identifier for the second reactant component.
	 * @throws SBOLValidationException If there is an error validating the SBOL document.
	 */
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

	/**
	 * Converts the gate parts to SBOL format and adds them to the given SBOL document.
	 *
	 * @param document      the SBOL document to which gate parts will be added
	 * @param gate_partsArr a set containing JSONObjects representing the gate parts
	 * @param gatesMap      a map containing gate information with the gate name as the key
	 * @param responseMap   a map containing response function information with the gate name as the key
	 * @param functionMap   a map containing function information used to generate response functions
	 * @throws SBOLValidationException if there is an error in the SBOL validation process
	 */
	private static void convertGatePartsToSBOL(SBOLDocument document,HashSet<JSONObject> gate_partsArr,
			HashMap<String,JSONObject> gatesMap,HashMap<String,JSONObject> responseMap, HashMap<String,JSONObject> functionMap) throws SBOLValidationException {
		for (JSONObject gate : gate_partsArr) {
			boolean v2 = (functionMap != null);
			String gate_name = (String)gate.get("gate_name");
			String respfxn = null;
			
			if(v2) {
				gate_name = (String)gate.get("name");
				gate_name = gate_name.substring(0, gate_name.length()-10);
				respfxn = (String) ((functionMap.get((String)(((JSONObject)responseMap.get(gate_name).get("functions")).get("response_function")))).get("equation"));
			}
			
			else {
				respfxn = (String)responseMap.get(gate_name).get("equation");
			}
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
	        		(String)gatesMap.get(gate_name).get(v2 ? "group" : "group_name"));
	        componentDefinition.createAnnotation(new QName(celloNS,"color_hexcode","cello"), 
	        		(String)gatesMap.get(gate_name).get(v2 ? "color" : "color_hexcode"));
	        componentDefinition.createAnnotation(new QName(celloNS,"response_function","cello"), 
	        		respfxn);
	        if (((JSONObject)responseMap.get(gate_name).get("functions")).get("tandem_interference_factor") != null) {
	        	componentDefinition.createAnnotation(new QName(celloNS,"tandem_efficiency_factor","cello"), 
	        			(String) ((functionMap.get((String)(((JSONObject)responseMap.get(gate_name).get("functions")).get("tandem_interference_factor")))).get("equation")));
	        }
	        if (((JSONObject)responseMap.get(gate_name).get("functions")).get("tandem_input_composition") != null) {
	        	componentDefinition.createAnnotation(new QName(celloNS,"tandem_input_composition","cello"), 
	        			(String) ((functionMap.get((String)(((JSONObject)responseMap.get(gate_name).get("functions")).get("tandem_interference_factor")))).get("equation")));
	        }
	        
	        JSONArray parameters = (JSONArray)responseMap.get(gate_name).get("parameters");
	        for (Object obj : parameters) {
	        	String name = (String)((JSONObject)obj).get("name");
	        	componentDefinition.createAnnotation(new QName(celloNS,name,"cello"), 
	        			(Double)((JSONObject)obj).get("value"));
	        }
//	        JSONArray variables = (JSONArray)responseMap.get(gate_name).get("variables");
//	        for (Object obj : variables) {
//	        	String name = (String)((JSONObject)obj).get("name");
//	        	componentDefinition.createAnnotation(new QName(celloNS,name+"_off_threshold","cello"), 
//	        			(Double)((JSONObject)obj).get("off_threshold"));
//	        	componentDefinition.createAnnotation(new QName(celloNS,name+"_on_threshold","cello"), 
//	        			(Double)((JSONObject)obj).get("on_threshold"));
//	        }
	        JSONArray expression_cassettes = v2 ? (JSONArray) gate.get("devices") : (JSONArray) gate.get("expression_cassettes");
			String seq = "";
			for (Object obj : expression_cassettes) {
				int annotationCount = 0;
				int start = 1;
				//int constraintCount = 0;
				//Component previousComponent = null;
//		        Component currentComponent = null;
		        
				JSONObject expression_cassette = (JSONObject) obj;
				JSONArray cassette_parts = v2 ? (JSONArray)expression_cassette.get("components") : (JSONArray)expression_cassette.get("cassette_parts");
				
				if (((String)cassette_parts.get(0)).startsWith("#in")) { // && ((String)cassette_parts.get(1)).startsWith("#in")) {
					continue;
				}
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
						String promoter = v2 ? (String)((JSONArray)gate.get("outputs")).get(0) : (String)gate.get("promoter");
						if (document.getModuleDefinition(partId+"_protein_"+promoter+"_repression", version)==null) {
							createInhibition(document,partId+"_protein",promoter,null,null,null,null,null,null);
						}
					}
					if (partComponentDefinition.getRoles().contains(URI.create(so + "SO:0001264"))) {
						String promoter = v2 ? (String)((JSONArray)gate.get("outputs")).get(0) : (String)gate.get("promoter");
						createComplex(document,partId+"_rna","dCAS9_Mxi1_protein");
						if (document.getModuleDefinition(partId+"_rna_dCAS9_Mxi1_protein_"+promoter+"_repression", version)==null) {
							createInhibition(document,partId+"_rna_dCAS9_Mxi1_protein",promoter,null,null,null,null,null,null);
						}
					}
				}
				break;
				
			}
			
			Sequence sequence = document.createSequence(gate_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(gate_name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}

	/**
	 * Converts input sensor data to SBOL (Synthetic Biology Open Language) format and adds it to the SBOLDocument.
	 *
	 * @param document The SBOLDocument to which the converted input sensors will be added.
	 * @param input_sensorsArr A HashSet of JSONObjects, each representing an input sensor.
	 * @param responseMap A HashMap where the key is the sensor name and the value is a JSONObject representing the response for each sensor.
	 * @param functionMap A HashMap where the key is the function name and the value is a JSONObject representing the function details.
	 * @throws SBOLValidationException If there is a validation error when creating or modifying SBOL objects.
	 */
	private static void convertInputSensorsToSBOL(SBOLDocument document,HashSet<JSONObject> input_sensorsArr,HashMap<String,JSONObject> responseMap, HashMap<String,JSONObject> functionMap) throws SBOLValidationException {
		for (JSONObject sensor : input_sensorsArr) {
			String sensor_name = (String)sensor.get("name");
			boolean v2 = (responseMap != null);
			String respfxn = null;
			
			if(v2) {
				sensor_name = sensor_name.substring(0, sensor_name.length()-10);
				respfxn = (String) ((functionMap.get((String)(((JSONObject)responseMap.get(sensor_name).get("functions")).get("response_function")))).get("equation"));
			}
			
			else {
				respfxn = (String)responseMap.get(sensor_name).get("equation");
			}
			
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(sensor_name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(sensor_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
	        componentDefinition.createAnnotation(new QName(celloNS,"gateType","cello"), "input_sensor");
	        componentDefinition.createAnnotation(new QName(celloNS,"response_function","cello"), 
	        		respfxn);
	        
	        if (((JSONObject)responseMap.get(sensor_name).get("functions")).get("tandem_interference_factor") != null) {
	        	componentDefinition.createAnnotation(new QName(celloNS,"tandem_efficiency_factor","cello"), 
	        			(String) ((functionMap.get((String)(((JSONObject)responseMap.get(sensor_name).get("functions")).get("tandem_interference_factor")))).get("equation")));
	        }
	        
	        JSONArray parameters = (JSONArray)responseMap.get(sensor_name).get("parameters");
	        if (parameters != null) {
	        	for (Object obj : parameters) {
		        	String name = (String)((JSONObject)obj).get("name");
		        	componentDefinition.createAnnotation(new QName(celloNS,name,"cello"), 
		        			(Double)((JSONObject)obj).get("value"));
		        }
	        }
			JSONArray parts = (JSONArray)sensor.get(v2 ? "outputs" : "parts");
			JSONObject oldsensor = sensor;
			sensor = (v2 ? responseMap.get(sensor_name) : sensor);
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
					Double tau_on = null;
					Double tau_off = null;

			        if (parameters != null) {
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
			        }

					document.createComponentDefinition(input_molecule, version, ComponentDefinition.SMALL_MOLECULE);
					createComplex(document,input_molecule,partId+"_protein");
					if (((String)sensor.get("type")).equals("complex_stimulator")) {
						createActivation(document,input_molecule+"_"+partId+"_protein",promoter,signal_low,signal_high,alpha,beta,tau_on,tau_off);
					} else if (((String)sensor.get("type")).equals("sequester_inhibitor")) {
						createInhibition(document,partId+"_protein",promoter,signal_low,signal_high,alpha,beta,tau_on,tau_off);
					}
				}
				
			}
			sensor = oldsensor;
			Sequence sequence = document.createSequence(sensor_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(sensor_name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}

	/**
	 * Converts output reporters to SBOL (Synthetic Biology Open Language) format components
	 * and adds them to the provided SBOL document.
	 *
	 * @param document The SBOLDocument object where the output reporters will be added.
	 * @param output_reportersArr A HashSet containing JSONObjects representing the output reporters.
	 * @param responseMap A HashMap mapping reporter names to JSONObjects that include the response functions and parameters.
	 * @param functionMap A HashMap mapping function names to JSONObjects that include the equations and other function details.
	 * @throws SBOLValidationException If there is an error in the SBOL validation process.
	 */
	private static void convertOutputReportersToSBOL(SBOLDocument document,HashSet<JSONObject> output_reportersArr,HashMap<String,JSONObject> responseMap, HashMap<String,JSONObject> functionMap) throws SBOLValidationException {
		for (JSONObject sensor : output_reportersArr) {
			String reporter_name = (String)sensor.get("name");
			boolean v2 = (responseMap != null);
			String respfxn = null;
			
			if(v2) {
				reporter_name = reporter_name.substring(0, reporter_name.length()-10);
				respfxn = (String) ((functionMap.get((String)(((JSONObject)responseMap.get(reporter_name).get("functions")).get("response_function")))).get("equation"));
			}
			
			else {
				respfxn = (String)responseMap.get(reporter_name).get("equation");
			}
			
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(reporter_name, version, ComponentDefinition.DNA_REGION);
			componentDefinition.setName(reporter_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasGeneratedBy(activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
	        componentDefinition.createAnnotation(new QName(celloNS,"gateType","cello"), "output_reporter");
	        componentDefinition.createAnnotation(new QName(celloNS,"response_function","cello"), 
	        		respfxn);
	        
	        if (((JSONObject)responseMap.get(reporter_name).get("functions")).get("tandem_interference_factor") != null) {
	        	componentDefinition.createAnnotation(new QName(celloNS,"tandem_efficiency_factor","cello"), 
	        			(String) ((functionMap.get((String)(((JSONObject)responseMap.get(reporter_name).get("functions")).get("tandem_interference_factor")))).get("equation")));
	        }
	        
	        JSONArray parameters = (JSONArray)responseMap.get(reporter_name).get("parameters");
	        if (parameters != null) {
	        	for (Object obj : parameters) {
		        	String name = (String)((JSONObject)obj).get("name");
		        	componentDefinition.createAnnotation(new QName(celloNS,name,"cello"), 
		        			(Double)((JSONObject)obj).get("value"));
		        }
	        }		        
			JSONArray parts = v2 ? (JSONArray) (((JSONObject) ((JSONArray)sensor.get("devices")).get(0)).get("components")): (JSONArray)sensor.get("parts");
			
			String seq = "";
			int annotationCount = 0;
			int start = 1;
			for (Object obj2 : parts) {
				if (((String)obj2).startsWith("#in")) {
					continue;
				}
				String partId = (String)obj2;
				ComponentDefinition currentComponent = document.getComponentDefinition(partId, version);
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
				
				String output_protein = partId.replace("_cassette", "");
				createProtein(document,output_protein,currentComponent);
			}
			
			Sequence sequence = document.createSequence(reporter_name+"_sequence", version, seq, Sequence.IUPAC_DNA);
			sequence.setName(reporter_name+"_sequence");
			sequence.addWasGeneratedBy(activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}
	
	/**
	 * Checks if a given JSON array contains both "models" and "structures" collections.
	 *
	 * @param a JSONArray to be checked.
	 * @return true if both "models" and "structures" are present, false otherwise.
	 */
	public static boolean isV2(JSONArray a) {
		boolean models = false;
		boolean structures = false;
		for (Object o : a)
		{
			JSONObject ucf = (JSONObject) o;
			
			String collection = (String) ucf.get("collection");
			if(collection.equals("models")) {
				models = true;
			}
			else if(collection.equals("structures")) {
				structures = true;
			}
			
			if(models && structures) {
				return true;
			}
			

		}
		return false;
	}
	
	/**
	 * The main method to convert UCF to SBOL
	 *
	 * @param args the array of string arguments where:
	 *        args[0] - login email
	 *        args[1] - password
	 *        args[2] - login user
	 *        args[3] - temporary directory
	 *        args[4] - database prefix
	 *        args[5] - path to UCF file
	 *        args[6] - path to input file (optional)
	 *        args[7] - path to output file (optional)
	 *        args[8] - database URL (optional)
	 *        args[9] - collection id (optional)
	 *        args[10] - collection version (optional)
	 *        args[11] - collection name (optional)
	 *        args[12] - collection description (optional)
	 *        args[13] - collection pubMedId (optional)
	 * @throws SBOLValidationException if an SBOL validation error occurs
	 * @throws SBOLConversionException if an SBOL conversion error occurs
	 * @throws SynBioHubException if a SynBioHub related error occurs
	 * @throws FileNotFoundException if a specified file is not found
	 * @throws IOException if an I/O error occurs
	 * @throws ParseException if a parsing error occurs
	 * @throws URISyntaxException if a URI syntax error occurs
	 */
	// args[0] - login email
	// args[1] - password
	// args[2] - login user
	// args[3] - temporary directory
	// args[4] - databasePrefix
	// args[5] - path to UCF file
	// args[6] - path to input file
	// args[7] - path to output file
	// args[8] - databaseURL
	// args[9] - collection id
	// args[10] - collection version
	// args[11] - collection name
	// args[12] - collection description
	// args[13] - collection pubMedId
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
			System.err.println(" path to UCF input file");
			System.err.println(" path to UCF outpuf file");
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
		String pathToInputFile = null;
		String pathToOutputFile = null;				
		String databaseURL = databasePrefix;
		String collectionId = "Eco1C1G1T1_Parts";
		String collectionVersion = "1";
		String collectionName = "Cello Parts";
		String collectionDescription = "These are the Cello parts";
		String collectionPubMedId = "27034378";
		if (args.length > 6) {
			pathToInputFile = args[6];
		}
		if (args.length > 7) {
			pathToOutputFile = args[7];
		}
		if (args.length > 8) {
			databaseURL = args[8];
		}
		if (args.length > 9) {
			collectionId = args[9];
		}
		if (args.length > 10) {
			collectionVersion = args[10];
		}
		if (args.length > 11) {
			collectionName = args[11];
		}
		if (args.length > 12) {
			collectionDescription = args[12];
		}
		if (args.length > 13) {
			collectionPubMedId = args[13];
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
		HashMap<String,JSONObject> functionMap = new HashMap<String,JSONObject>();

		JSONParser parser = new JSONParser();
		JSONArray a = (JSONArray) parser.parse(new FileReader(pathToUCFFile));
		JSONArray in = null;
		JSONArray out = null;
		boolean v2 = isV2(a);
		if (v2) {
			in = (JSONArray) parser.parse(new FileReader(pathToInputFile));
			out = (JSONArray) parser.parse(new FileReader(pathToOutputFile));
			for (Object o : in) {
				JSONObject ucf = (JSONObject) o;
				String collection = (String) ucf.get("collection");
				if (collection.equals("parts")) {
					partsMap.put((String)ucf.get("name"),ucf);
				} else if (collection.equals("functions")) {
					functionMap.put((String)ucf.get("name"),ucf);
				} else if (collection.equals("structures")) {
					input_sensorsArr.add(ucf);
				} else if (collection.equals("models")) {
					String name = (String)ucf.get("name");
					responseMap.put(name.substring(0, name.length()-6),ucf);
				}
			}
			for (Object o : out)
			{
				JSONObject ucf = (JSONObject) o;

				String collection = (String) ucf.get("collection");
				if (collection.equals("functions")) {
					functionMap.put((String)ucf.get("name"),ucf);
				} else if (collection.equals("parts")) {
					partsMap.put((String)ucf.get("name"),ucf);
				} else if (collection.equals("structures")) {
					output_reportersArr.add(ucf);
				} else if (collection.equals("models")) {
					String name = (String)ucf.get("name");
					responseMap.put(name.substring(0, name.length()-6),ucf);
				}
			}

			for (Object o : a)
			{
				JSONObject ucf = (JSONObject) o;

				String collection = (String) ucf.get("collection");
				if (collection.equals("functions")) {
					functionMap.put((String)ucf.get("name"),ucf);
				}
			}
			for (Object o : a)
			{
				JSONObject ucf = (JSONObject) o;

				String collection = (String) ucf.get("collection");

				if (collection.equals("parts")) {
					partsMap.put((String)ucf.get("name"),ucf);
				}
				else if (collection.equals("structures")) {
					gate_partsArr.add(ucf);
				}
				else if (collection.equals("gates")) {
					gatesMap.put((String)ucf.get("name"),ucf);
				}
				else if (collection.equals("models")) {
					String name = (String)ucf.get("name");
					responseMap.put(name.substring(0, name.length()-6),ucf);
				}
			}
		}
		else {
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
		if (v2) {
	        convertGatePartsToSBOL(document,gate_partsArr,gatesMap,responseMap,functionMap);
	        convertInputSensorsToSBOL(document,input_sensorsArr,responseMap,functionMap);
	        convertOutputReportersToSBOL(document,output_reportersArr,responseMap,functionMap);
		}
		else {
	        convertGatePartsToSBOL(document,gate_partsArr,gatesMap,responseMap,null);
	        convertInputSensorsToSBOL(document,input_sensorsArr,null,null);
	        convertOutputReportersToSBOL(document,output_reportersArr,null,null);
		}

        
        //createSensorsReporters(document);
        
//        Collection cdsCollection = document.createCollection("cdsCollection", version);
//        for (ComponentDefinition cd : document.getComponentDefinitions()) {
//        	if (cd.containsRole(SequenceOntology.CDS)) {
//        		cdsCollection.addMember(cd.getIdentity());
//        	}
//        }
        
        // Validate
        SBOLValidate.validateSBOL(document,true,true,true);
//        document.write(collectionId + ".SBOL");
        int i = 0;
        if (SBOLValidate.getNumErrors()>0) {
        	for (String error : SBOLValidate.getErrors()) {
        		System.out.println(error);
        	}
        } else {   
        	// Upload to SynBioHub
        	SynBioHubFrontend sbh = new SynBioHubFrontend(databaseURL,databasePrefix);
        	sbh.login(loginEmail, password);
        	sbh.createCollection(collectionId, collectionVersion, collectionName, collectionDescription,
        			collectionPubMedId, true);
        	sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), pathToUCFFile);

			// Define the base URI
			String baseURI = databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion;

			// Print the URL for attaching the main file
			System.out.println("Attaching file to URL: " + baseURI);



			if (v2) {
        		sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), pathToInputFile);
            	sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), pathToOutputFile);

        	}
        	SBOLDocument doc = sbh.getSBOL(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion));
        	for (Attachment attachment : doc.getAttachments()) {
        		switch(i) {
        		case 0: activity.createUsage("UCF_file", attachment.getIdentity());
        		i += 1;
        		
        		case 1: activity.createUsage("input_file", attachment.getIdentity());
        		i += 1;
        		
        		case 2: activity.createUsage("output_file", attachment.getIdentity());
        		i += 1;
        		}
        		
        		
        	}
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
    			} else if (collection.equals("circuit_rules")) {
    				File file = new File(args[3] + "circuit_rules.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "circuit_rules.json");
    			} else if (collection.equals("functions")) {
    				String gateName = (String)ucf.get("name");
    				if (gateName.endsWith("cytometry")) {
        				File file = new File(args[3] + gateName+".json");
        				FileOutputStream stream = new FileOutputStream(file);
        				BufferedOutputStream buffer = new BufferedOutputStream(stream);
        				stream.write(ucf.toJSONString().getBytes());
        				stream.close();
        				buffer.close();
        				String SBHgateName = gateName.replace("_cytometry", "");
        				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + SBHgateName + "/" + collectionVersion), 
        						tmpDir + gateName+".json");
    				} else if (gateName.endsWith("toxicity")) {
        				File file = new File(args[3] + gateName+".json");
        				FileOutputStream stream = new FileOutputStream(file);
        				BufferedOutputStream buffer = new BufferedOutputStream(stream);
        				stream.write(ucf.toJSONString().getBytes());
        				stream.close();
        				buffer.close();
        				String SBHgateName = gateName.replace("_toxicity", "");
        				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + SBHgateName + "/" + collectionVersion), 
        						tmpDir + gateName+".json");
    				}
    			} else if (collection.equals("device_rules")) {
    				File file = new File(args[3] + "device_rules.json");
    				FileOutputStream stream = new FileOutputStream(file);
    				BufferedOutputStream buffer = new BufferedOutputStream(stream);
    				stream.write(ucf.toJSONString().getBytes());
    				stream.close();
    				buffer.close();
    				sbh.attachFile(URI.create(databasePrefix + "/user/" + loginUser + "/" + collectionId + "/" + collectionId + "_collection/" + collectionVersion), 
    						tmpDir + "device_rules.json");
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
    			} else if (collection.equals("functions")) {
    			} else if (collection.equals("response_functions")) {
    			} else if (collection.equals("models")) {
    			} else if (collection.equals("structures")) {
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
