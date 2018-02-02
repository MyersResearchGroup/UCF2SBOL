package SBOLExamples.UCF2SBOL;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TimeZone;

import javax.xml.namespace.QName;

import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormatter;
import org.joda.time.format.ISODateTimeFormat;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.Collection;
import org.sbolstandard.core2.Component;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.DirectionType;
import org.sbolstandard.core2.GenericTopLevel;
import org.sbolstandard.core2.Interaction;
import org.sbolstandard.core2.ModuleDefinition;
import org.sbolstandard.core2.OrientationType;
import org.sbolstandard.core2.RestrictionType;
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
	static String celloNS = "https://github.com/CIDARLAB/cello/wiki/Cello-SBOL-Description#";

	static URI activityURI;
	static String createdDate;
	
	public static URI getRole(String type) {
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
	
	public static void convertPartsToSBOL(SBOLDocument document,HashMap<String,JSONObject> partsMap) throws SBOLValidationException {
		for (JSONObject part : partsMap.values()) {
			String name = (String)part.get("name");
			String dnasequence = (String)part.get("dnasequence");
			Sequence sequence = document.createSequence(name + "_sequence", version, dnasequence, Sequence.IUPAC_DNA);
			sequence.setName(name+"_sequence");
			sequence.addWasDerivedFrom(derivedFrom);
			sequence.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(name, version, ComponentDefinition.DNA);
			componentDefinition.setName(name);
			componentDefinition.addWasDerivedFrom(derivedFrom);
			componentDefinition.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
			componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			String partType = (String)part.get("type");
			componentDefinition.addRole(getRole(partType));
			componentDefinition.addSequence(sequence);
			
			if (partType.equals("cds")) {
				ComponentDefinition proteinComponentDefinition =
						document.createComponentDefinition(name+"_protein", version, ComponentDefinition.PROTEIN);
				proteinComponentDefinition.setName(name+"_protein");
				proteinComponentDefinition.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
				proteinComponentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

				ModuleDefinition moduleDefinition = 
						document.createModuleDefinition(name+"_protein_production", version);
				moduleDefinition.setName(name+"_protein_production");
				moduleDefinition.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
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
			HashMap<String,JSONObject> gatesMap) throws SBOLValidationException {
		for (JSONObject gate : gate_partsArr) {
			String gate_name = (String)gate.get("gate_name");
			ComponentDefinition componentDefinition = 
					document.createComponentDefinition(gate_name, version, ComponentDefinition.DNA);
			componentDefinition.setName(gate_name);
			componentDefinition.addRole(SequenceOntology.ENGINEERED_REGION);
			componentDefinition.addWasDerivedFrom(derivedFrom);
			componentDefinition.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
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
							moduleDefinition.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
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
			sequence.createAnnotation(new QName(provNS,"wasGeneratedBy","prov"), activityURI);
			sequence.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);
			componentDefinition.addSequence(sequence);
			
		}
	}
	
	public static void main( String[] args ) throws SBOLValidationException, SBOLConversionException, SynBioHubException, FileNotFoundException, IOException, ParseException
    {
		// Create an SBOLDocument
		SBOLDocument document = new SBOLDocument(); 
		document.setDefaultURIprefix(uriPrefix); 
		document.setComplete(true); 
		document.setCreateDefaults(true);
		
		// Create an Activity
		GenericTopLevel genericTopLevel = document.createGenericTopLevel("cello2sbol", version, 
				new QName(provNS, "Activity", "prov"));
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

		JSONParser parser = new JSONParser();
		JSONArray a = (JSONArray) parser.parse(new FileReader(
				"/Users/myers/git/CelloUCF2SynbioHub/ucf/Eco1C1G1T0.UCF.json"));

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
		}

		convertPartsToSBOL(document,partsMap);
        convertGatePartsToSBOL(document,gate_partsArr,gatesMap);
        
        Collection cdsCollection = document.createCollection("cdsCollection", version);
        for (ComponentDefinition cd : document.getComponentDefinitions()) {
        	if (cd.containsRole(SequenceOntology.CDS)) {
        		cdsCollection.addMember(cd.getIdentity());
        	}
        }
        
        // Validate
        SBOLValidate.validateSBOL(document,true,true,true);
        if (SBOLValidate.getNumErrors()>0) {
        	for (String error : SBOLValidate.getErrors()) {
        		System.out.println(error);
        	}
        } else {   
        	// Upload to SynBioHub
        	SynBioHubFrontend sbh = new SynBioHubFrontend("http://localhost:7777","https://synbiohub.org");
        	sbh.login("myers@ece.utah.edu", "test");
        	sbh.createCollection("CelloParts", "1", "Cello Parts", "These are the Cello parts", "27034378", true, document);
        	System.out.println("Conversion, validation, and upload successful");
        }

        //document.write(System.out);
    }
}
