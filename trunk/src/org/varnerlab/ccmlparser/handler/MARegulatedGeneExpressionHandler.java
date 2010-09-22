package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MARegulatedGeneExpressionHandler extends CCMLMAObject implements IReactionHandler {

	private void init(Document ccmlTree) throws Exception
	{
		// Process keynames in the adapter_complex tags -
		String strACXPath = "//listOfGlobalSymbols/global_symbol/@key";
		NodeList nodeACList = (NodeList)_xpath.evaluate(strACXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_ACS= nodeACList.getLength();
		for (int index = 0;index<NUMBER_OF_ACS;index++)
		{
			// Get the key -
			Node tmpNode = nodeACList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = "//listOfGlobalSymbols/global_symbol[@key='"+strKeyName+"']/@symbol";
			String strSymbol = queryCCMLTree(ccmlTree,strSymbolXPath);
			
			// store in prop -
			setProperty(strKeyName,strSymbol);
		}
		
		// Get the list of prefixes -
		String strPrefixXPath = "//listOfSymbolPrefixes/symbol_prefix";
		populateProperties(strPrefixXPath,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize -
		init(ccmlTree);
		
		// Build the list -
		processRegulatedExpressionReactionBlock(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	private void processRegulatedExpressionReactionBlock(ArrayList<String> arrRxnList,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGenePrefixSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		String strXPath = "//regulated_gene/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_REGULATED_GENES = nodeList.getLength();
		for (int gene_index=0;gene_index<NUMBER_OF_REGULATED_GENES;gene_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(gene_index);
			String strRawGeneSymbol = tmpNode.getNodeValue();
			
			// Send the mRNA out of the nucleus -
			transportMRNAFromNucleus(arrRxnList,doc,strRawGeneSymbol);
			
			// Translation the mRNA -
			translateMRNA(arrRxnList,doc,strRawGeneSymbol);
			
			// Degrdation mRNA -
			degradationReactionMRNA(arrRxnList,doc,strRawGeneSymbol);
			
			// Degrade protein -
			degradationReactionProtein(arrRxnList,doc,strRawGeneSymbol);
			
			
			// Get the list of activators -
			String strActivatorsXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfActivators/activator/@symbol";
			ArrayList<String> listOfActivators = getRegulatorList(strActivatorsXPath,doc);
			int NUMBER_OF_ACTIVATORS = listOfActivators.size();
			for (int regulator_index=0;regulator_index<NUMBER_OF_ACTIVATORS;regulator_index++)
			{
				// Get the activator -
				String strActivator = listOfActivators.get(regulator_index);
				
				// Binding of TF to the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strActivator,strCompartment);
				
				// Binding of RNAP with the gene_TF complex -
				String strTFComplex = strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strActivator;
				generateReversiblebBindingReactions(arrRxnList,doc,strTFComplex,strRNAP,strCompartment);
				
				// Binding of the TF to the nuclear transporter (in) -
				arrReactants.add(strActivator+"_"+strTranslationCompartment);
				arrReactants.add(strImport+"_"+strCompartment);
				arrProducts.add(strActivator+"_"+strImport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction -
				arrReactants.add(strActivator+"_"+strImport+"_"+strCompartment);
				arrProducts.add(strImport+"_"+strCompartment);
				arrProducts.add(strActivator+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Binding of the TF to the nuclear transporter (out) -
				arrReactants.add(strActivator+"_"+strCompartment);
				arrReactants.add(strExport+"_"+strCompartment);
				arrProducts.add(strActivator+"_"+strExport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction (transport from the nucleus)
				arrReactants.add(strActivator+"_"+strExport+"_"+strCompartment);
				arrProducts.add(strExport+"_"+strCompartment);
				arrProducts.add(strActivator+"_"+strTranslationCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Generation of the MRNA -
				arrReactants.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strActivator+"_"+strRNAP+"_"+strCompartment);
				arrProducts.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				arrProducts.add(strActivator+"_"+strCompartment);
				arrProducts.add(strRNAP+"_"+strCompartment);
				arrProducts.add(strMRNAPrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			// Get the list of repressors -
			String strRepressorXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfRepressors/repressor/@symbol";
			ArrayList<String> listOfRepressors = getRegulatorList(strRepressorXPath,doc);
			int NUMBER_OF_REPRESSORS = listOfRepressors.size();
			for (int regulator_index=0;regulator_index<NUMBER_OF_REPRESSORS;regulator_index++)
			{
				// Get the activator -
				String strRepressor = listOfRepressors.get(regulator_index);
				
				// Binding of TF to the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strRepressor,strCompartment);
				
				// Ok, we need to get the repressors in and out of the nucleus -
				// Binding of the TF to the nuclear transporter (in) -
				arrReactants.add(strRepressor+"_"+strTranslationCompartment);
				arrReactants.add(strImport+"_"+strCompartment);
				arrProducts.add(strRepressor+"_"+strImport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport of TF into the nucleus 
				arrReactants.add(strRepressor+"_"+strImport+"_"+strCompartment);
				arrProducts.add(strImport+"_"+strCompartment);
				arrProducts.add(strRepressor+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Binding of the TF to the nuclear transporter (out) -
				arrReactants.add(strRepressor+"_"+strCompartment);
				arrReactants.add(strExport+"_"+strCompartment);
				arrProducts.add(strRepressor+"_"+strExport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport TF *out* of the nucleus
				arrReactants.add(strRepressor+"_"+strExport+"_"+strCompartment);
				arrProducts.add(strExport+"_"+strCompartment);
				arrProducts.add(strRepressor+"_"+strTranslationCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			// Process complex activator -
			processComplexActivator(arrRxnList,doc,strRawGeneSymbol);
			
			// Process complex deactivator -
			processComplexDeactivator(arrRxnList,doc,strRawGeneSymbol);
		}
	}
	
	private void processComplexDeactivator(ArrayList<String> arrRxnList,Document doc,String strRawGeneSymbol) throws Exception
	{
		
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGenePrefixSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		
		// Process the complex repressors -
		String strComplexRepressorXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexRepressors/complex_repressor/@symbol";
		ArrayList<String> listOfComplexRepressors = getRegulatorList(strComplexRepressorXPath,doc);
		int NUMBER_OF_COMPLEX_REPRESSORS = listOfComplexRepressors.size();
		for (int regulator_index=0;regulator_index<NUMBER_OF_COMPLEX_REPRESSORS;regulator_index++)
		{
			// Get symbol for this complex activator -
			String strComplexSymbol = listOfComplexRepressors.get(regulator_index);
			
			// Do we have any co-repressors?
			String strCoregulatorXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexRepressors/complex_repressor[@symbol='"+strComplexSymbol+"']/@corepressor";
			String strCorepressor = queryCCMLTree(doc,strCoregulatorXPath);
			
			if (strCorepressor.isEmpty())
			{
				// Bind the repressor with the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strComplexSymbol,strCompartment);
			}
			else
			{
				// Ok, crazy bitches ... I need to bind the co-repressor first and then make sure there is a transport reaction -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strCorepressor,strCompartment);
				
				// Transport -
				// Ok, I need to get the co-activator into the nucleus -
				// Binding of the TF to the nuclear transporter (in) -
				arrReactants.add(strCorepressor+"_"+strTranslationCompartment);
				arrReactants.add(strImport+"_"+strCompartment);
				arrProducts.add(strCorepressor+"_"+strImport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport of TF into the nucleus 
				arrReactants.add(strCorepressor+"_"+strImport+"_"+strCompartment);
				arrProducts.add(strImport+"_"+strCompartment);
				arrProducts.add(strCorepressor+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Binding of the TF to the nuclear transporter (out) -
				arrReactants.add(strCorepressor+"_"+strCompartment);
				arrReactants.add(strExport+"_"+strCompartment);
				arrProducts.add(strCorepressor+"_"+strExport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport TF *out* of the nucleus
				arrReactants.add(strCorepressor+"_"+strExport+"_"+strCompartment);
				arrProducts.add(strExport+"_"+strCompartment);
				arrProducts.add(strCorepressor+"_"+strTranslationCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			// Get where this complex assembles -
			String strAssemblyComplexCompartmentXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexRepressors/complex_repressor[@symbol='"+strComplexSymbol+"']/@compartment";
			String strAssemblyComplexCompartment = queryCCMLTree(doc,strAssemblyComplexCompartmentXPath);
		
			// Formulate the xpath to get the new the list of components = 
			String strComplexComponentsXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexRepressors/complex_repressor[@symbol='"+strComplexSymbol+"']/repressor/@symbol";
			ArrayList<String> listOfRepressorComponents = getRegulatorList(strComplexComponentsXPath,doc);
			
			int NUMBER_OF_COMPONENTS = listOfRepressorComponents.size();
			ArrayList<String> arrListReactants = new ArrayList<String>();
			ArrayList<String> arrListProducts = new ArrayList<String>();
			for (int complex_component=0;complex_component<NUMBER_OF_COMPONENTS;complex_component++)
			{
				// Get symbol for this complex activator -
				String strComplexComponent = listOfRepressorComponents.get(complex_component);
				arrListReactants.add(strComplexComponent+"_"+strAssemblyComplexCompartment);
			}
			
			// Ok, when I get here I have everything I need to assemble the complex activator -
			// Build the product -
			arrListProducts.add(strComplexSymbol+"_"+strAssemblyComplexCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrListReactants,arrListProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrListProducts,arrListReactants,ReactionType.REVERSE_RATE);
			arrListReactants.clear();
			arrListProducts.clear();
			
			// Transport the complex into the nucleus -
			
			// Ok, we need to get the repressors in and out of the nucleus -
			// Binding of the TF to the nuclear transporter (in) -
			arrReactants.add(strComplexSymbol+"_"+strTranslationCompartment);
			arrReactants.add(strImport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strImport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Finish the reaction - transport of TF into the nucleus 
			arrReactants.add(strComplexSymbol+"_"+strImport+"_"+strCompartment);
			arrProducts.add(strImport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Binding of the TF to the nuclear transporter (out) -
			arrReactants.add(strComplexSymbol+"_"+strCompartment);
			arrReactants.add(strExport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strExport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Finish the reaction - transport TF *out* of the nucleus
			arrReactants.add(strComplexSymbol+"_"+strExport+"_"+strCompartment);
			arrProducts.add(strExport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strTranslationCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
		}		
	}
	
	
	private void processComplexActivator(ArrayList<String> arrRxnList,Document doc,String strRawGeneSymbol) throws Exception
	{
	
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGenePrefixSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		
		// Process the list of complex activators -
		String strComplexActivatorXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexActivators/complex_activator/@symbol";
		ArrayList<String> listOfComplexActivators = getRegulatorList(strComplexActivatorXPath,doc);
		int NUMBER_OF_COMPLEX_ACTIVATORS = listOfComplexActivators.size();
		for (int regulator_index=0;regulator_index<NUMBER_OF_COMPLEX_ACTIVATORS;regulator_index++)
		{
			// Get symbol for this complex activator -
			String strComplexSymbol = listOfComplexActivators.get(regulator_index);
			
			// Do we have any co-activators?
			String strCoregulatorXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexActivators/complex_activator[@symbol='"+strComplexSymbol+"']/@coactivator";
			String strCoactivator = queryCCMLTree(doc,strCoregulatorXPath);
			
			// Get where this complex assembles -
			String strAssemblyComplexCompartmentXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexActivators/complex_activator[@symbol='"+strComplexSymbol+"']/@compartment";
			String strAssemblyComplexCompartment = queryCCMLTree(doc,strAssemblyComplexCompartmentXPath);
		
			// Formulate the xpath to get the new the list of components = 
			String strComplexComponentsXPath = "//regulated_gene[@symbol='"+strRawGeneSymbol+"']/listOfComplexActivators/complex_activator[@symbol='"+strComplexSymbol+"']/activator/@symbol";
			ArrayList<String> listOfActivatorComponents = getRegulatorList(strComplexComponentsXPath,doc);
			
			if (strCoactivator.isEmpty())
			{
				// Binding of TF to the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strComplexSymbol,strCompartment);
			
				// Binding of RNAP with the gene_TF complex -
				String strTFComplex = strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strComplexSymbol;
				generateReversiblebBindingReactions(arrRxnList,doc,strTFComplex,strRNAP,strCompartment);
				
				// Generation of mRNA -
				arrReactants.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strComplexSymbol+"_"+strRNAP+"_"+strCompartment);
				arrProducts.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				arrProducts.add(strComplexSymbol+"_"+strCompartment);
				arrProducts.add(strRNAP+"_"+strCompartment);
				arrProducts.add(strMRNAPrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			else
			{
				// If I get here then I have a co-activator that must bind first -
				// Binding of co-activator to the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol,strCoactivator,strCompartment);
			
				// Binding of the TF to the gene -
				generateReversiblebBindingReactions(arrRxnList,doc,strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCoactivator,strComplexSymbol,strCompartment);
				
				// Binding of RNAP with the gene_TF complex -
				String strTFComplex = strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCoactivator+"_"+strComplexSymbol;
				generateReversiblebBindingReactions(arrRxnList,doc,strTFComplex,strRNAP,strCompartment);
				
				// Generation of mRNA with coactivator -
				arrReactants.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCoactivator+"_"+strComplexSymbol+"_"+strRNAP+"_"+strCompartment);
				arrProducts.add(strGenePrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				arrProducts.add(strCoactivator+"_"+strCompartment);
				arrProducts.add(strComplexSymbol+"_"+strCompartment);
				arrProducts.add(strRNAP+"_"+strCompartment);
				arrProducts.add(strMRNAPrefixSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Ok, I need to get the co-activator into the nucleus -
				// Binding of the TF to the nuclear transporter (in) -
				arrReactants.add(strCoactivator+"_"+strTranslationCompartment);
				arrReactants.add(strImport+"_"+strCompartment);
				arrProducts.add(strCoactivator+"_"+strImport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport of TF into the nucleus 
				arrReactants.add(strCoactivator+"_"+strImport+"_"+strCompartment);
				arrProducts.add(strImport+"_"+strCompartment);
				arrProducts.add(strCoactivator+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Binding of the TF to the nuclear transporter (out) -
				arrReactants.add(strCoactivator+"_"+strCompartment);
				arrReactants.add(strExport+"_"+strCompartment);
				arrProducts.add(strCoactivator+"_"+strExport+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Finish the reaction - transport TF *out* of the nucleus
				arrReactants.add(strCoactivator+"_"+strExport+"_"+strCompartment);
				arrProducts.add(strExport+"_"+strCompartment);
				arrProducts.add(strCoactivator+"_"+strTranslationCompartment);
				encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
			}
			
			int NUMBER_OF_COMPONENTS = listOfActivatorComponents.size();
			ArrayList<String> arrListReactants = new ArrayList<String>();
			ArrayList<String> arrListProducts = new ArrayList<String>();
			for (int complex_component=0;complex_component<NUMBER_OF_COMPONENTS;complex_component++)
			{
				// Get symbol for this complex activator -
				String strComplexComponent = listOfActivatorComponents.get(complex_component);
				arrListReactants.add(strComplexComponent+"_"+strAssemblyComplexCompartment);
			}
			
			// Ok, when I get here I have everything I need to assemble the complex activator -
			// Build the product -
			arrListProducts.add(strComplexSymbol+"_"+strAssemblyComplexCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrListReactants,arrListProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrListProducts,arrListReactants,ReactionType.REVERSE_RATE);
			arrListReactants.clear();
			arrListProducts.clear();
			
			
			// Ok, we need to get the repressors in and out of the nucleus -
			// Binding of the TF to the nuclear transporter (in) -
			arrReactants.add(strComplexSymbol+"_"+strTranslationCompartment);
			arrReactants.add(strImport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strImport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Finish the reaction - transport of TF into the nucleus 
			arrReactants.add(strComplexSymbol+"_"+strImport+"_"+strCompartment);
			arrProducts.add(strImport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Binding of the TF to the nuclear transporter (out) -
			arrReactants.add(strComplexSymbol+"_"+strCompartment);
			arrReactants.add(strExport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strExport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Finish the reaction - transport TF *out* of the nucleus
			arrReactants.add(strComplexSymbol+"_"+strExport+"_"+strCompartment);
			arrProducts.add(strExport+"_"+strCompartment);
			arrProducts.add(strComplexSymbol+"_"+strTranslationCompartment);
			encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Binding of RNAP -
			
			
		}
	}
	
	private void degradationReactionProtein(ArrayList<String> arrRxnList,Document doc,String strMRNASymbol) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		
		// Encode the mRNA degrdation reaction -
		arrReactants.add(strMRNASymbol+"_"+strTranslationCompartment);
		arrProducts.add("[]");
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.DEGRADATION);
		arrReactants.clear();
		arrProducts.clear();	
	}
	
	private void degradationReactionMRNA(ArrayList<String> arrRxnList,Document doc,String strMRNASymbol) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGeneSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		// Encode the mRNA degrdation reaction -
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationCompartment);
		arrProducts.add("[]");
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.DEGRADATION);
		arrReactants.clear();
		arrProducts.clear();
	}
	
	private void translateMRNA(ArrayList<String> arrRxnList,Document doc,String strMRNASymbol) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGeneSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		// Ok, encode ribsome binding, scanning and generation of the protein -
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationCompartment);
		arrReactants.add(strTranslationRibosome+"_"+strTranslationCompartment);
		arrProducts.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Encode the scanning step =
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
		arrProducts.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strActiveSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Generate the protein -
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strActiveSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
		arrProducts.add(strTranslationRibosome+"_"+strTranslationCompartment);
		arrProducts.add(strMRNASymbol+"_"+strTranslationCompartment);
		arrProducts.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationCompartment);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}
	
	private void transportMRNAFromNucleus(ArrayList<String> arrRxnList,Document doc,String strMRNASymbol) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = (String)getProperty("TRANSCRIPTION_COMPARTMENT");
		String strTranslationCompartment = (String)getProperty("TRANSLATION_COMPARTMENT");
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNAPrefixSymbol = (String)getProperty("MRNA_PREFIX");
		String strGeneSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
		
		// Ok, encode this mofo - get on the transporter -
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strCompartment);
		arrReactants.add(strExport+"_"+strCompartment);
		arrProducts.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strExport+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Out of the nucleus -
		arrReactants.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strExport+"_"+strCompartment);
		arrProducts.add(strMRNAPrefixSymbol+"_"+strMRNASymbol+"_"+strTranslationCompartment);
		arrProducts.add(strExport+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);	
 	}
	
	private void generateReversiblebBindingReactions(ArrayList<String> arrRxnList,Document doc, String strGeneSymbol,String strOnSymbol,String strCompartment) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrListReactants = new ArrayList<String>();
		ArrayList<String> arrListProducts = new ArrayList<String>();
		
		// Ok build the reactants and products -
		String strComplex = strGeneSymbol+"_"+strOnSymbol+"_"+strCompartment;
		String strGeneSymbolWCompartment = strGeneSymbol+"_"+strCompartment;
		String strOnSymbolWCompartment = strOnSymbol+"_"+strCompartment;
		
		
		// Generate the ON-step -
		arrListReactants.add(strGeneSymbolWCompartment);
		arrListReactants.add(strOnSymbolWCompartment);
		arrListProducts.add(strComplex);
		
		// Call the reaction routine on-rate
		encodeMassActionSBMLReaction(arrRxnList,doc,arrListReactants,arrListProducts,ReactionType.FORWARD_RATE);
		
		// Call the reaction routine off-rate
		encodeMassActionSBMLReaction(arrRxnList,doc,arrListProducts,arrListReactants,ReactionType.REVERSE_RATE);
		
		// clear out the lists -
		arrListReactants.clear();
		arrListProducts.clear();
	}
}
