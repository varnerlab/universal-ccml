package org.varnerlab.ccmlparser;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public abstract class CCMLMAObject {
	protected XPathFactory  _xpFactory = XPathFactory.newInstance();
	protected XPath _xpath = _xpFactory.newXPath();
	protected ArrayList<String> _arrListSBMLReactions = new ArrayList<String>();
	protected ArrayList<ReactionType> _arrListReactionType = new ArrayList<ReactionType>();
	protected Hashtable<String,Object> _propTable = new Hashtable<String,Object>();
	
	public Object getProperty(String key) {
		return(_propTable.get(key));
	}

	public void setProperty(String key, Object value) {
		_propTable.put(key, value);
	}
	
	protected void buildInterfaceReactions(String strBlockName,ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		// Get prefix -
		String strPhosphorylationPrefix = (String)getProperty("PHOSPHORYLATION_PREFIX");
		String strDoublePhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strActivatedPrefix = (String)getProperty("ACTIVATED_PREFIX");
		String strDeativatedPrefix = (String)getProperty("DEACTIVATED_PREFIX");
		String strUBPrefix = (String)getProperty("UBIQUITIN_PREFIX");
		
		
		// Ok, let's get the list interface blocks -
		String strBlockXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfInterfaces/interface/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strBlockXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_INTERFACES= nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_INTERFACES;index++)
		{
			// Get the interface symbol -
			Node tmpNode = nodeList.item(index);
			String strInterfaceSymbol = tmpNode.getNodeValue();
			
			// Ok now that I have the interface symbol - process the targets (activate) 
			String strTargetXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_phosphorylate/@symbol";
			NodeList nodeTargetList = (NodeList)_xpath.evaluate(strTargetXPath,ccmlTree,XPathConstants.NODESET);
			int NUMBER_OF_TARGETS= nodeTargetList.getLength();
			for (int target_index=0;target_index<NUMBER_OF_TARGETS;target_index++)
			{
				// Get the target symbol -
				Node tmpTargetNode = nodeTargetList.item(target_index);
				String strTargetSymbol = tmpTargetNode.getNodeValue();
				
				// Process the targets -
				// Binding -
				arrReactants.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrReactants.add(strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Phosphorylate that bitch ..
				arrReactants.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strPhosphorylationPrefix+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			// Ok, process the Ub interface -
			// Ok now that I have the interface symbol - process the targets (activate) 
			String strTargetActivateUbXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_ubiquitylation/@symbol";
			NodeList nodeTargetActivateUbList = (NodeList)_xpath.evaluate(strTargetActivateUbXPath,ccmlTree,XPathConstants.NODESET);
			NUMBER_OF_TARGETS= nodeTargetActivateUbList.getLength();
			for (int target_index=0;target_index<NUMBER_OF_TARGETS;target_index++)
			{
				// Get the target symbol -
				Node tmpTargetNode = nodeTargetActivateUbList.item(target_index);
				String strTargetSymbol = tmpTargetNode.getNodeValue();
				
				// Process the targets -
				// Binding -
				arrReactants.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrReactants.add(strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Phosphorylate that bitch ..
				arrReactants.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strUBPrefix+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			// Ok now that I have the interface symbol - process the targets (activate) 
			String strTargetActivateXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_activate/@symbol";
			NodeList nodeTargetActivateList = (NodeList)_xpath.evaluate(strTargetActivateXPath,ccmlTree,XPathConstants.NODESET);
			NUMBER_OF_TARGETS= nodeTargetActivateList.getLength();
			for (int target_index=0;target_index<NUMBER_OF_TARGETS;target_index++)
			{
				// Get the target symbol -
				Node tmpTargetNode = nodeTargetActivateList.item(target_index);
				String strTargetSymbol = tmpTargetNode.getNodeValue();
				
				// Process the targets -
				// Binding -
				arrReactants.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrReactants.add(strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Phosphorylate that bitch ..
				arrReactants.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strActivatedPrefix+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
			
			String strTargetDeactivateXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_dephosphorylate/@symbol";
			NodeList nodeTargetDeactivateList = (NodeList)_xpath.evaluate(strTargetDeactivateXPath,ccmlTree,XPathConstants.NODESET);
			int NUMBER_OF_DEACTIVATE_TARGETS= nodeTargetDeactivateList.getLength();
			for (int target_index=0;target_index<NUMBER_OF_DEACTIVATE_TARGETS;target_index++)
			{
				// Get the target symbol -
				Node tmpTargetNode = nodeTargetDeactivateList.item(target_index);
				String strTargetSymbol = tmpTargetNode.getNodeValue();
				
				// Find the *last* _
				int INDEX_LAST_UNDERSCORE = strTargetSymbol.indexOf("_");
				int LENGTH = strTargetSymbol.length();
				
				// Get the initiator -
				String strDeactivatedTarget = strTargetSymbol.substring(INDEX_LAST_UNDERSCORE+1,LENGTH);
				
				// Process the targets -
				// Binding -
				arrReactants.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrReactants.add(strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// De-Phosphorylate that bitch ..
				arrReactants.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strReceptorCompartment);
				arrProducts.add(strDeactivatedTarget+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
		}
		
	}
	
	protected void init(String strBlockName, Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
		
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
		
		// List of prefixes -
		String strPrefixXPath = "//listOfSymbolPrefixes/symbol_prefix";
		populateProperties(strPrefixXPath,ccmlTree);
		
		// Load the siganling components -
		String strXPathSignalingComponents = "//signaling_block[@block_class='"+strBlockName+"']/listOfSignalingComponents/signaling_component";
		populateProperties(strXPathSignalingComponents,ccmlTree);
		
		// Get the initiator symbols -
		String strXPathSignalingIniators = "//signaling_block[@block_class='"+strBlockName+"']/listOfInitiators/initiator";
		populateProperties(strXPathSignalingIniators,ccmlTree);
	}
	
	// Logic to build the degradation -
	protected void buildDegradationReactions(String strBlockName,ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//signaling_block[@block_class='"+strBlockName+"']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//signaling_block[@block_class='"+strBlockName+"']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(ccmlTree,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	
	// Code to load key-symbol pairs
	protected void populateProperties(String strXPathBase,Document ccmlTree) throws Exception
	{
		String strSignalingXPathKey = strXPathBase+"/@key";
		NodeList nodeACList = (NodeList)_xpath.evaluate(strSignalingXPathKey,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_ACS= nodeACList.getLength();
		for (int index = 0;index<NUMBER_OF_ACS;index++)
		{
			// Get the key -
			Node tmpNode = nodeACList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = strXPathBase+"[@key='"+strKeyName+"']/@symbol";
			String strSymbol = queryCCMLTree(ccmlTree,strSymbolXPath);
			
			// store in prop -
			setProperty(strKeyName,strSymbol);
		}
	}
	
	// Get a string -
	protected String queryCCMLTree(Document ccmlTree,String strXPath)
	{
		// Method attributes -
		String strProp = "";
		
		try {
			Node propNode = (Node) _xpath.evaluate(strXPath, ccmlTree, XPathConstants.NODE);
			strProp = propNode.getNodeValue();
		}
		catch (Exception error)
		{
			error.printStackTrace();
			System.out.println("ERROR: Property lookup failed on CCMLTree. The following XPath "+strXPath+" resuled in an error - "+error.toString());
		}
		
		return(strProp);
	}
	
	protected ArrayList<String> getRegulatorList(String strXPath,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> tmpList = new ArrayList<String>();
		
		NodeList nodeListRegulators = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_REGULATORS = nodeListRegulators.getLength();
		for (int regulator_index=0;regulator_index<NUMBER_OF_REGULATORS;regulator_index++)
		{
			// Get the activator symbols -
			Node tmpNodeRegulator = nodeListRegulators.item(regulator_index);
			String strRawGeneSymbol = tmpNodeRegulator.getNodeValue();
			
			// Add the regulator to the list -
			tmpList.add(strRawGeneSymbol);
		}
		
		// return -
		return(tmpList);
	}
	
	protected void encodeMassActionSBMLReaction(ArrayList<String> buffer,Document doc,ArrayList<String> arrReactants,ArrayList<String> arrProducts,ReactionType rxnType) throws Exception
	{
		// Method attributes -
		StringBuffer local_buffer = new StringBuffer();
		
		
		// ok - build the MRNA generation reaction (forward) -
		local_buffer.append("\t\t<reaction id=\"REACTION");
		
		switch (rxnType)
		{
			case FORWARD_RATE:
				local_buffer.append("_FORWARD");
				break;
			case REVERSE_RATE:
				local_buffer.append("_REVERSE");
				break;
			case CATALYTIC_RATE:
				local_buffer.append("_CATALYTIC");
				break;
			case DEGRADATION:
				local_buffer.append("_DEGRADATION");
				break;
			case GENERATION:
				local_buffer.append("_GENERATION");
				break;
		}
		
		//local_buffer.append(_intReactionCounter);
		local_buffer.append("\" name=\"");
		
		int NUMBER_OF_REACTANTS = arrReactants.size();
		for (int reactant_index=0;reactant_index<NUMBER_OF_REACTANTS;reactant_index++)
		{
			String strTmpSymbol = arrReactants.get(reactant_index);
			
			// Add the symbol to the buffer -
			local_buffer.append(strTmpSymbol);
			
			if (reactant_index!=NUMBER_OF_REACTANTS-1)
			{
				local_buffer.append(" + ");
			}
		}
		
		local_buffer.append(" = ");
		
		// Get the list of products -
		int NUMBER_OF_PRODUCTS = arrProducts.size();
		for (int product_index=0;product_index<NUMBER_OF_PRODUCTS;product_index++)
		{
			String strTmpSymbolRaw = arrProducts.get(product_index);
			local_buffer.append(strTmpSymbolRaw);
			
			if (product_index!=NUMBER_OF_PRODUCTS-1)
			{
				local_buffer.append(" + ");
			}
		}
		
		local_buffer.append("\" reversible=\"false\">\n");
		local_buffer.append("\t\t\t<listOfReactants>\n");
		
		for (int reactant_index=0;reactant_index<NUMBER_OF_REACTANTS;reactant_index++)
		{
			String strTmpSymbol = arrReactants.get(reactant_index);
			local_buffer.append("\t\t\t\t<speciesReference species=\"");
			local_buffer.append(strTmpSymbol);
			local_buffer.append("\" stoichiometry=\"1\"/>\n");	
		}
		
		local_buffer.append("\t\t\t</listOfReactants>\n");
		local_buffer.append("\t\t\t<listOfProducts>\n");
		
		for (int product_index=0;product_index<NUMBER_OF_PRODUCTS;product_index++)
		{
			String strTmpSymbolRaw = arrProducts.get(product_index);			
			local_buffer.append("\t\t\t\t<speciesReference species=\"");
			local_buffer.append(strTmpSymbolRaw);
			local_buffer.append("\" stoichiometry=\"1\"/>\n");
			
		}
		
		local_buffer.append("\t\t\t</listOfProducts>\n");
		local_buffer.append("\t\t</reaction>\n");
		
		// Reaction string -
		String strTmpString = local_buffer.toString();
		
		// Check to see if this string is in the main buffer already -
		if (!buffer.contains(strTmpString))
		{
			// Ok, add the reaction string to the buffer -
			//_arrListSBMLReactions.add(strTmpString);
			
			// Add the local buffer to the main buffer -
			buffer.add(local_buffer.toString());
			
			// Ok, if we keep this reaction then I need to store its type -
			_arrListReactionType.add(rxnType);
		}	
	}
	
	protected ArrayList<String> generateReactionTypeList(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrTypes = new ArrayList<String>();
		StringBuffer local_buffer = new StringBuffer();
		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
    	dbFactory.setNamespaceAware(true);
    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
		
		// Ok, so we are going to turn the reaction list into an xml tree and then use xpath to extract all the species -
		local_buffer.append("<?xml version=\"1.0\"?>\n");
		local_buffer.append("<listOfReactions>\n");
		int NUMBER_OF_REACTIONS = arrListReactions.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = arrListReactions.get(index);
			local_buffer.append(tmpReaction);
		}
		local_buffer.append("</listOfReactions>\n");
		String strTmp = local_buffer.toString();
		
		// Ok, so now we need to create a document -
		Document reaction_dom_tree = dBuilder.parse(new InputSource(new StringReader(local_buffer.toString())));
		String strXPath = "//reaction/@id";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,reaction_dom_tree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strID = tmpNode.getNodeValue();
			
			// put in the arrTypes -
			arrTypes.add(strID);
		}
		
		// return the array -
		return(arrTypes);
	}
	
	protected ArrayList<String> generateSpeciesList(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		StringBuffer local_buffer = new StringBuffer();
		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
    	dbFactory.setNamespaceAware(true);
    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
		
		// Ok, so we are going to turn the reaction list into an xml tree and then use xpath to extract all the species -
		local_buffer.append("<?xml version=\"1.0\"?>\n");
		local_buffer.append("<listOfReactions>\n");
		int NUMBER_OF_REACTIONS = arrListReactions.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = arrListReactions.get(index);
			local_buffer.append(tmpReaction);
		}
		local_buffer.append("</listOfReactions>\n");
		String strTmp = local_buffer.toString();
		
		// Ok, so now we need to create a document -
		Document reaction_dom_tree = dBuilder.parse(new InputSource(new StringReader(local_buffer.toString())));
		
		// Create the XPath -
		String strXPath = "//speciesReference/@species";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,reaction_dom_tree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// check to see if we already have the species in the list -
			if (!arrSpecies.contains(strSpecies))
			{
				arrSpecies.add(strSpecies);
			}
		}
		
		// return the array -
		return(arrSpecies);
	}
	
}
