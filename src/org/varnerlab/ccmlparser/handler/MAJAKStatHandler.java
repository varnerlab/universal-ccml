package org.varnerlab.ccmlparser.handler;

import java.io.InputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.varnerlab.ccmlparser.IReceptorNetworkHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public class MAJAKStatHandler implements IReceptorNetworkHandler {
	// Class/instance attributes -
	private Hashtable<String,Object> _propTable = new Hashtable<String,Object>();
	private XPathFactory  _xpFactory = XPathFactory.newInstance();
	private XPath _xpath = _xpFactory.newXPath();
	private ArrayList<String> _arrListSBMLReactions = new ArrayList<String>();
	private int _intReactionCounter = 0;
	private ArrayList<ReactionType> _arrListReactionType = new ArrayList<ReactionType>();

	// Set and get properties -
	public void setProperty(String key,Object value)
	{
		_propTable.put(key, value);
	}
	
	public Object getProperty(String key)
	{
		return(_propTable.get(key));
	}
	
	public void setReactionIndex(int index)
	{
		_intReactionCounter = index;
	}
	

	private void initHandler(Document ccmlTree) throws Exception
	{
		// Ok, get the name of the receptor -
		String strReceptorXPath = "//receptor_block[@block_class='JAK_STAT']/@symbol";
		String strReceptor = queryCCMLTree(ccmlTree,strReceptorXPath);
		setProperty("RECEPTOR_SYMBOL",strReceptor);
		
		// OK, get the name of the ligand -
		String strLigandXPath = "//receptor_block[@block_class='JAK_STAT']/listOfLigands/ligand/@symbol";
		String strLigand = queryCCMLTree(ccmlTree,strLigandXPath);
		setProperty("LIGAND_SYMBOL",strLigand);
		
		// OK, get the name of the ligand -
		String strLigandCXPath = "//receptor_block[@block_class='JAK_STAT']/listOfLigands/ligand/@compartment";
		String strLigandC = queryCCMLTree(ccmlTree,strLigandCXPath);
		setProperty("LIGAND_COMPARTMENT",strLigandC);
		
		// Get the compartment of the receptor -
		String strReceptorCompartmentXPath = "//receptor_block[@block_class='JAK_STAT']/@compartment";
		String strReceptorCompartment = queryCCMLTree(ccmlTree,strReceptorCompartmentXPath);
		setProperty("RECEPTOR_COMPARTMENT",strReceptorCompartment);
		
		// Get the adaptor symbol -
		String strAdaptorXPath = "//receptor_block[@block_class='JAK_STAT']/listOfAdapters/adapter/@symbol";
		String strAdapterSymbol = queryCCMLTree(ccmlTree,strAdaptorXPath);
		setProperty("ADAPTER_SYMBOL",strAdapterSymbol);
		
		// Get the messenger symbol -
		String strMessengerXPath = "//receptor_block[@block_class='JAK_STAT']/listOfMessengers/messenger/@symbol";
		String strMessengerSymbol = queryCCMLTree(ccmlTree,strMessengerXPath);	
		setProperty("MESSENGER_SYMBOL",strMessengerSymbol);
		
		// Ok, we need to process the list of regulators and store the info in properties -
		String strXPath = "//Receptor_network_block/listOfReceptors/receptor_block[@block_class='JAK_STAT']/listOfRegulators/regulator/@key";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_REGULATORS = nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_REGULATORS;index++)
		{
			// Get the key -
			Node tmpNode = nodeList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = "//Receptor_network_block/listOfReceptors/receptor_block[@block_class='JAK_STAT']/listOfRegulators/regulator[@key='"+strKeyName+"']/@symbol";
			String strSymbol = queryCCMLTree(ccmlTree,strSymbolXPath);
			
			// store in prop -
			setProperty(strKeyName,strSymbol);
		}
	}
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrListLocal = new ArrayList<String>();
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize - get specific names from tree 
		initHandler(ccmlTree);
	
		// Ok, so here is the specific logic for this network -
		buildReceptorBindingReactions(arrRxnList,ccmlTree);
		buildReceptorActivation(arrRxnList,ccmlTree);
		buildStat3Activation(arrRxnList,ccmlTree);
		buildCytosolicPIAS3Binding(arrRxnList,ccmlTree);
		buildReceptorSOCS3Binding(arrRxnList,ccmlTree);
		buildBackgroundMessengerRegulation(arrRxnList,ccmlTree);
		buildDegradationReactions(arrRxnList,ccmlTree);
		
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);
	}
	
	private void buildDegradationReactions(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//receptor_block[@block_class='JAK_STAT']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//receptor_block[@block_class='JAK_STAT']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(doc,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	private ArrayList<String> generateSpeciesList(ArrayList<String> arrListReactions,Document doc) throws Exception
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
	
	private void buildBackgroundMessengerRegulation(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrFinalProducts = new ArrayList<String>();
		
		// Get some stuff from properties -
		String strMessenger = (String)getProperty("MESSENGER_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strMessengerKinase = (String)getProperty("MESSENGER_KINASE");
		String strMessengerPASE = (String)getProperty("MESSENGER_PASE");
		
		// Ok, we need to encode the activation and deactivation of the messenger -
		arrReactants.add(strMessenger+"_2P_"+strCompartment);
		arrReactants.add(strMessengerKinase+"_"+strCompartment);
		arrProducts.add(strMessenger+"_2P_"+strMessengerKinase+"_"+strCompartment);
		arrFinalProducts.add(strMessenger+"_4P_"+strCompartment);
		arrFinalProducts.add(strMessengerKinase+"_"+strCompartment);
		encodeMassActionSBMLCatalyticReaction(arrListReactions,doc,arrReactants,arrProducts,arrFinalProducts);
		arrReactants.clear();
		arrProducts.clear();
		arrFinalProducts.clear();
		
		// Deactivation -
		arrReactants.add(strMessenger+"_4P_"+strCompartment);
		arrReactants.add(strMessengerPASE+"_"+strCompartment);
		arrProducts.add(strMessenger+"_4P_"+strMessengerPASE+"_"+strCompartment);
		arrFinalProducts.add(strMessenger+"_2P_"+strCompartment);
		arrFinalProducts.add(strMessengerPASE+"_"+strCompartment);
		encodeMassActionSBMLCatalyticReaction(arrListReactions,doc,arrReactants,arrProducts,arrFinalProducts);
		arrReactants.clear();
		arrProducts.clear();
		arrFinalProducts.clear();
		
		// Add a reaction to de-phosphorylate the monomer -
		arrReactants.add(strMessenger+"_P_"+strCompartment);
		arrReactants.add(strMessengerPASE+"_"+strCompartment);
		arrProducts.add(strMessenger+"_P_"+strMessengerPASE+"_"+strCompartment);
		arrFinalProducts.add(strMessenger+"_"+strCompartment);
		arrFinalProducts.add(strMessengerPASE+"_"+strCompartment);
		encodeMassActionSBMLCatalyticReaction(arrListReactions,doc,arrReactants,arrProducts,arrFinalProducts);
		arrReactants.clear();
		arrProducts.clear();
		arrFinalProducts.clear();
	}
	
	private void buildReceptorSOCS3Binding(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get some stuff from properties list --
		String strComplex = (String)getProperty("ACTIVATED_RECEPTOR_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strRepressor = (String)getProperty("RECEPTOR_BINDING");
		
		// Ok, encode reversible receptor binding -
		arrReactants.add(strComplex+"_"+strCompartment);
		arrReactants.add(strRepressor+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strRepressor+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
	}
	
	private void buildCytosolicPIAS3Binding(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get some stuff from the properties list -
		String strMessenger = (String)getProperty("MESSENGER_SYMBOL");
		String strRepressor = (String)getProperty("CYTOSOL_MESSENGER_BINDING");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		// Ok, encode reversible messenger binding -
		arrReactants.add(strMessenger+"_2P_"+strCompartment);
		arrReactants.add(strRepressor+"_"+strCompartment);
		arrProducts.add(strMessenger+"_2P_"+strRepressor+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		
	}
	
	private void buildStat3Activation(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrFinalProducts = new ArrayList<String>();
		
		// Get some stuff from the properties list -
		String strComplex = (String)getProperty("ACTIVATED_RECEPTOR_SYMBOL");
		String strAdapter = (String)getProperty("ADAPTER_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strKinase = (String)getProperty("RECEPTOR_KINASE");
		String strPASE = (String)getProperty("RECEPTOR_PASE");
		String strMessenger = (String)getProperty("MESSENGER_SYMBOL");
		
		// Binding of first STAT3 -
		arrReactants.add(strComplex+"_"+strCompartment);
		arrReactants.add(strMessenger+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strMessenger+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Binding of second STAT3 -
		String strFirstStat3 = strComplex+"_"+strMessenger;
		arrReactants.add(strFirstStat3+"_"+strCompartment);
		arrReactants.add(strMessenger+"_"+strCompartment);
		arrProducts.add(strFirstStat3+"_2_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Standard activation (binding and activation step) -
		arrReactants.add(strFirstStat3+"_2_"+strCompartment);
		arrProducts.add(strFirstStat3+"_2P_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
			
		// Encode the release of first Stat3_2P -
		arrReactants.add(strFirstStat3+"_2P_"+strCompartment);
		arrProducts.add(strMessenger+"_P_"+strCompartment);
		arrProducts.add(strFirstStat3+"_P_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Encode the release of the second Stat3_P -
		arrReactants.add(strFirstStat3+"_P_"+strCompartment);
		arrProducts.add(strMessenger+"_P_"+strCompartment);
		arrProducts.add(strComplex+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// OK, so now we need to dimerize the stat3 -
		String strReactant = strMessenger+"_P_"+strCompartment;
		String strProduct = strMessenger+"_2P_"+strCompartment;
		encodeMassActionSBMLDimerization(arrListReactions,strReactant,"2",strProduct,"1",ReactionType.FORWARD_RATE);
		encodeMassActionSBMLDimerization(arrListReactions,strProduct,"1",strReactant,"2",ReactionType.REVERSE_RATE);
		
	}
	
	private void buildReceptorActivation(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrFinalProducts = new ArrayList<String>();
		
		// Get some stuff from the properties list -
		String strComplex = (String)getProperty("RECEPTOR_LIGAND_COMPLEX_SYMBOL");
		String strAdapter = (String)getProperty("ADAPTER_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strKinase = (String)getProperty("RECEPTOR_KINASE");
		String strPASE = (String)getProperty("RECEPTOR_PASE");
		
		// Ok, construct the reactants and products -
		arrReactants.add(strComplex+"_"+strCompartment);
		arrReactants.add(strAdapter+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strAdapter+"_"+strCompartment);
		
		// Encode the forward and reverse rates - of adapter binding 
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		
		// Ok, so we need to encode the *activation* and *deactivation* of the receptor -
		arrReactants.clear();
		arrProducts.clear();
		
		// Ok, we need to bind *two* JAKS -
		// JAK 1
		arrReactants.add(strComplex+"_"+strAdapter+"_"+strCompartment);
		arrReactants.add(strKinase+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strAdapter+"_"+strKinase+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// JAK 2
		String strFirstJAK = strComplex+"_"+strAdapter+"_"+strKinase;
		arrReactants.add(strFirstJAK+"_"+strCompartment);
		arrReactants.add(strKinase+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strAdapter+"_"+strKinase+"_2_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Autophopshporylation of JAKs
		String strSecondJAK = strComplex+"_"+strAdapter+"_"+strKinase+"_2";
		arrReactants.add(strSecondJAK+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strAdapter+"_"+strKinase+"_2P_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylation of the receptor -
		String strActivatedJAKs = strComplex+"_"+strAdapter+"_"+strKinase+"_2P";
		arrReactants.add(strActivatedJAKs+"_"+strCompartment);
		arrProducts.add(strComplex+"_"+strAdapter+"_2P_"+strKinase+"_2P_"+strCompartment);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Set the symbol of the activated complex -
		setProperty("ACTIVATED_RECEPTOR_SYMBOL",strComplex+"_"+strAdapter+"_2P_"+strKinase+"_2P");
	}
	
	
	private void encodeMassActionSBMLCatalyticReaction(ArrayList<String> arrListReactions,Document doc,ArrayList<String> arrReactants,ArrayList<String> arrProducts,ArrayList<String> arrFinalProducts) throws Exception
	{
		
		// forward and reverse on - binding step -
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		
		// catalytic step -
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrFinalProducts,ReactionType.CATALYTIC_RATE);
	}
	
	// Ok, so here we go - L binding with receptor -
	private void buildReceptorBindingReactions(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Ok, get the receptor and ligand symbols -
		String strReceptorSymbol = (String)getProperty("RECEPTOR_SYMBOL");
		String strLigandSymbol = (String)getProperty("LIGAND_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strLCompartment = (String)getProperty("LIGAND_COMPARTMENT");
		
   		// Encode ligand-receptor binding -
		arrReactants.add(strLigandSymbol+"_"+strLCompartment);
		arrReactants.add(strReceptorSymbol+"_"+strCompartment);
		arrProducts.add(strLigandSymbol+"_"+strReceptorSymbol+"_"+strCompartment);
		
		// Encode the forward and reverse rates -
		encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrListReactions,doc,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		
		// Ok, we'll need the complex symbol for later -
		String strComplex = strLigandSymbol+"_"+strReceptorSymbol;
		setProperty("RECEPTOR_LIGAND_COMPLEX_SYMBOL",strComplex);
	}
	
	
	private void encodeMassActionSBMLDimerization(ArrayList<String> buffer,String strReactant,String strCoeffReactant,String strProduct,String strCoeffProduct,ReactionType rxnType) throws Exception
	{
		// Method attributes -
		StringBuffer local_buffer = new StringBuffer();
		
		
		// ok - build the MRNA generation reaction (forward) -
		local_buffer.append("\t\t<reaction id=\"R_");
		local_buffer.append(_intReactionCounter);
		local_buffer.append("\" name=\""+strCoeffReactant+"*");
		local_buffer.append(strReactant);
		local_buffer.append(" = "+strCoeffProduct+"*");
		local_buffer.append(strProduct);
		local_buffer.append("\" reversible=\"false\">\n");
		local_buffer.append("\t\t\t<listOfReactants>\n");
		local_buffer.append("\t\t\t\t<speciesReference species=\"");
		local_buffer.append(strReactant);
		local_buffer.append("\" stoichiometry=\"");
		local_buffer.append(strCoeffReactant);
		local_buffer.append("\"/>\n");	
		local_buffer.append("\t\t\t</listOfReactants>\n");
		local_buffer.append("\t\t\t<listOfProducts>\n");
		local_buffer.append("\t\t\t\t<speciesReference species=\"");
		local_buffer.append(strProduct);
		local_buffer.append("\" stoichiometry=\"");
		local_buffer.append(strCoeffProduct);
		local_buffer.append("\"/>\n");	
		local_buffer.append("\t\t\t</listOfProducts>\n");
		local_buffer.append("\t\t</reaction>\n");
		
		// Reaction string -
		String strTmpString = local_buffer.toString();
		
		// Check to see if this string is in the main buffer already -
		if (!_arrListSBMLReactions.contains(strTmpString))
		{
			// Ok, add the reaction string to the buffer -
			_arrListSBMLReactions.add(strTmpString);
			
			// Add the local buffer to the main buffer -
			buffer.add(local_buffer.toString());
			
			// Ok, if we keep this reaction then I need to store its type -
			_arrListReactionType.add(rxnType);
			
			// update the counter - 
			_intReactionCounter++;
		}	
	}
	
	private void encodeMassActionSBMLReaction(ArrayList<String> buffer,Document doc,ArrayList<String> arrReactants,ArrayList<String> arrProducts,ReactionType rxnType) throws Exception
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
		if (!_arrListSBMLReactions.contains(strTmpString))
		{
			// Ok, add the reaction string to the buffer -
			_arrListSBMLReactions.add(strTmpString);
			
			// Add the local buffer to the main buffer -
			buffer.add(local_buffer.toString());
			
			// Ok, if we keep this reaction then I need to store its type -
			_arrListReactionType.add(rxnType);
			
			// update the counter - 
			_intReactionCounter++;
		}	
	}
	
	
	
	// Get a string -
	private String queryCCMLTree(Document ccmlTree,String strXPath)
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
}
