package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MANotchSignalingHandler extends CCMLMAObject implements IReactionHandler {

	
	private void init(Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
		
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
		
		// Load the siganling components -
		String strXPathSignalingComponents = "//signaling_block[@block_class='NOTCH_NETWORK']/listOfSignalingComponents/signaling_component";
		populateProperties(strXPathSignalingComponents,ccmlTree);
		
		// Get the initiator symbols -
		String strXPathSignalingIniators = "//signaling_block[@block_class='NOTCH_NETWORK']/listOfInitiators/initiator";
		populateProperties(strXPathSignalingIniators,ccmlTree);
	}
	

	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler -
		init(ccmlTree);
		
		// Encode the NCID regulation -
		encodeNCIDActivation(arrRxnList,ccmlTree);
		
		// Encode ACTIVE_NCID UB -
		encodeActNCIDUbReactions(arrRxnList,ccmlTree);
		
		// Encode DELTEX_NEDD4 activity -
		encodeDELTEXNEDD4Activity(arrRxnList,ccmlTree);
		
		// Encode degradation reactions -
		buildDegradationReactions(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	// Encode DELTEX_NEDD4 activity -
	private void encodeDELTEXNEDD4Activity(ArrayList<String >arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strDELTEXSymbol = (String)getProperty("DELTEX_NEDD4_SYMBOL");
		String strDELTA4NCED = (String)getProperty("DELTA4_NCED_SYMBOL");
		String strJAG1NCEDComplex = (String)getProperty("JAG1_NCED_SYMBOL");
		String strJag1UBSymbol = (String)getProperty("UB_JAG1_NCED_SYMBOL");
		String strDelta4UBSymbol = (String)getProperty("UB_DELTA4_NCED_SYMBOL");
		
		// Encode tmp lists -
		ArrayList<String> tmpReactants = new ArrayList<String>();
		tmpReactants.add(strJAG1NCEDComplex);
		tmpReactants.add(strDELTA4NCED);
		
		ArrayList<String> tmpProducts = new ArrayList<String>();
		tmpProducts.add(strJag1UBSymbol);
		tmpProducts.add(strDelta4UBSymbol);
		
		// How many reactants?
		int NUMBER_OF_REACTANTS = tmpReactants.size();
		for (int index=0;index<NUMBER_OF_REACTANTS;index++)
		{
			// Get the symbol -
			String strTmpReactant = tmpReactants.get(index);
			String strTmpProduct = tmpProducts.get(index);
			
			
			
			// Populate the lists of arrR and arrP lists -
			// Binding -
			arrReactants.add(strTmpReactant+"_"+strCompartment);
			arrReactants.add(strDELTEXSymbol+"_"+strCompartment);
			arrProducts.add(strTmpReactant+"_"+strDELTEXSymbol+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Reaction -
			arrReactants.add(strTmpReactant+"_"+strDELTEXSymbol+"_"+strCompartment);
			arrProducts.add(strDELTEXSymbol+"_"+strCompartment);
			arrProducts.add(strTmpProduct+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	// Log to add Ub groups to ACT_NCID -
	private void encodeActNCIDUbReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strActNCIDSymbol = (String)getProperty("ACT_NCID_SYMBOL");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strSEL1Symbol = (String)getProperty("SEL1_SYMBOL");
		String strUBActNCIDSymbol = (String)getProperty("UB_ACT_NCID_SYMBOL");
		
		// Binding of SEL1 to ACT_NCID -
		arrReactants.add(strActNCIDSymbol+"_"+strCompartment);
		arrReactants.add(strSEL1Symbol+"_"+strCompartment);
		arrProducts.add(strActNCIDSymbol+"_"+strSEL1Symbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Addition of UB grps -
		arrReactants.add(strActNCIDSymbol+"_"+strSEL1Symbol+"_"+strCompartment);
		arrProducts.add(strSEL1Symbol+"_"+strCompartment);
		arrProducts.add(strUBActNCIDSymbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}
	
	// Logic to build the degradation -
	private void buildDegradationReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//signaling_block[@block_class='NOTCH_NETWORK']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//signaling_block[@block_class='NOTCH_NETWORK']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(ccmlTree,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	private void encodeNCIDActivation(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		
		// Ok, so let's party ...
		String strJAGComplex = (String)getProperty("JAG1_NOTCH1_TACE_COMPLEX");
		String strDeltaComplex = (String)getProperty("DELTA4_NOTCH1_TACE_COMPLEX");
		String strCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strJAG1NCEDComplex = (String)getProperty("JAG1_NCED_SYMBOL");
		String strDelta4NCEDComplex = (String)getProperty("DELTA4_NCED_SYMBOL");
		String strActNCIDSymbol = (String)getProperty("ACT_NCID_SYMBOL");
		
		// Generation of active JAG1 NCID -
		arrReactants.add(strJAGComplex+"_"+strCompartment);
		arrProducts.add(strJAG1NCEDComplex+"_"+strCompartment);
		arrProducts.add(strActNCIDSymbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Generation of active DELTA4 NCID -
		arrReactants.add(strDeltaComplex+"_"+strCompartment);
		arrProducts.add(strDelta4NCEDComplex+"_"+strCompartment);
		arrProducts.add(strActNCIDSymbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
	}
}
