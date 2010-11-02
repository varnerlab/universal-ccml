package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MAMAPKHandler extends CCMLMAObject implements IReactionHandler {
	// Class/attributes -
	
	
	private void init(Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
		
		// Get the siganling components -
		strXPathBase = "//signaling_block[@block_class='MAPK_NETWORK']/listOfSignalingComponents/signaling_component";
		populateProperties(strXPathBase,ccmlTree);
		
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
		
		// Get the list of prefixes -
		String strPrefixXPath = "//listOfSymbolPrefixes/symbol_prefix";
		populateProperties(strPrefixXPath,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize -
		init(ccmlTree);
		
		// Build the MAPK cascade -
		buildMAPKActivation(arrRxnList,ccmlTree);
		
		// Build the PASE interactions -
		buildPASEReactions(arrRxnList,ccmlTree);
		
		// Build the ERK feedback -
		buildERKFeedback(arrRxnList,ccmlTree);
		
		// Process the interface block -
		buildInterfaceReactions("MAPK_NETWORK",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	private void buildERKFeedback(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrInitiators = new ArrayList<String>();
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strERKSymbol = (String)getProperty("ERK_SYMBOL");

		// Get prefix -
		String strPhosphorylationPrefix = (String)getProperty("PHOSPHORYLATION_PREFIX");
		String strDoublePhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strActivatedPrefix = (String)getProperty("ACTIVATED_PREFIX");
		
		// Get the list of initiator species -
		String strInitXPath = "//signaling_block[@block_class='MAPK_NETWORK']/listOfInitiators/initiator/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strInitXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_INITIATORS= nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_INITIATORS;index++)
		{
			// Get the key -
			Node tmpNode = nodeList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// store these in the array -
			arrInitiators.add(strKeyName);
		}
		
		// Ok, so now - let's go through, bind ERK_PP and then knock off SOS?
		for (int index=0;index<NUMBER_OF_INITIATORS;index++)
		{		
			// Get the initiator -
			String strReceptor = arrInitiators.get(index);
			
			// Find the *last* _
			int INDEX_LAST_UNDERSCORE = strReceptor.lastIndexOf("_");
			
			// Get the initiator -
			String strDeactivatedReceptor = strReceptor.substring(0,INDEX_LAST_UNDERSCORE);
			String strFragment = strReceptor.substring(INDEX_LAST_UNDERSCORE+1,strReceptor.length());
			
			// encode -- binding -
			arrReactants.add(strReceptor+"_"+strReceptorCompartment);
			arrReactants.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strReceptor+"_"+strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Deactivate the receptor -
			arrReactants.add(strReceptor+"_"+strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strDeactivatedReceptor+"_"+strReceptorCompartment);
			arrProducts.add(strFragment+"_"+strReceptorCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	private void buildPASEReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		// Ok, get the PASE symbols -
		String strRAFPASESymbol = (String)getProperty("RAF_PASE_SYMBOL");
		String strMEKPASESymbol = (String)getProperty("MEK_PASE_SYMBOL");
		String strERKPASESymbol = (String)getProperty("ERK_PASE_SYMBOL");
		
		// Get the MAPK symbols -
		String strRAFSymbol = (String)getProperty("RAF_SYMBOL");
		String strMEKSymbol = (String)getProperty("MEK_SYMBOL");
		String strERKSymbol = (String)getProperty("ERK_SYMBOL");
		
		// Get prefix -
		String strPhosphorylationPrefix = (String)getProperty("PHOSPHORYLATION_PREFIX");
		String strDoublePhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strActivatedPrefix = (String)getProperty("ACTIVATED_PREFIX");
		
		// Code the interactions -
		// Binding -
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strRAFPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strRAFPASESymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// De-phosphorylate that bitch -
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strRAFPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strRAFPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strRAFSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// MEK ----------------- //
		// Binding - MEK
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strMEKPASESymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// De-phosphorylate that bitch -
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Binding -
		arrReactants.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strMEKPASESymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// De-phosphorylate that bitch -
		arrReactants.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strMEKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// ERK ----------------- //
		// Binding - ERK
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strERKPASESymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// De-phosphorylate that bitch -
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Binding -
		arrReactants.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strERKPASESymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// De-phosphorylate that bitch -
		arrReactants.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strERKPASESymbol+"_"+strReceptorCompartment);
		arrProducts.add(strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();	
	}
	
	private void buildInterfaceReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		// Get prefix -
		String strPhosphorylationPrefix = (String)getProperty("PHOSPHORYLATION_PREFIX");
		String strDoublePhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strActivatedPrefix = (String)getProperty("ACTIVATED_PREFIX");
		
		// Ok, let's get the list interface blocks -
		String strBlockXPath = "//signaling_block[@block_class='MAPK_NETWORK']/listOfInterfaces/interface/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strBlockXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_INTERFACES= nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_INTERFACES;index++)
		{
			// Get the interface symbol -
			Node tmpNode = nodeList.item(index);
			String strInterfaceSymbol = tmpNode.getNodeValue();
			
			// Ok now that I have the interface symbol - process the targets (activate) 
			String strTargetXPath = "//signaling_block[@block_class='MAPK_NETWORK']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_activate/@symbol";
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
			
			String strTargetDeactivateXPath = "//signaling_block[@block_class='MAPK_NETWORK']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_deactivate/@symbol";
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
	
	private void buildMAPKActivation(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrInitiators = new ArrayList<String>();
		
		// Get the symbols for the MAPK components -
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strRASSymbol = (String)getProperty("RAS_SYMBOL");
		String strRAFSymbol = (String)getProperty("RAF_SYMBOL");
		String strMEKSymbol = (String)getProperty("MEK_SYMBOL");
		String strERKSymbol = (String)getProperty("ERK_SYMBOL");
		String strGTPSymbol = (String)getProperty("GTP_SYMBOL");
		String strGDPSymbol = (String)getProperty("GDP_SYMBOL");
		
		// Get prefix -
		String strPhosphorylationPrefix = (String)getProperty("PHOSPHORYLATION_PREFIX");
		String strDoublePhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strActivatedPrefix = (String)getProperty("ACTIVATED_PREFIX");
		
		// Get the list of initiator species -
		String strInitXPath = "//signaling_block[@block_class='MAPK_NETWORK']/listOfInitiators/initiator/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strInitXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_INITIATORS= nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_INITIATORS;index++)
		{
			// Get the key -
			Node tmpNode = nodeList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// store these in the array -
			arrInitiators.add(strKeyName);
		}
			
		// Ok, so we need to have RAS bind with GDP -
		arrReactants.add(strRASSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strGDPSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strRASSymbol+"_"+strGDPSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Ok, so we need this reaction to repeat for every receptor -
		for (int receptor_index=0;receptor_index<NUMBER_OF_INITIATORS;receptor_index++)
		{
			// Get the initiator symbol -
			String strActiveMAPKReceptor = arrInitiators.get(receptor_index);
			
			// Activation of RAS by the activated receptor -
			arrReactants.add(strActiveMAPKReceptor+"_"+strReceptorCompartment);
			arrReactants.add(strRASSymbol+"_"+strGDPSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strActiveMAPKReceptor+"_"+strRASSymbol+"_"+strGDPSymbol+"_"+strReceptorCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Exchange of GDP w/GTP -
			arrReactants.add(strActiveMAPKReceptor+"_"+strRASSymbol+"_"+strGDPSymbol+"_"+strReceptorCompartment);
			arrReactants.add(strGTPSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strActiveMAPKReceptor+"_"+strReceptorCompartment);
			arrProducts.add(strRASSymbol+"_"+strGTPSymbol+"_"+strReceptorCompartment);
			arrProducts.add(strGDPSymbol+"_"+strReceptorCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
		}
		
		// Activation of RAF by RAS - binding
		arrReactants.add(strRASSymbol+"_"+strGTPSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strRAFSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strRASSymbol+"_"+strGTPSymbol+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Activation of RAF by RAS - phosphorylation -
		arrReactants.add(strRASSymbol+"_"+strGTPSymbol+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strRASSymbol+"_"+strGTPSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// First phospphorylation of MEK by RAS - binding
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// First MEK phosphorylation by RAS - catalytic
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Second phospphorylation of MEK - binding
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Second MEK phosphorylation by RAS - catalytic
		arrReactants.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strPhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strActivatedPrefix+"_"+strRAFSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// First phospphorylation of ERK by MEK - binding
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strERKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// First ERK phosphorylation by MEK - catalytic
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strERKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Second phospphorylation of ERK - binding
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrReactants.add(strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Second ERK phosphorylation by MEK - catalytic
		arrReactants.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strPhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strMEKSymbol+"_"+strReceptorCompartment);
		arrProducts.add(strDoublePhosphorylationPrefix+"_"+strERKSymbol+"_"+strReceptorCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// 
	}
	
	


}
