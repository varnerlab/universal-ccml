package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MABCATENINNetworkHandler extends CCMLMAObject implements IReactionHandler {

	
	
	private void init(Document ccmlTree) throws Exception
	{
		// Get the global symbols -
		String strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize -
		init(ccmlTree);
		
		// Implement -
		encodeBCATENINBinding(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}

	private void encodeBCATENINBinding(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		
		// Ok, we need to go through all the SMAD blocks -
		String strSMADBlocksXPath = "//signaling_block[@block_class='BCATENIN_NETWORK']/@symbol";
		NodeList nodeSMADBlocksList = (NodeList)_xpath.evaluate(strSMADBlocksXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_BLOCKS = nodeSMADBlocksList.getLength();
		for (int block_index=0;block_index<NUMBER_OF_BLOCKS;block_index++)
		{
			// Initialize the ArrayLists -
			ArrayList<String> arrReactants = new ArrayList<String>();
			ArrayList<String> arrProducts = new ArrayList<String>();
			ArrayList<String> arrInitiators = new ArrayList<String>();
			ArrayList<String> arrSignalingComponents = new ArrayList<String>();
			
			// Get the key -
			Node tmpBlockNode = nodeSMADBlocksList.item(block_index);
			String strBlockName = tmpBlockNode.getNodeValue();
			
			// Get the list of initiator species -
			String strInitXPath = "//signaling_block[@symbol='"+strBlockName+"']/listOfInitiators/initiator/@symbol";
			System.out.println(strInitXPath);
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
			
			// Ok, get the list of siganling components -
			String strSigXPath = "//signaling_block[@symbol='"+strBlockName+"']/listOfSignalingComponents/signaling_component/@symbol";
			NodeList nodeSigList = (NodeList)_xpath.evaluate(strSigXPath,ccmlTree,XPathConstants.NODESET);
			int NUMBER_OF_SIGNALS= nodeSigList.getLength();
			for (int index = 0;index<NUMBER_OF_SIGNALS;index++)
			{
				// Get the key -
				Node tmpNode = nodeSigList.item(index);
				String strKeyName = tmpNode.getNodeValue();
				
				// store these in the array -
				arrSignalingComponents.add(strKeyName);
			}
			
			// OK, when I get here I will cross the number of initiators with the signals -
			for (int initiator_index=0;initiator_index<NUMBER_OF_INITIATORS;initiator_index++)
			{
				
				// Get the inititar symbol =
				String strInitSymbol = arrInitiators.get(initiator_index);
				
				for (int signal_index=0;signal_index<NUMBER_OF_SIGNALS;signal_index++)
				{
					// Get the signal -
					String strSignalSymbol = arrSignalingComponents.get(signal_index);
					
					// Encode the reaction ... bitches..
					// Binding -
					arrReactants.add(strInitSymbol+"_"+strReceptorCompartment);
					arrReactants.add(strSignalSymbol+"_"+strReceptorCompartment);
					arrProducts.add(strInitSymbol+"_"+strSignalSymbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
					encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
					arrReactants.clear();
					arrProducts.clear();
				}		
			}
		}
	}
}
