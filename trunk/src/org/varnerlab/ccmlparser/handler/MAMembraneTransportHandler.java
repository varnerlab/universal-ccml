package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MAMembraneTransportHandler extends CCMLMAObject implements IReactionHandler {

	// Override the parent init method -
	private void init(Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
		
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
		
		// List of prefixes -
		String strPrefixXPath = "//listOfSymbolPrefixes/symbol_prefix";
		populateProperties(strPrefixXPath,ccmlTree);
	}
	
	
	// Ok, this guy processes the membrane transport block -
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
	
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize -
		init(ccmlTree);
		
		// Build the list -
		processTransportReactions(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
		
	}
	
	// Build the reactions for transport - simple, on, off, transport -
	private void processTransportReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Grab the transport symbols -
		String strExport = (String)getProperty("EXTRACELLULAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("EXTRACELLULAR_TRANSPORTER_IN");
		
		// Ok, let's get the list interface blocks -
		String strBlockXPath = "//Membrane_transport_block[@block_class='DEFAULT_MEMBRANE_TRANSPORT']/listOfInterfaces/interface/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strBlockXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_INTERFACES= nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_INTERFACES;index++)
		{
			// Get the interface symbol -
			Node tmpNode = nodeList.item(index);
			String strInterfaceSymbol = tmpNode.getNodeValue();
			
			// What compartment is this interface symbol in?
			String strISymbolComparmentXPath = "//Membrane_transport_block[@block_class='DEFAULT_MEMBRANE_TRANSPORT']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/@compartment";
			String strCompartment = queryCCMLTree(ccmlTree,strISymbolComparmentXPath);
			
			// Ok now that I have the interface symbol - process the targets (activate) 
			String strTargetXPath = "//Membrane_transport_block[@block_class='DEFAULT_MEMBRANE_TRANSPORT']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_transport/@symbol";
			NodeList nodeTargetList = (NodeList)_xpath.evaluate(strTargetXPath,ccmlTree,XPathConstants.NODESET);
			int NUMBER_OF_TARGETS= nodeTargetList.getLength();
			for (int target_index=0;target_index<NUMBER_OF_TARGETS;target_index++)
			{
				// Get the target symbol -
				Node tmpTargetNode = nodeTargetList.item(target_index);
				String strTargetSymbol = tmpTargetNode.getNodeValue();
				
				// Have to do just two more xpath hits to get the to and from compartments -
				// to compartment 
				String strToCompartmentXPath = "//Membrane_transport_block[@block_class='DEFAULT_MEMBRANE_TRANSPORT']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_transport[@symbol='"+strTargetSymbol+"']/@to_compartment";
				String strToCompartment = queryCCMLTree(ccmlTree,strToCompartmentXPath);
				
				// From compartment -
				String strFromCompartmentXPath = "//Membrane_transport_block[@block_class='DEFAULT_MEMBRANE_TRANSPORT']/listOfInterfaces/interface[@symbol='"+strInterfaceSymbol+"']/target_transport[@symbol='"+strTargetSymbol+"']/@from_compartment";
				String strFromCompartment = queryCCMLTree(ccmlTree,strFromCompartmentXPath);
				
				// Process the targets -
				// Binding -
				arrReactants.add(strInterfaceSymbol+"_"+strCompartment);
				arrReactants.add(strTargetSymbol+"_"+strFromCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Transport that bitch ..
				arrReactants.add(strInterfaceSymbol+"_"+strTargetSymbol+"_"+strCompartment);
				arrProducts.add(strInterfaceSymbol+"_"+strCompartment);
				arrProducts.add(strTargetSymbol+"_"+strToCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
				arrReactants.clear();
				arrProducts.clear();
			}
		}
	}
}
