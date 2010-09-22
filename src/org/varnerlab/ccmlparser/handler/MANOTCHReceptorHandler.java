package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.IReceptorNetworkHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.varnerlab.ccmlparser.CCMLMAObject;

public class MANOTCHReceptorHandler extends CCMLMAObject implements IReceptorNetworkHandler {

	
	private void init(Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
		
		
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception 
	{
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler -
		init(ccmlTree);
		
		// Generate the notch pathway -
		buildNotchPathway(arrRxnList,ccmlTree);
		
		// Generate the degradation reactions -
		buildDegradationReactions(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	// Logic to build the degradation -
	private void buildDegradationReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//receptor_block[@block_class='NOTCH_RECEPTOR']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//receptor_block[@block_class='NOTCH_RECEPTOR']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(ccmlTree,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
	
	// Logic for the notch pathway -
	private void buildNotchPathway(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		String strLigandCompartment = (String)getProperty("LIGAND_COMPARTMENT");
		
		// Get the receptor list -
		String strReceptorXPath = "//receptor_block[@block_class='NOTCH_RECEPTOR']/@symbol";
		NodeList nodeReceptorList = (NodeList)_xpath.evaluate(strReceptorXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_RECEPTORS = nodeReceptorList.getLength();
		for (int receptor_index=0;receptor_index<NUMBER_OF_RECEPTORS;receptor_index++)
		{
			// Get the receptor -
			Node tmpNodeReceptor = nodeReceptorList.item(receptor_index);
			String strReceptor = tmpNodeReceptor.getNodeValue();
		
			// Process the ligand binding with the receptor -
			String strLigandXPath = "//receptor_block[@symbol='"+strReceptor+"']/listOfLigands/ligand/@symbol";
			NodeList nodeList = (NodeList)_xpath.evaluate(strLigandXPath,ccmlTree,XPathConstants.NODESET);
			int NUMBER_OF_LIGANDS = nodeList.getLength();
			for (int ligand_index=0;ligand_index<NUMBER_OF_LIGANDS;ligand_index++)
			{
				// Get the ligand -
				Node tmpNode = nodeList.item(ligand_index);
				String strLigandSymbol = tmpNode.getNodeValue();
				
				// Encode the binding of the ligand with the receptor -
				arrReactants.add(strLigandSymbol+"_"+strLigandCompartment);
				arrReactants.add(strReceptor+"_"+strReceptorCompartment);
				arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strReceptorCompartment);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
				arrReactants.clear();
				arrProducts.clear();
				
				// Ok, let's process the binding of the adapter proteins (standard adapter)
				String strAdapterXPath = "//receptor_block[@symbol='"+strReceptor+"']/listOfLigands/ligand[@symbol='"+strLigandSymbol+"']/listOfAdapters/adapter/@symbol";
				NodeList nodeListAdapters = (NodeList)_xpath.evaluate(strAdapterXPath,ccmlTree,XPathConstants.NODESET);
				int NUMBER_OF_ADAPTERS = nodeListAdapters.getLength();
				for (int adapter_index=0;adapter_index<NUMBER_OF_ADAPTERS;adapter_index++)
				{
					// Get the ligand -
					Node tmpNodeAdapter = nodeListAdapters.item(adapter_index);
					String strAdapterSymbol = tmpNodeAdapter.getNodeValue();
					
					// Adapter binding
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strReceptorCompartment);
					arrReactants.add(strAdapterSymbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
					encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
					arrReactants.clear();
					arrProducts.clear();
				}
			}
		}	
	}
}
