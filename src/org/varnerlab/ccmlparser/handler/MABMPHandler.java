package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;
import java.util.Hashtable;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.IReceptorNetworkHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MABMPHandler extends CCMLMAObject implements IReceptorNetworkHandler {
	
	
	// Main hook method - this gets called from SBMLCCMLModel -
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
	
		// Method attributes -
		ArrayList<String> arrListLocal = new ArrayList<String>();
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize - get specific names from tree 
		initHandler(ccmlTree);
	
		// Process the network blocks -
		buildReceptorComplex(arrRxnList,ccmlTree);
		buildDegradationReactions(arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);
	}
	
	
	

	// Initialize the handler -
	private void initHandler(Document ccmlTree) throws Exception
	{
		
		// OK, get the name of the ligand -
		String strLigandCXPath = "//receptor_block[@block_class='BMP_SMAD']/listOfLigands/ligand/@compartment";
		String strLigandC = queryCCMLTree(ccmlTree,strLigandCXPath);
		setProperty("LIGAND_COMPARTMENT",strLigandC);
		
		// Get the compartment of the receptor -
		String strReceptorCompartmentXPath = "//receptor_block[@block_class='BMP_SMAD']/@compartment";
		String strReceptorCompartment = queryCCMLTree(ccmlTree,strReceptorCompartmentXPath);
		setProperty("RECEPTOR_COMPARTMENT",strReceptorCompartment);
		
		// Ok, we need to process the list of regulators and store the info in properties -
		String strXPath = "//Receptor_signaling_network_block/listOfReceptors/receptor_block[@block_class='BMP_SMAD']/listOfComponents/component/@key";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_REGULATORS = nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_REGULATORS;index++)
		{
			// Get the key -
			Node tmpNode = nodeList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = "//Receptor_signaling_network_block/listOfReceptors/receptor_block[@block_class='BMP_SMAD']/listOfComponents/component[@key='"+strKeyName+"']/@symbol";
			String strSymbol = queryCCMLTree(ccmlTree,strSymbolXPath);
			
			// store in prop -
			setProperty(strKeyName,strSymbol);
		}
		
		// Process keynames in the adapter_complex tags -
		String strACXPath = "//adapter_complex/@key";
		NodeList nodeACList = (NodeList)_xpath.evaluate(strACXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_ACS= nodeACList.getLength();
		for (int index = 0;index<NUMBER_OF_ACS;index++)
		{
			// Get the key -
			Node tmpNode = nodeACList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = "//adapter_complex[@key='"+strKeyName+"']/@symbol";
			String strSymbol = queryCCMLTree(ccmlTree,strSymbolXPath);
			
			// store in prop -
			setProperty(strKeyName,strSymbol);
		}
	}
	
	private void buildReceptorComplex(ArrayList<String> arrListReactions,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the symbol for the receptor -
		String strLigandCompartment = (String)getProperty("LIGAND_COMPARTMENT");
		String strReceptorCompartment = (String)getProperty("RECEPTOR_COMPARTMENT");
		
		// Ok, get the list of receptors -
		String strReceptorXPath = "//receptor_block[@block_class='BMP_SMAD']/@symbol";
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
				encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
				encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
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
					
					// Adapter binding - (ALKs)
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strReceptorCompartment);
					arrReactants.add(strAdapterSymbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
					arrReactants.clear();
					arrProducts.clear();
					
					/*
					// Bind with SMAD1 -
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					arrReactants.add(strSMAD1Symbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strSMAD1Symbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
					arrReactants.clear();
					arrProducts.clear();
					
					// Convert to PSMAD1
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strSMAD1Symbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					arrProducts.add("P"+strSMAD1Symbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
					arrReactants.clear();
					arrProducts.clear();
					
					// Bind with SMAD5 -
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					arrReactants.add(strSMAD5Symbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strSMAD5Symbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
					arrReactants.clear();
					arrProducts.clear();
					
					// Convert to PSMAD5 -
					arrReactants.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strSMAD5Symbol+"_"+strReceptorCompartment);
					arrProducts.add(strLigandSymbol+"_"+strReceptor+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
					arrProducts.add("P"+strSMAD5Symbol+"_"+strReceptorCompartment);
					encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
					arrReactants.clear();
					arrProducts.clear();
					*/
				}
				
				// Ok, we need to process the complex list of adapters -
				String strAdapterComplexXPath = "//receptor_block[@symbol='"+strReceptor+"']/listOfLigands/ligand[@symbol='"+strLigandSymbol+"']/listOfAdapters/adapter_complex/@symbol";
				NodeList nodeListAdapterComplex = (NodeList)_xpath.evaluate(strAdapterComplexXPath,ccmlTree,XPathConstants.NODESET);
				int NUMBER_OF_ADAPTERS_COMPLEX = nodeListAdapterComplex.getLength();
				for (int complex_index=0;complex_index<NUMBER_OF_ADAPTERS_COMPLEX;complex_index++)
				{
					// Ok, get the next adapter -
					Node tmpNodeAdapter = nodeListAdapterComplex.item(complex_index);
					String strAdapterComplexSymbol = tmpNodeAdapter.getNodeValue();
					
					// Ok, so let's pull the individual adapter symbols -
					String strAdapterComplexSymbolXPath = "//receptor_block[@symbol='"+strReceptor+"']/listOfLigands/ligand[@symbol='"+strLigandSymbol+"']/listOfAdapters/adapter_complex[@symbol='"+strAdapterComplexSymbol+"']/adapter/@symbol";
					NodeList nodeListAdapterComplexSymbols = (NodeList)_xpath.evaluate(strAdapterComplexSymbolXPath,ccmlTree,XPathConstants.NODESET);
					int NUMBER_OF_ADAPTERS_COMPLEX_SYMBOLS = nodeListAdapterComplexSymbols.getLength();
					String strPreviousSymbol = strLigandSymbol+"_"+strReceptor;
					for (int symbol_index=0;symbol_index<NUMBER_OF_ADAPTERS_COMPLEX_SYMBOLS;symbol_index++)
					{
						// Ok, get the next adapter -
						Node tmpNodeAdapterSymbol = nodeListAdapterComplexSymbols.item(symbol_index);
						String strAdapterSymbol = tmpNodeAdapterSymbol.getNodeValue();
						
						// Does this adapter have a key?
						String strKeyXPath = "//receptor_block[@symbol='"+strReceptor+"']/listOfLigands/ligand[@symbol='"+strLigandSymbol+"']/listOfAdapters/adapter_complex/adapter[@symbol='"+strAdapterSymbol+"']/@key";
						String strKeyName = queryCCMLTree(ccmlTree,strKeyXPath);
						
						// Ok - if we don't have a keyname then just add the adapter to the growing complex -
						if (strKeyName.isEmpty())
						{
							// Encode the interaction -
							arrReactants.add(strPreviousSymbol+"_"+strReceptorCompartment);
							arrReactants.add(strAdapterSymbol+"_"+strReceptorCompartment);
							arrProducts.add(strPreviousSymbol+"_"+strAdapterSymbol+"_"+strReceptorCompartment);
							encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
							encodeMassActionSBMLReaction(arrListReactions,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
							arrReactants.clear();
							arrProducts.clear();
							
							// Update the previous symbol -
							strPreviousSymbol = strPreviousSymbol+"_"+strAdapterSymbol;
						}
						else if (strKeyName.equalsIgnoreCase("SHC_SYMBOL"))
						{
							// Do nothing for now ...
						}
					}
				}
			}
		}
	}
	
	private void buildDegradationReactions(ArrayList<String> arrListReactions,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//receptor_block[@block_class='TGFB_SMAD']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//receptor_block[@block_class='BMP_SMAD']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(doc,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrListReactions,doc,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}


}
