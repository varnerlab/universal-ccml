package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReceptorNetworkHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MAReceptorHandler extends CCMLMAObject implements IReceptorNetworkHandler {
	
	@Override
	public void constructNetworkReactions(ArrayList<String> arrRxnList,
			Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler (use the super object --)
		initHandler(ccmlTree);
		
		// Process the interface block (use the inherited method)
		buildInterfaceReactions("DEFAULT_RECEPTOR",arrRxnList,ccmlTree);
		
		// Process the degradation block (use CCML parent's method)
		buildDegradationReactions("DEFAULT_RECEPTOR",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
		
	}

	private void initHandler(Document ccmlTree) throws Exception
	{
		
		// OK, get the name of the ligand -
		String strLigandCXPath = "//receptor_block[@block_class='DEFAULT_RECEPTOR']/listOfLigands/ligand/@compartment";
		String strLigandC = queryCCMLTree(ccmlTree,strLigandCXPath);
		setProperty("LIGAND_COMPARTMENT",strLigandC);
		
		// Get the compartment of the receptor -
		String strReceptorCompartmentXPath = "//receptor_block[@block_class='DEFAULT_RECEPTOR']/@compartment";
		String strReceptorCompartment = queryCCMLTree(ccmlTree,strReceptorCompartmentXPath);
		setProperty("RECEPTOR_COMPARTMENT",strReceptorCompartment);
		
		// Ok, we need to process the list of regulators and store the info in properties -
		String strXPath = "//Receptor_network_block/listOfReceptors/receptor_block[@block_class='DEFAULT_RECEPTOR']/listOfComponents/component/@key";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_REGULATORS = nodeList.getLength();
		for (int index = 0;index<NUMBER_OF_REGULATORS;index++)
		{
			// Get the key -
			Node tmpNode = nodeList.item(index);
			String strKeyName = tmpNode.getNodeValue();
			
			// Get the symbol -
			String strSymbolXPath = "//Receptor_network_block/listOfReceptors/receptor_block[@block_class='DEFAULT_RECEPTOR']/listOfComponents/component[@key='"+strKeyName+"']/@symbol";
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
	
	protected void buildDegradationReactions(String strBlockName,ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the list of receptor species that degrade -
		String strXPath = "//receptor_block[@block_class='"+strBlockName+"']/listOfDegradation/degrade/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// Get the compartment that this species is in -
			String strXPathCompartment = "//receptor_block[@block_class='"+strBlockName+"']/listOfDegradation/degrade[@symbol='"+strSpecies+"']/@compartment";
			String strCompartment = queryCCMLTree(ccmlTree,strXPathCompartment);
			
			// Encode the degrdation reaction -
			arrReactants.add(strSpecies+"_"+strCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
		}
	}
}
