package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;
import java.util.Hashtable;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MABasalGeneExpressionHandler extends CCMLMAObject implements IReactionHandler {
	
	
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
		String strPrefixXPath = "//listOfSymbolPrefixes/prefix";
		populateProperties(strPrefixXPath,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// initialize -
		init(ccmlTree);
		
		// Ok, get the compartment names etc for gene expression -
		String strCompartment = doCCMLCompartmentLookup(ccmlTree,"NUCLEAR_KEY");
		String strTranslationCompartment = doCCMLCompartmentLookup(ccmlTree,"CYTOSOL_KEY");
		
		
		String strTranslationRibosome = (String)getProperty("RIBOSOME_SYMBOL");
		String strRNAP = (String)getProperty("RNA_POLYMERASE_SYMBOL");
		String strExport = (String)getProperty("NUCLEAR_TRANSPORTER_OUT");
		String strImport = (String)getProperty("NUCLEAR_TRANSPORTER_IN");
		
		// Get the MRNA and G symbols -
		String strMRNASymbol = (String)getProperty("MRNA_PREFIX");
		String strGeneSymbol = (String)getProperty("GENE_PREFIX");
		String strActiveSymbol = (String)getProperty("START_PREFIX");
	
		
		// Grab the symbol from basal gene block -
		String strXPath = "//basal_gene/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,ccmlTree,XPathConstants.NODESET);
		int NUMBER_OF_BASAL_GENES = nodeList.getLength();
		for (int gene_index=0;gene_index<NUMBER_OF_BASAL_GENES;gene_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(gene_index);
			String strRawGeneSymbol = tmpNode.getNodeValue();
			
			// Encode the RNAP binding -
			arrReactants.add(strGeneSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
			arrReactants.add(strRNAP+"_"+strCompartment);
			arrProducts.add(strGeneSymbol+"_"+strRawGeneSymbol+"_"+strRNAP+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Encode mRNA synthesis -
			arrReactants.add(strGeneSymbol+"_"+strRawGeneSymbol+"_"+strRNAP+"_"+strCompartment);
			arrProducts.add(strGeneSymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
			arrProducts.add(strRNAP+"_"+strCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Export of the MRNA from the nucleus (binding w/export)
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strCompartment);
			arrReactants.add(strExport+"_"+strCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strExport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Export the MRNA from the nucleus (transport step)
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strExport+"_"+strCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationCompartment);
			arrProducts.add(strExport+"_"+strCompartment);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Binding of the ribosome -
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationCompartment);
			arrReactants.add(strTranslationRibosome+"_"+strTranslationCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
			encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Scan to active config -
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strActiveSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Generate the protein -
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strActiveSymbol+"_"+strTranslationRibosome+"_"+strTranslationCompartment);
			arrProducts.add(strTranslationRibosome+"_"+strTranslationCompartment);
			arrProducts.add(strRawGeneSymbol+"_"+strTranslationCompartment);
			arrProducts.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationCompartment);
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
			arrReactants.clear();
			arrProducts.clear();
			
			// Degrade the mRNA -
			arrReactants.add(strMRNASymbol+"_"+strRawGeneSymbol+"_"+strTranslationCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();
			
			// Degrade the protein -
			arrReactants.add(strRawGeneSymbol+"_"+strTranslationCompartment);
			arrProducts.add("[]");
			encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.DEGRADATION);
			arrReactants.clear();
			arrProducts.clear();			
		}
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
}
