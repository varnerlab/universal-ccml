package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import javax.xml.xpath.XPathConstants;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MAWNTNetworkHandler extends CCMLMAObject implements IReactionHandler {

	private void init(Document ccmlTree) throws Exception
	{
		// Method attributes -
		String strXPathBase = "";
			
		// Get the global symbols -
		strXPathBase = "//listOfGlobalSymbols/global_symbol";
		populateProperties(strXPathBase,ccmlTree);
		
		// Get the siganling components -
		strXPathBase = "//signaling_block[@block_class='WNT_NETWORK']/listOfSignalingComponents/signaling_component";
		populateProperties(strXPathBase,ccmlTree);
		
		// Get the list of prefixes -
		String strPrefixXPath = "//listOfSymbolPrefixes/symbol_prefix";
		populateProperties(strPrefixXPath,ccmlTree);
	}
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler -
		init(ccmlTree);
		
		// Encode the destruction complex interactions -
		encodeDestructionComplex(arrRxnList,ccmlTree);
		
		// Encode the Ub interactions -
		encodeBCatUbReactions(arrRxnList,ccmlTree);
		
		// Build the interface reactions (use the parent) -
		buildInterfaceReactions("WNT_NETWORK",arrRxnList,ccmlTree);
		
		// Encode degradation reactions -
		buildDegradationReactions("WNT_NETWORK",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	
	// Encode the UB of BCat -
	private void encodeBCatUbReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get the symbols involved with WNT -
		String strPhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strSignalingCompartment = (String)getProperty("SIGNALING_COMPARTMENT");
		String strUBPrefix = (String)getProperty("UBIQUITIN_PREFIX");
		String strBCatSymbol = (String)getProperty("BCATENIN_SYMBOL");
		String strUbLigaseSymbol = (String)getProperty("BCTRCP_SYMBOL");
		
		// Bind BCTRCP -
		arrReactants.add(strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		arrReactants.add(strUbLigaseSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strUbLigaseSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Add the Ub group --
		arrReactants.add(strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strUbLigaseSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strUbLigaseSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strUBPrefix+"_"+strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}
	
	// Encode the destruction complex and phosphorylation of BCat -
	private void encodeDestructionComplex(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();

		// Get the symbols involved with WNT -
		String strPhosphorylationPrefix = (String)getProperty("DOUBLE_PHOSPHORYLATION_PREFIX");
		String strUBPrefix = (String)getProperty("UBIQUITIN_PREFIX");
		String strCKISymbol = (String)getProperty("CK1_SYMBOL");
		String strGSK3BSymbol = (String)getProperty("GSK3B_SYMBOL");
		String strBCatSymbol = (String)getProperty("BCATENIN_SYMBOL");
		String strAPCSymbol = (String)getProperty("APC_SYMBOL");
		String strAXINSymbol = (String)getProperty("AXIN_SYMBOL");
		String strSignalingCompartment = (String)getProperty("SIGNALING_COMPARTMENT");
		
		// Ok - so we need to encode the formation of destruction complex
		// Bind AXIN -
		arrReactants.add(strAXINSymbol+"_"+strSignalingCompartment);
		arrReactants.add(strBCatSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Bind APC -
		arrReactants.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		arrReactants.add(strAPCSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylate Bcat - CK1 binding 
		arrReactants.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strSignalingCompartment);
		arrReactants.add(strCKISymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strCKISymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylate Bcat - GSK3B binding 
		arrReactants.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strSignalingCompartment);
		arrReactants.add(strGSK3BSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strGSK3BSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylate with CK1 -
		arrReactants.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strCKISymbol+"_"+strSignalingCompartment);
		arrProducts.add(strCKISymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAPCSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylate with GSK3B -
		arrReactants.add(strAXINSymbol+"_"+strBCatSymbol+"_"+strAPCSymbol+"_"+strGSK3BSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strGSK3BSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAXINSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strAPCSymbol+"_"+strSignalingCompartment);
		arrProducts.add(strPhosphorylationPrefix+"_"+strBCatSymbol+"_"+strSignalingCompartment);
		encodeMassActionSBMLReaction(arrRxnList, ccmlTree, arrReactants, arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}

}
