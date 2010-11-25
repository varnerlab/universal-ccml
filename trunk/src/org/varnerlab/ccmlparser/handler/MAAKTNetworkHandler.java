package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;

public class MAAKTNetworkHandler extends CCMLMAObject implements IReactionHandler {

	
	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler -
		init("AKT_NETWORK",ccmlTree);
		
		// Ok, build the reactions that encode AKT activation -
		encodeAKTActivation(arrRxnList,ccmlTree);
		
		// Process the interface block (use the inhereted method)
		buildInterfaceReactions("AKT_NETWORK",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
		
	}
	
	
	// Encode AKT activation -
	private void encodeAKTActivation(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get component strings -
		// String strCompartment = (String)getProperty("SIGNALING_COMPARTMENT");
		String strCompartment = doCCMLCompartmentLookup(ccmlTree,"CYTOSOL_KEY");
		
		String strActivePI3KSymbol = (String)getProperty("ACTIVE_PI3K_SYMBOL");
		String strActivated = (String)getProperty("ACTIVATED_PREFIX");
		String strAKTComponent = (String)getProperty("AKT_SYMBOL");
		String strPIP2Component = (String)getProperty("PIP2_SYMBOL");
		String strPIP3Component = (String)getProperty("PIP3_SYMBOL");
		String strPDK1Component = (String)getProperty("PDK1_SYMBOL");
		String strPTENComponent = (String)getProperty("PTEN_SYMBOL");
		
		// Encode the conversion of PIP2 to PIP3 -
		// Binding step -
		arrReactants.add(strActivePI3KSymbol+"_"+strCompartment);
		arrReactants.add(strPIP2Component+"_"+strCompartment);
		arrProducts.add(strActivePI3KSymbol+"_"+strPIP2Component+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Conversion step -
		arrReactants.add(strActivePI3KSymbol+"_"+strPIP2Component+"_"+strCompartment);
		arrProducts.add(strActivePI3KSymbol+"_"+strCompartment);
		arrProducts.add(strPIP3Component+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Back reaction from PTEN -
		// Binding step -
		arrReactants.add(strPTENComponent+"_"+strCompartment);
		arrReactants.add(strPIP3Component+"_"+strCompartment);
		arrProducts.add(strPTENComponent+"_"+strPIP3Component+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Conversion step -
		arrReactants.add(strPTENComponent+"_"+strPIP3Component+"_"+strCompartment);
		arrProducts.add(strPTENComponent+"_"+strCompartment);
		arrProducts.add(strPIP2Component+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Binding of PIP3 and PKD1 -
		// Binding step -
		arrReactants.add(strPDK1Component+"_"+strCompartment);
		arrReactants.add(strPIP3Component+"_"+strCompartment);
		arrProducts.add(strPDK1Component+"_"+strPIP3Component+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Activation of AKT by PDK1 -
		// Binding -
		arrReactants.add(strPDK1Component+"_"+strPIP3Component+"_"+strCompartment);
		arrReactants.add(strAKTComponent+"_"+strCompartment);
		arrProducts.add(strPDK1Component+"_"+strPIP3Component+"_"+strAKTComponent+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Phosphorylation of AKT -
		arrReactants.add(strPDK1Component+"_"+strPIP3Component+"_"+strAKTComponent+"_"+strCompartment);
		arrProducts.add(strPDK1Component+"_"+strPIP3Component+"_"+strCompartment);
		arrProducts.add(strActivated+"_"+strAKTComponent+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.CATALYTIC_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}

}
