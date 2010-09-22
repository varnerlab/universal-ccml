package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.varnerlab.ccmlparser.ReactionType;
import org.w3c.dom.Document;

public class MANFKBNetworkHandler extends CCMLMAObject implements IReactionHandler {

	
	
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Use the CCMLObject's init method -
		init("NFKB_NETWORK",ccmlTree);
		
		// Ok, so we need to encode a few reactions dealing with the assembly and dissociation of IKB_NKKB
		encodeNFKBReactions(arrRxnList,ccmlTree);
		
		// Process the interface block (use CCML's method)
		buildInterfaceReactions("NFKB_NETWORK",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}
	
	private void encodeNFKBReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Get component symbols -
		String strCompartment = (String)getProperty("SIGNALING_COMPARTMENT");
		String strIKBSymbol = (String)getProperty("IKB_SYMBOL");
		String strNFKBSymbol = (String)getProperty("NFKB_SYMBOL");
		String strPhosphorylation = (String)getProperty("PHOSPHORYLATION_PREFIX");
		
		// Ok, make sure these guys are bound -
		arrReactants.add(strIKBSymbol+"_"+strCompartment);
		arrReactants.add(strNFKBSymbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		arrReactants.clear();
		arrProducts.clear();
		
		// Only allow the phosphorylated variant to come off -
		arrReactants.add(strPhosphorylation+"_"+strIKBSymbol+"_"+strNFKBSymbol+"_"+strCompartment);
		arrProducts.add(strPhosphorylation+"_"+strIKBSymbol+"_"+strCompartment);
		arrProducts.add(strNFKBSymbol+"_"+strCompartment);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrReactants,arrProducts,ReactionType.FORWARD_RATE);
		encodeMassActionSBMLReaction(arrRxnList,ccmlTree,arrProducts,arrReactants,ReactionType.REVERSE_RATE);
		arrReactants.clear();
		arrProducts.clear();
	}

}
