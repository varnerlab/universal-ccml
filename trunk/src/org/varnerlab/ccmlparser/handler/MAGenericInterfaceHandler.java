package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReactionHandler;
import org.w3c.dom.Document;

public class MAGenericInterfaceHandler extends CCMLMAObject implements IReactionHandler  {

	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Initialize the handler (use the super object --)
		init("DEFAULT_SIGNALING_NETWORK",ccmlTree);
		
		// Process the interface block (use the inhereted method)
		buildInterfaceReactions("DEFAULT_SIGNALING_NETWORK",arrRxnList,ccmlTree);
		
		// Process the degrdation block (use CCML parent's method)
		buildDegradationReactions("DEFAULT_SIGNALING_NETWORK",arrRxnList,ccmlTree);
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
	}

}
