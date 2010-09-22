package org.varnerlab.ccmlparser.handler;

import java.util.ArrayList;

import org.varnerlab.ccmlparser.CCMLMAObject;
import org.varnerlab.ccmlparser.IReceptorNetworkHandler;
import org.w3c.dom.Document;

public class MAWNTReceptorHandler extends CCMLMAObject implements IReceptorNetworkHandler {

	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception {
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		
		// Ok, now that we have done all the reactions, we need to determine the species list -
		arrSpecies = generateSpeciesList(arrRxnList,ccmlTree);
		
		// Ok, add the types of reactions to the properties -
		setProperty("REACTION_TYPE_LIST",_arrListReactionType);
		setProperty("RECEPTOR_SPECIES_LIST",arrSpecies);	
		
	}

}
