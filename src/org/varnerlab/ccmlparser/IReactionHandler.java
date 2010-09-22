package org.varnerlab.ccmlparser;

import java.util.ArrayList;

import org.w3c.dom.Document;

public interface IReactionHandler {

	// Method used for reactions -
	public void constructNetworkReactions(ArrayList<String> arrRxnList,Document ccmlTree) throws Exception;
	public void setProperty(String key,Object value);
	public Object getProperty(String key);
	
	
}
