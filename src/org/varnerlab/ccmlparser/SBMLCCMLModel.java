package org.varnerlab.ccmlparser;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.varnerlab.ccmlparser.handler.MAJAKStatHandler;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public class SBMLCCMLModel extends CCMLMAObject {
	// Class/instance attributes -
	private ArrayList<String> _arrListSBMLSpecies = new ArrayList<String>();
	private int _intReactionCounter = 0;
	//private ArrayList<ReactionType> _arrListReactionType = new ArrayList<ReactionType>();
	
	
	public void populateHeader(StringBuffer buffer,Document doc) throws Exception
	{
		buffer.append("<?xml version=\"1.0\"?>\n");
		buffer.append("<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">\n");
		
		// Grab the name of the model -
		String strModelNameXPath = "/Model/@name";
		String strModelName = queryCCMLTree(doc,strModelNameXPath);
		
		buffer.append("\t<model id=\"");
		buffer.append(strModelName);
		buffer.append("\">\n");
	}
	
	public void populateFooter(StringBuffer buffer,Document doc) throws Exception
	{
		buffer.append("\t</model>\n");
		buffer.append("</sbml>");
	}
	
	public void populateParameterBufferFromReactions(int offset,StringBuffer buffer,ArrayList<String> rxnList,ArrayList<String> oldParametersList,Document ccmlTree,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> newReactionList = new ArrayList<String>();
		
		// Copy the rxnList -
		newReactionList.addAll(rxnList);
		
		// Get a list of types -
		ArrayList<String> arrListTypes = generateReactionTypeList(rxnList,ccmlTree);
		
		// Ok, first we need to put the old parameter in the first whatever rows -
		int NUMBER_OLD_PARAMETERS = oldParametersList.size();
		for (int old_parameter_index=0;old_parameter_index<NUMBER_OLD_PARAMETERS;old_parameter_index++)
		{
			
			// Get the parameter value -
			String strOldParameterValue = oldParametersList.get(old_parameter_index);
				
			// Construct the parameter line -
			buffer.append("\t\t\t<parameter id=\"PARAMETER_R_");
			buffer.append(old_parameter_index+offset);
			buffer.append("\" name=\"k_");
			buffer.append(old_parameter_index+offset);
			buffer.append("\" value=\"");
			buffer.append(strOldParameterValue);
			buffer.append("\"/>\n");
			
			// Remove the elements from the newReactionList -
			newReactionList.remove(old_parameter_index);
		}
		
		// Ok, so know we need to process the rest of the list -
		populateParameterBufferFromReactions(NUMBER_OLD_PARAMETERS,buffer,newReactionList,ccmlTree,doc);
	}
	
	public void populateParameterBufferFromReactions(int offset,StringBuffer buffer,ArrayList<String> rxnList,Document ccmlTree,Document doc) throws Exception
	{
		// Method attributes -
		double dblValue = 0.0;
		
		// Get a list of types -
		ArrayList<String> arrListTypes = generateReactionTypeList(rxnList,ccmlTree);
		
		// Ok, grab the list of reaction types and then generate random parameters on the proper scale
		int NUMBER_OF_REACTIONS = arrListTypes.size();
		for (int parameter_index=0;parameter_index<NUMBER_OF_REACTIONS;parameter_index++)
		{
			// Get the reaction type -
			String strRxnTmp = arrListTypes.get(parameter_index);
			
			if (strRxnTmp.equalsIgnoreCase("REACTION_FORWARD"))
			{
				dblValue = 1*Math.random();
			}
			else if (strRxnTmp.equalsIgnoreCase("REACTION_REVERSE"))
			{
				dblValue = 0.1*Math.random();
			}
			else if (strRxnTmp.equalsIgnoreCase("REACTION_CATALYTIC"))
			{
				dblValue = 1*Math.random();
			}
			else if (strRxnTmp.equalsIgnoreCase("REACTION_DEGRADATION"))
			{
				dblValue = 0.5*Math.random();
			}
			else if (strRxnTmp.equalsIgnoreCase("REACTION_GENERATION"))
			{
				dblValue = 1*Math.random();
			}
			
			// Construct the parameter line -
			buffer.append("\t\t\t<parameter id=\"PARAMETER_R_");
			buffer.append(parameter_index+offset);
			buffer.append("\" name=\"k_");
			buffer.append(parameter_index+offset);
			buffer.append("\" value=\"");
			buffer.append(dblValue);
			buffer.append("\"/>\n");
		}
		
	}
	
	/*
	public void populateParameterBuffer(StringBuffer buffer,Document ccmlTree,Document doc) throws Exception
	{
		// Method attributes -
		double dblValue = 0.0;
		
		// Ok, grab the list of reaction types and then generate random parameters on the proper scale
		int NUMBER_OF_REACTIONS = _arrListReactionType.size();
		for (int parameter_index=0;parameter_index<NUMBER_OF_REACTIONS;parameter_index++)
		{
			// Get the reaction type -
			ReactionType rxnType = _arrListReactionType.get(parameter_index);
			switch (rxnType)
			{
				case FORWARD_RATE:
					dblValue = 1*Math.random();
					break;
				case REVERSE_RATE:
					dblValue = 0.1*Math.random();
					break;
				case CATALYTIC_RATE:
					dblValue = 1*Math.random();
					break;
				case DEGRADATION:
					dblValue = Math.random();
					break;
				case GENERATION:
					dblValue = Math.random();
					break;
			}
			
			// Construct the parameter line -
			buffer.append("\t\t\t<parameter id=\"PARAMETER_R_");
			buffer.append(parameter_index);
			buffer.append("\" name=\"k_");
			buffer.append(parameter_index);
			buffer.append("\" value=\"");
			buffer.append(dblValue);
			buffer.append("\"/>\n");
		}
	}*/
	
	public void populateReactionList(ArrayList<String> reactionArray, Document ccmlTree,Document doc) throws Exception
	{
		// Method attributes -
		StringBuffer buffer = new StringBuffer();
		
		// Just to make sure all is ok - let's clear out the reaction list -
		_arrListSBMLReactions.clear();
		
		// Put comments in the SBML so we can see which block the reactions come from -
		processMembraneTransportBlock(buffer,ccmlTree,doc);
		processBasalExpressionReactionBlock(buffer, ccmlTree,doc);
		processRegulatedExpressionReactionBlock(buffer, ccmlTree,doc);
		processInfrastructureSynthesisReactionsBlock(buffer, ccmlTree);
		processReceptorNetworkBlock(buffer, ccmlTree, doc);
		processSignalingNetworkBlock(buffer, ccmlTree,doc);
		
		// Ok, so we need to add these to the buffer -
		int NUMBER_OF_REACTIONS = _arrListSBMLReactions.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = _arrListSBMLReactions.get(index);
			//final_buffer.append(tmpReaction);
			//final_buffer.append("\n");
			reactionArray.add(tmpReaction);
		}
	}
	
	public void convertReactionListToStringBuffer(ArrayList<String> reactionArray,StringBuffer final_buffer) throws Exception
	{
		int NUMBER_OF_REACTIONS = reactionArray.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = reactionArray.get(index);
			
			final_buffer.append("\t\t<!-- REACTION ");
			final_buffer.append(index);
			final_buffer.append(" -->\n");
			final_buffer.append(tmpReaction);
			final_buffer.append("\n");
		}
	}
	
	
	private void processMembraneTransportBlock(StringBuffer buffer,Document ccmlTree,Document configTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrLBlockClassList = new ArrayList<String>();
		ArrayList<String> rxnList = _arrListSBMLReactions;
		
		// Ok, so I need to get the list of block classes -
		String strBlockClassXPath = "//Membrane_transport_block/@block_class";
		String strBlockClassName = queryCCMLTree(ccmlTree, strBlockClassXPath);
		
		// Ok, so we need to load the matching java classname -
		String strClassNameXPath = "//mapping[@keyname='"+strBlockClassName+"']/@classname";
		String strClassName = queryCCMLTree(configTree,strClassNameXPath);
		
		
		// Initiate the handler -
		IReactionHandler handler = (IReactionHandler)Class.forName(strClassName).newInstance();
		
		// Ok, for now hard-code the handler (for testing)
		//ArrayList<String> rxnList = new ArrayList<String>();
		//ArrayList<String> rxnList = new ArrayList<String>();
		
		// Build the reactions -
		handler.constructNetworkReactions(rxnList,ccmlTree);
		
		/*
		// Ok, so we need to add these to the buffer -
		int NUMBER_OF_REACTIONS = rxnList.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = rxnList.get(index);
			buffer.append(tmpReaction);
		}
		*/
		
		// Grab the list of reaction types and add to the main list -
		ArrayList<ReactionType> localTypeList = (ArrayList<ReactionType>)handler.getProperty("REACTION_TYPE_LIST");
		_arrListReactionType.addAll(localTypeList);
		
		// Grab the list of species list -
		ArrayList<String> localSpeciesList = (ArrayList<String>)handler.getProperty("RECEPTOR_SPECIES_LIST");
		int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
		for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
		{
			// Get the species -
			String strLocalSpecies = localSpeciesList.get(species_index);
			
			// Check for the null species - do not include in the list of species 
			if (!strLocalSpecies.equalsIgnoreCase("[]"))
			{
				// Add it to the list -
				addSBMLSpeciesToList(strLocalSpecies,"0.0");
			}
		}
	}
	
	private void processSignalingNetworkBlock(StringBuffer buffer,Document ccmlTree,Document configTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrLBlockClassList = new ArrayList<String>();
		ArrayList<String> rxnList = _arrListSBMLReactions;
		
		// Ok, so I need to get the list of block classes -
		String strBlockClassXPath = "//signaling_block/@block_class";
		NodeList nodeList = (NodeList)_xpath.evaluate(strBlockClassXPath,ccmlTree,XPathConstants.NODESET);
		
		// Go through the list of the block_class names and get a unique list -
		int NUMBER_NAMES_TMP = nodeList.getLength();
		for (int tmp_index=0;tmp_index<NUMBER_NAMES_TMP;tmp_index++)
		{
			// Ok, get the block class -
			Node tmpNode = nodeList.item(tmp_index);
			String strBlockClassName = tmpNode.getNodeValue();
			
			// Check to see if this is already in the list -
			if (!arrLBlockClassList.contains(strBlockClassName))
			{
				// Ok, we do *not* have this block name - store it
				arrLBlockClassList.add(strBlockClassName);
			}
		}
		
		
		int NUMBER_OF_SIGNAL_BLOCKS = arrLBlockClassList.size();
		for (int signal_index=0;signal_index<NUMBER_OF_SIGNAL_BLOCKS;signal_index++)
		{
			
			// Get the class name -
			String strBlockClassName = arrLBlockClassList.get(signal_index);
			
			// Get the java classname -
			String strClassNameXPath = "//mapping[@keyname='"+strBlockClassName+"']/@classname";
			String strClassName = queryCCMLTree(configTree,strClassNameXPath);
			
			// Initiate the handler -
			IReactionHandler handler = (IReactionHandler)Class.forName(strClassName).newInstance();
			
			// Ok, for now hard-code the handler (for testing)
			//ArrayList<String> rxnList = new ArrayList<String>();
			
			// Build the reactions -
			handler.constructNetworkReactions(rxnList,ccmlTree);
			
			/*
			// Ok, so we need to add these to the buffer -
			int NUMBER_OF_REACTIONS = rxnList.size();
			for (int index=0;index<NUMBER_OF_REACTIONS;index++)
			{
				String tmpReaction = rxnList.get(index);
				buffer.append(tmpReaction);
			}*/
			
			// Grab the list of reaction types and add to the main list -
			ArrayList<ReactionType> localTypeList = (ArrayList<ReactionType>)handler.getProperty("REACTION_TYPE_LIST");
			_arrListReactionType.addAll(localTypeList);
			
			// Grab the list of species list -
			ArrayList<String> localSpeciesList = (ArrayList<String>)handler.getProperty("RECEPTOR_SPECIES_LIST");
			int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
			for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
			{
				// Get the species -
				String strLocalSpecies = localSpeciesList.get(species_index);
				
				// Check for the null species - do not include in the list of species 
				if (!strLocalSpecies.equalsIgnoreCase("[]"))
				{
					// Add it to the list -
					addSBMLSpeciesToList(strLocalSpecies,"0.0");
				}
			}
		}
	}
	
	private void processRegulatedExpressionReactionBlock(StringBuffer buffer,Document ccmlTree,Document configTree) throws Exception
	{
		// Method attributes -
		
		// Ok, so I need to get the list of block classes -
		String strBlockClassXPath = "//Regulated_expression_block/@block_class";
		String strBlockClassName = queryCCMLTree(ccmlTree, strBlockClassXPath);
		
		// Ok, so we need to load the matching java classname -
		String strClassNameXPath = "//mapping[@keyname='"+strBlockClassName+"']/@classname";
		String strClassName = queryCCMLTree(configTree,strClassNameXPath);
		
		
		// Initiate the handler -
		IReactionHandler handler = (IReactionHandler)Class.forName(strClassName).newInstance();
		
		// Ok, for now hard-code the handler (for testing)
		//ArrayList<String> rxnList = new ArrayList<String>();
		ArrayList<String> rxnList = _arrListSBMLReactions;
		
		// Build the reactions -
		handler.constructNetworkReactions(rxnList,ccmlTree);
		
		/*
		// Ok, so we need to add these to the buffer -
		int NUMBER_OF_REACTIONS = rxnList.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = rxnList.get(index);
			buffer.append(tmpReaction);
		}*/
		
		// Grab the list of reaction types and add to the main list -
		ArrayList<ReactionType> localTypeList = (ArrayList<ReactionType>)handler.getProperty("REACTION_TYPE_LIST");
		_arrListReactionType.addAll(localTypeList);
		
		// Grab the list of species list -
		ArrayList<String> localSpeciesList = (ArrayList<String>)handler.getProperty("RECEPTOR_SPECIES_LIST");
		int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
		for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
		{
			// Get the species -
			String strLocalSpecies = localSpeciesList.get(species_index);
			
			// Check for the null species - do not include in the list of species 
			if (!strLocalSpecies.equalsIgnoreCase("[]"))
			{
				// Add it to the list -
				addSBMLSpeciesToList(strLocalSpecies,"0.0");
			}
		}
	}
	
	// Add the synthesis and degrdation reactions for the infrastruture species -
	private void processInfrastructureSynthesisReactionsBlock(StringBuffer buffer,Document doc) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrReactants = new ArrayList<String>();
		ArrayList<String> arrProducts = new ArrayList<String>();
		
		// Execute XPath call to get the list of infrastructure species -
		String strInfrastructureXPath = "//Infrastructure_synthesis_block/infrastructure/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strInfrastructureXPath,doc,XPathConstants.NODESET);
		int NUMBER_INFRASTRUCTURE_COMPONENTS = nodeList.getLength();
		for (int component_index=0;component_index<NUMBER_INFRASTRUCTURE_COMPONENTS;component_index++)
		{
			// Ok, get the symbol node -
			Node tmpNode = nodeList.item(component_index);
			String strComponentSymbol = tmpNode.getNodeValue();
			
			// Ok, grab the compartment -
			String strComponentXPath = "//Infrastructure_synthesis_block/infrastructure[@symbol='"+strComponentSymbol+"']/@compartment_key";
			String strCompartmentKey = queryCCMLTree(doc,strComponentXPath);
			String strCompartment = doCCMLCompartmentLookup(doc,strCompartmentKey);
			
			// Encode the reactions -
			arrReactants.add(strComponentSymbol+"_"+strCompartment);
			arrProducts.add("[]");
			generateIrreversibleGeneralReaction(buffer,doc,arrReactants,arrProducts,ReactionType.DEGRADATION);
			generateIrreversibleGeneralReaction(buffer,doc,arrProducts,arrReactants,ReactionType.GENERATION);
			
			// Clear out the lists -
			arrReactants.clear();
			arrProducts.clear();
		}
		
	}
	
	// Process the Network receptor block -
	private void processReceptorNetworkBlock(StringBuffer buffer, Document ccmlTree,Document configTree) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrLBlockClassList = new ArrayList<String>();
		ArrayList<String> rxnList = _arrListSBMLReactions;
		
		// Ok, so I need to get the list of block classes -
		String strBlockClassXPath = "//receptor_block/@block_class";
		NodeList nodeList = (NodeList)_xpath.evaluate(strBlockClassXPath,ccmlTree,XPathConstants.NODESET);
		
		// Go through the list of the block_class names and get a unique list -
		int NUMBER_NAMES_TMP = nodeList.getLength();
		for (int tmp_index=0;tmp_index<NUMBER_NAMES_TMP;tmp_index++)
		{
			// Ok, get the block class -
			Node tmpNode = nodeList.item(tmp_index);
			String strBlockClassName = tmpNode.getNodeValue();
			
			// Check to see if this is already in the list -
			if (!arrLBlockClassList.contains(strBlockClassName))
			{
				// Ok, we do *not* have this block name - store it
				arrLBlockClassList.add(strBlockClassName);
			}
		}
		
		
		
		// Ok, so I need to get the list of block classes -
		int NUMBER_OF_BLOCKS = arrLBlockClassList.size();
		for (int block_index=0;block_index<NUMBER_OF_BLOCKS;block_index++)
		{
			// Get the class name -
			String strBlockClassName = arrLBlockClassList.get(block_index);
			
			// Ok, so we need to load the matching java classname -
			String strClassNameXPath = "//mapping[@keyname='"+strBlockClassName+"']/@classname";
			String strClassName = queryCCMLTree(configTree,strClassNameXPath);
			
			// Initiate the handler -
			IReceptorNetworkHandler handler = (IReceptorNetworkHandler)Class.forName(strClassName).newInstance();
			
			// Ok, for now hard-code the handler (for testing)
			//ArrayList<String> rxnList = new ArrayList<String>();
			
			// Build the reactions -
			handler.constructNetworkReactions(rxnList, ccmlTree);
			
			/*
			// Ok, so we need to add these to the buffer -
			int NUMBER_OF_REACTIONS = rxnList.size();
			for (int index=0;index<NUMBER_OF_REACTIONS;index++)
			{
				String tmpReaction = rxnList.get(index);
				buffer.append(tmpReaction);
			}*/
			
			// Grab the list of reaction types and add to the main list -
			ArrayList<ReactionType> localTypeList = (ArrayList<ReactionType>)handler.getProperty("REACTION_TYPE_LIST");
			_arrListReactionType.addAll(localTypeList);
			
			// Grab the list of species list -
			ArrayList<String> localSpeciesList = (ArrayList<String>)handler.getProperty("RECEPTOR_SPECIES_LIST");
			int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
			for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
			{
				// Get the species -
				String strLocalSpecies = localSpeciesList.get(species_index);
				
				// Check for the null species - do not include in the list of species 
				if (!strLocalSpecies.equalsIgnoreCase("[]"))
				{
					// Add it to the list -
					addSBMLSpeciesToList(strLocalSpecies,"0.0");
				}
			}
		}
	}
	
	
	public void populateSpeciesList(StringBuffer buffer,ArrayList<String> rxnList,Document ccmlTree,Document doc) throws Exception
	{
		
		// Generate the list of species from the reaction list -
		ArrayList<String> localSpeciesList = generateSpeciesList(rxnList);
		
		int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
		for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
		{
			// Get the species -
			String strLocalSpecies = localSpeciesList.get(species_index);
			
			// Check for the null species -- do not include in the list of species 
			
			if (!strLocalSpecies.equalsIgnoreCase("[]"))
			{
				// Add it to the list -
				addSBMLSpeciesToList(strLocalSpecies,"0.0");
			}
		}
		
		int NUMBER_OF_SPECIES = _arrListSBMLSpecies.size();
		for (int species_index = 0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			String strTmp = _arrListSBMLSpecies.get(species_index);
			buffer.append(strTmp);
		}
	}
	
	
	// Determine the reactions from basal expression -
	private void processBasalExpressionReactionBlock(StringBuffer buffer,Document doc,Document configTree) throws Exception
	{
		// Method attributes -
		
		// Ok, so I need to get the list of block classes -
		String strBlockClassXPath = "//Basal_expression_block/@block_class";
		String strBlockClassName = queryCCMLTree(doc, strBlockClassXPath);
		
		// Ok, so we need to load the matching java classname -
		String strClassNameXPath = "//mapping[@keyname='"+strBlockClassName+"']/@classname";
		String strClassName = queryCCMLTree(configTree,strClassNameXPath);
		
		
		// Initiate the handler -
		IReactionHandler handler = (IReactionHandler)Class.forName(strClassName).newInstance();
		
		// Ok, for now hard-code the handler (for testing)
		//ArrayList<String> rxnList = new ArrayList<String>();
		ArrayList<String> rxnList = _arrListSBMLReactions;
		
		// Build the reactions -
		handler.constructNetworkReactions(rxnList,doc);
		
		// Ok, so we need to add these to the buffer -
		int NUMBER_OF_REACTIONS = rxnList.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = rxnList.get(index);
			buffer.append(tmpReaction);
		}
		
		// Grab the list of reaction types and add to the main list -
		ArrayList<ReactionType> localTypeList = (ArrayList<ReactionType>)handler.getProperty("REACTION_TYPE_LIST");
		_arrListReactionType.addAll(localTypeList);
		
		// Grab the list of species list -
		ArrayList<String> localSpeciesList = (ArrayList<String>)handler.getProperty("RECEPTOR_SPECIES_LIST");
		int NUMBER_OF_LOCAL_SPECIES = localSpeciesList.size();
		for (int species_index=0;species_index<NUMBER_OF_LOCAL_SPECIES;species_index++)
		{
			// Get the species -
			String strLocalSpecies = localSpeciesList.get(species_index);
			
			// Check for the null species - do not include in the list of species 
			if (!strLocalSpecies.equalsIgnoreCase("[]"))
			{
				// Add it to the list -
				addSBMLSpeciesToList(strLocalSpecies,"0.0");
			}
		}
		
		// put the reaction list in properties 
		this.setProperty("LOCAL_REACTION_LIST",rxnList);
	}
	
	private void generateIrreversibleGeneralReaction(StringBuffer buffer,Document doc,ArrayList<String> arrReactants,ArrayList<String> arrProducts,ReactionType rxnType) throws Exception
	{
		// Method attributes -
		StringBuffer local_buffer = new StringBuffer();
		
		
		// ok - build the MRNA generation reaction (forward) -
		local_buffer.append("\t\t<reaction id=\"REACTION");
		
		switch (rxnType)
		{
			case FORWARD_RATE:
				local_buffer.append("_FORWARD");
				break;
			case REVERSE_RATE:
				local_buffer.append("_REVERSE");
				break;
			case CATALYTIC_RATE:
				local_buffer.append("_CATALYTIC");
				break;
			case DEGRADATION:
				local_buffer.append("_DEGRADATION");
				break;
			case GENERATION:
				local_buffer.append("_GENERATION");
				break;
		}
		
		local_buffer.append("\" name=\"");
		
		int NUMBER_OF_REACTANTS = arrReactants.size();
		for (int reactant_index=0;reactant_index<NUMBER_OF_REACTANTS;reactant_index++)
		{
			String strTmpSymbol = arrReactants.get(reactant_index);
			
			// Add the symbol to the buffer -
			local_buffer.append(strTmpSymbol);
			
			if (reactant_index!=NUMBER_OF_REACTANTS-1)
			{
				local_buffer.append(" + ");
			}
		}
		
		local_buffer.append(" = ");
		
		// Get the list of products -
		int NUMBER_OF_PRODUCTS = arrProducts.size();
		for (int product_index=0;product_index<NUMBER_OF_PRODUCTS;product_index++)
		{
			String strTmpSymbolRaw = arrProducts.get(product_index);
			local_buffer.append(strTmpSymbolRaw);
			
			if (product_index!=NUMBER_OF_PRODUCTS-1)
			{
				local_buffer.append(" + ");
			}
		}
		
		local_buffer.append("\" reversible=\"false\">\n");
		local_buffer.append("\t\t\t<listOfReactants>\n");
		
		for (int reactant_index=0;reactant_index<NUMBER_OF_REACTANTS;reactant_index++)
		{
			String strTmpSymbol = arrReactants.get(reactant_index);
			local_buffer.append("\t\t\t\t<speciesReference species=\"");
			local_buffer.append(strTmpSymbol);
			local_buffer.append("\"/>\n");	
		}
		
		local_buffer.append("\t\t\t</listOfReactants>\n");
		local_buffer.append("\t\t\t<listOfProducts>\n");
		
		for (int product_index=0;product_index<NUMBER_OF_PRODUCTS;product_index++)
		{
			String strTmpSymbolRaw = arrProducts.get(product_index);			
			local_buffer.append("\t\t\t\t<speciesReference species=\"");
			local_buffer.append(strTmpSymbolRaw);
			local_buffer.append("\"/>\n");
			
		}
		
		local_buffer.append("\t\t\t</listOfProducts>\n");
		local_buffer.append("\t\t</reaction>\n");
		
		// Reaction string -
		String strTmpString = local_buffer.toString();
		
		// Check to see if this string is in the main buffer already -
		if (!_arrListSBMLReactions.contains(strTmpString))
		{
			// Ok, add the reaction string to the buffer -
			_arrListSBMLReactions.add(strTmpString);
			
			// Add the local buffer to the main buffer -
			buffer.append(local_buffer);
			
			// Ok, if we keep this reaction then I need to store its type -
			_arrListReactionType.add(rxnType);
			
			// update the counter - 
			_intReactionCounter++;
		}	
	}
	
	private ArrayList<String> generateSpeciesList(ArrayList<String> arrListReactions) throws Exception
	{
		// Method attributes -
		ArrayList<String> arrSpecies = new ArrayList<String>();
		StringBuffer local_buffer = new StringBuffer();
		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
    	dbFactory.setNamespaceAware(true);
    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
		
		// Ok, so we are going to turn the reaction list into an xml tree and then use xpath to extract all the species -
		local_buffer.append("<?xml version=\"1.0\"?>\n");
		local_buffer.append("<listOfReactions>\n");
		int NUMBER_OF_REACTIONS = arrListReactions.size();
		for (int index=0;index<NUMBER_OF_REACTIONS;index++)
		{
			String tmpReaction = arrListReactions.get(index);
			local_buffer.append(tmpReaction);
		}
		local_buffer.append("</listOfReactions>\n");
		String strTmp = local_buffer.toString();
		
		// Ok, so now we need to create a document -
		Document reaction_dom_tree = dBuilder.parse(new InputSource(new StringReader(local_buffer.toString())));
		
		// Create the XPath -
		String strXPath = "//speciesReference/@species";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,reaction_dom_tree,XPathConstants.NODESET);
		int NUMBER_OF_SPECIES = nodeList.getLength();
		for (int species_index=0;species_index<NUMBER_OF_SPECIES;species_index++)
		{
			// Get the gene symbol -
			Node tmpNode = nodeList.item(species_index);
			String strSpecies = tmpNode.getNodeValue();
			
			// check to see if we already have the species in the list -
			if (!arrSpecies.contains(strSpecies))
			{
				arrSpecies.add(strSpecies);
			}
		}
		
		// return the array -
		return(arrSpecies);
	}
	
	// Determine the network compartments -
	public void populateListOfCompartments(StringBuffer buffer,Document doc) throws Exception
	{
		// Method attributes -
		
		// Populate --
		buffer.append("\t\t<listOfCompartments>\n");
	
		// XPath string -
		String strXPath = "//compartment/@symbol";
		NodeList cList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int INT_COMPARTMENT_LIST = cList.getLength();
		for (int index=0;index<INT_COMPARTMENT_LIST;index++)
		{
			// Get node -
			Node tmpNode = cList.item(index);
			
			// Get the compartment -
			String strCompartmentID = tmpNode.getNodeValue();
			
			// populate the buffer -
			buffer.append("\t\t\t<compartment id=\"");
			buffer.append(strCompartmentID);
			buffer.append("\" name=\"");
			buffer.append(strCompartmentID);
			buffer.append("\" />\n");
		}
		
		// list close -
		buffer.append("\t\t</listOfCompartments>\n");
	}
	
	private void addSBMLSpeciesToList(String strGeneSymbol,String strInitialCondition) throws Exception
	{
		// Method attributes -
		StringBuffer local_buffer = new StringBuffer(); 
		
		// Ok, so this is sort of a HACK .. but I need to get the compartment some how - 
		// we have the convention that the name ends in _<COMPARTMENT>. Grad that name -
		int last_slash = strGeneSymbol.lastIndexOf("_");
		String strCompartment = strGeneSymbol.substring(last_slash+1, strGeneSymbol.length());
		
		
		local_buffer.append("\t\t\t<species id=\"");
		local_buffer.append(strGeneSymbol);
		local_buffer.append("\" name=\"");
		local_buffer.append(strGeneSymbol);
		local_buffer.append("\" compartment=\"");
		local_buffer.append(strCompartment);
		local_buffer.append("\" initialAmount=\"");
		local_buffer.append(strInitialCondition);
		local_buffer.append("\"/>\n");
		
		String strBufferEntry = local_buffer.toString();
		
		// Ok, so I have a local_buffer with an sbml species in it - have we done this one before?
		if (!_arrListSBMLSpecies.contains(strBufferEntry))
		{
			// Ok, so we have *not* seen this species before -
			// Add the sbml line to the list -
			_arrListSBMLSpecies.add(strBufferEntry);
		}
	}
	
	
	
}
