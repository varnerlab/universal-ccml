package org.varnerlab.ccmlparser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

public class ProcessCCMLTree implements ILogicHandler {

	// Class/instance attributes -
	private XPathFactory  _xpFactory = XPathFactory.newInstance();
	private XPath _xpath = _xpFactory.newXPath();
	private StringBuffer _sbmlBuffer = new StringBuffer();
	
	private Document loadBlueprintFile(String fileName) throws Exception {
		
		// Load the XML prop file -
    	File configFile = new File(fileName);
    	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
    	dbFactory.setNamespaceAware(true);
    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
  	  	Document doc = dBuilder.parse(configFile);
  	  	doc.getDocumentElement().normalize();
  	  	
  	  	// Return the document -
  	  	return(doc);
	}
	
	private ArrayList<String> loadParameterArrayFile(String strFileName) throws Exception
	{
		// Method attributes -
		ArrayList<String> parameterList = new ArrayList<String>();
		
		// Ok, so we need to load the parameter file, parse the contents -
		FileReader parameterFile = new FileReader(strFileName);
		BufferedReader inReader=new BufferedReader(parameterFile);
		String dataRecord=null;
		while ((dataRecord=inReader.readLine())!=null)
		{
			// Ok, add to the list -
			parameterList.add(dataRecord);
		}
		
		
		// return the list of parameters - 
		return(parameterList);
	}
	
	
	public void processDocumentTree(Document doc) throws Exception {
		// Method attributes -
		SBMLCCMLModel model_wrapper = new SBMLCCMLModel();
		StringBuffer speciesBuffer = new StringBuffer();
		StringBuffer reactionBuffer = new StringBuffer();
		StringBuffer parameterBuffer = new StringBuffer();
		
		ArrayList<String> currentReactionList = new ArrayList<String>();
		ArrayList<String> previousReactionList = new ArrayList<String>();
		ArrayList<String> newReactionList = new ArrayList<String>();
		
		Document ccmlTree = null;
		Document ccmlTreePrevious = null;
		
		// Ok, get the name of the input CCML file -
		String strArrayFileName = "/JobConfiguration/input_file_name/@name";
		Node fileNode = (Node)_xpath.evaluate(strArrayFileName, doc, XPathConstants.NODE);
		String strBlueprintFile = fileNode.getNodeValue();
		
		// Get the filename that I'm going to write to -
		String strOutputFileNameXP = "/JobConfiguration/output_file_name/@name";
		Node outNode = (Node)_xpath.evaluate(strOutputFileNameXP, doc, XPathConstants.NODE);
		String strOutputFileName = outNode.getNodeValue();
		
		// Get the name of the reference CCML (previous CCML version that we'll use to order the reactions)
		String strPrevArrayFileNameXP = "/JobConfiguration/previous_input_file_name/@name";
		Node prevFileNode = (Node)_xpath.evaluate(strPrevArrayFileNameXP, doc, XPathConstants.NODE);
		String strPrevBlueprintFile = prevFileNode.getNodeValue();
		
		// Get the name of the parameter input file (if there is one...)
		String strParaArrayFileNameXP = "/JobConfiguration/parameter_input_file_name/@name";
		Node paraFileNode = (Node)_xpath.evaluate(strParaArrayFileNameXP, doc, XPathConstants.NODE);
		String strParaFile = paraFileNode.getNodeValue();
		
		// Load the blueprint xml file -
		ccmlTree = loadBlueprintFile(strBlueprintFile);
		
		// Ok, we need to check to see if I have a previous blueprint file, if so load it
		if (!strPrevBlueprintFile.isEmpty())
		{
			ccmlTreePrevious = loadBlueprintFile(strPrevBlueprintFile);
		}
		
		
		// Put in the sbml header -
		model_wrapper.populateHeader(_sbmlBuffer, ccmlTree);
		
		// Build list of compartments -
		model_wrapper.populateListOfCompartments(_sbmlBuffer, ccmlTree);
		
		// Build list of reactions -
		if (ccmlTreePrevious!=null)
		{
			// Ok, so I get here, then I have a previous CCML. I'm going to use this previous file to establish the order of the reactions.
			// Any new reactions get appended to the *bottom* of the reaction list
			model_wrapper.populateReactionList(currentReactionList, ccmlTree,doc);
			model_wrapper.populateReactionList(previousReactionList, ccmlTreePrevious,doc);
			
			// Ok, so I should have both reactions lists - now what??
			if (currentReactionList.size()>previousReactionList.size())
			{
				// Ok, if I get here I have *new* reactions that I need to append to the bottom of the reaction list -
				int NUMBER_PREVIOUS_REACTIONS = previousReactionList.size();
				for (int previous_reaction_index=0;previous_reaction_index<NUMBER_PREVIOUS_REACTIONS;previous_reaction_index++)
				{
					// Get a reaction from the previous list -
					String strPreviousReaction = previousReactionList.get(previous_reaction_index);
					
					// Check to see if this reaction is contained in the new list -
					// If it is, remove from the currentList and add to the new list -
					if (currentReactionList.contains(strPreviousReaction))
					{
						// Ok, we *have* the reaction in the current list - remove and add to the new list
						currentReactionList.remove(strPreviousReaction);
						
						// Add to the new list -
						newReactionList.add(strPreviousReaction);
					}
				}
				
				// Ok, so when I get here I *should* only have reactions in the currentList that are *different* than the old list.
				// I need to add the *different* reactions to the newReactionList -
				int NUMBER_DIFFERENT_REACTIONS = currentReactionList.size();
				for (int different_reaction_index=0;different_reaction_index<NUMBER_DIFFERENT_REACTIONS;different_reaction_index++)
				{
					newReactionList.add(currentReactionList.get(different_reaction_index));
				}
				
				// Ok, so now I have the new reaction list - use it to populate parameters, species etc.
				model_wrapper.convertReactionListToStringBuffer(newReactionList, reactionBuffer);
				model_wrapper.populateSpeciesList(speciesBuffer,newReactionList,ccmlTree,doc);
				
				if (strParaFile.isEmpty())
				{
					// No previous parameter need to be loaded -- generate new random parameters -
					model_wrapper.populateParameterBufferFromReactions(0,parameterBuffer,newReactionList,ccmlTree, doc);
				}
				else
				{
					// Ok, so we have some *old* parameters that need to be loaded. Load these up.
					ArrayList<String> oldParameterList = new ArrayList<String>();
					oldParameterList = loadParameterArrayFile(strParaFile);
					
					// Ok, so now we need to populate the parameter buffer *including* the old parameters -
					model_wrapper.populateParameterBufferFromReactions(0,parameterBuffer,newReactionList,oldParameterList,ccmlTree, doc);
				}
			}
			else
			{
				// Ok, for some reason the current reaction list is *shorter* than the old one. Just process normally...
				model_wrapper.convertReactionListToStringBuffer(currentReactionList, reactionBuffer);
				model_wrapper.populateSpeciesList(speciesBuffer,currentReactionList,ccmlTree,doc);
				model_wrapper.populateParameterBufferFromReactions(0,parameterBuffer,currentReactionList,ccmlTree, doc);
			}
		}
		else
		{
			// Populate the species and parameter list -
			model_wrapper.convertReactionListToStringBuffer(currentReactionList, reactionBuffer);
			model_wrapper.populateSpeciesList(speciesBuffer,currentReactionList,ccmlTree,doc);
			model_wrapper.populateParameterBufferFromReactions(0,parameterBuffer,currentReactionList,ccmlTree, doc);
		}
		
		
		// Add the species to the sbmlBuffer -
		_sbmlBuffer.append("\n");
		_sbmlBuffer.append("\t\t<listOfSpecies>\n");
		_sbmlBuffer.append(speciesBuffer);
		_sbmlBuffer.append("\t\t</listOfSpecies>\n");
		_sbmlBuffer.append("\n");
		
		// Add the parameters to the sbmlBuffer -
		_sbmlBuffer.append("\n");
		_sbmlBuffer.append("\t\t<listOfParameters>\n");
		_sbmlBuffer.append(parameterBuffer);
		_sbmlBuffer.append("\t\t</listOfParameters>\n");
		_sbmlBuffer.append("\n");
		
		// Add the reactions to the sbmlBuffer -
		_sbmlBuffer.append("\n");
		_sbmlBuffer.append("\t\t<listOfReactions>\n");
		_sbmlBuffer.append(reactionBuffer);
		_sbmlBuffer.append("\t\t</listOfReactions>\n");
		_sbmlBuffer.append("\n");
		
		
		// Put in the sbml footer -
		model_wrapper.populateFooter(_sbmlBuffer, ccmlTree);
		
		// dump the buffer to disk -
		GIOL.write(strOutputFileName,_sbmlBuffer);
	}
	
	// Generate 
}
