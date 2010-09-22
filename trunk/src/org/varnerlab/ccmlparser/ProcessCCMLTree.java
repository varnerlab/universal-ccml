package org.varnerlab.ccmlparser;

import java.io.File;

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
	
	
	public void processDocumentTree(Document doc) throws Exception {
		// Method attributes -
		SBMLCCMLModel model_wrapper = new SBMLCCMLModel();
		StringBuffer speciesBuffer = new StringBuffer();
		StringBuffer reactionBuffer = new StringBuffer();
		StringBuffer parameterBuffer = new StringBuffer();
		
		// Ok, get the name of the *array* file -
		String strArrayFileName = "/JobConfiguration/input_file_name/@name";
		Node fileNode = (Node)_xpath.evaluate(strArrayFileName, doc, XPathConstants.NODE);
		String strBlueprintFile = fileNode.getNodeValue();
		
		// Get the filename that I'm going to write to -
		String strOutputFileNameXP = "/JobConfiguration/output_file_name/@name";
		Node outNode = (Node)_xpath.evaluate(strOutputFileNameXP, doc, XPathConstants.NODE);
		String strOutputFileName = outNode.getNodeValue();
		
		// Load the blueprint xml file -
		Document ccmlTree = loadBlueprintFile(strBlueprintFile);
		
		// Put in the sbml header -
		model_wrapper.populateHeader(_sbmlBuffer, ccmlTree);
		
		// Build list of compartments -
		model_wrapper.populateListOfCompartments(_sbmlBuffer, ccmlTree);
		
		// Build list of species -
		model_wrapper.populateReactionList(reactionBuffer, ccmlTree,doc);
		model_wrapper.populateParameterBuffer(parameterBuffer, ccmlTree,doc);
		model_wrapper.populateSpeciesList(speciesBuffer,ccmlTree,doc);
		
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
