package org.varnerlab.ccmlparser;

import java.io.File;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

public class CCMLParser {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Method attributes -
        String strPropPath = args[0];
        
        try 
        {
        	// Load the XML prop file -
        	File configFile = new File(strPropPath);
        	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        	dbFactory.setNamespaceAware(true);
        	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
      	  	Document doc = dBuilder.parse(configFile);
      	  	doc.getDocumentElement().normalize();
        	
      	  	// Ok, so we need to get the name of the logic handler - get the class name of the handler -
      	  	XPathFactory  _xpFactory = XPathFactory.newInstance();
      	  	XPath _xpath = _xpFactory.newXPath();
      	  	
      	  	String strXPath = ".//CCMLParserConfiguration/classname/@name";
      	  	Node classNode = (Node)_xpath.evaluate(strXPath,doc,XPathConstants.NODE);
      	  	
      	  	// Create an instance of the logic handler -
      	  	ILogicHandler handler = (ILogicHandler)Class.forName(classNode.getNodeValue()).newInstance();
      	  	
        	// Ok, so we have loaded the file. Nows lets read it and dump the file to disk
      	  	handler.processDocumentTree(doc);	
        }
        catch (Exception error)
        {
        	System.out.println("ERROR: Hey now, the message server did not start -"+error.toString());
        	error.printStackTrace();
        }
	}
	
	
}
