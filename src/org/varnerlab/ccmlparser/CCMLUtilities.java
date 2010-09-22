package org.varnerlab.ccmlparser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Vector;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;


public class CCMLUtilities {

	
	public static void processNetworkImportFiles(StringBuffer buffer,Document doc,Document jobdoc) throws Exception 
	{
		// Method attributes -
		XPathFactory  _xpFactory = XPathFactory.newInstance();
		XPath _xpath = _xpFactory.newXPath();

		
		// Put a comment into -
		buffer.append("// ------- IMPORTED TEMPLATE BLOCKS FROM LIBRARY (CHECK NOMENCLATURE...) ------------ //\n");
		
		
		String strLibXPath = "//library/@path_prefix";
		Node tmpNodePrefix = (Node)_xpath.evaluate(strLibXPath, jobdoc,XPathConstants.NODE);
		String strLibPathPrefix = tmpNodePrefix.getNodeValue();

		String strNameXPath = "//SignalingNetworks/import/@name";
		NodeList tmpNodeList = (NodeList)_xpath.evaluate(strNameXPath, doc,XPathConstants.NODESET);
		int NUMBER_OF_IMPORTS = tmpNodeList.getLength();
		for (int index=0;index<NUMBER_OF_IMPORTS;index++)
		{
			// Get the name -
			String strName = tmpNodeList.item(index).getNodeValue(); 

			// Setup to get filename 
			String strFileXPath="//SignalingNetworks//import[@name='"+strName+"']/@file";
			Node tmpNode = (Node)_xpath.evaluate(strFileXPath,doc,XPathConstants.NODE);
			String strFileName = tmpNode.getNodeValue();

			// Put a comment -
			buffer.append("// ");
			buffer.append(strName);
			buffer.append("\n");

			// Load the import file -
			String strTemplateFile = strLibPathPrefix+strFileName;
			File iFile=new File(strTemplateFile);
			BufferedReader reader=new BufferedReader(new FileReader(iFile));
			String s="";
			while ((s=reader.readLine())!=null)
			{
				buffer.append(s);
				buffer.append("\n");
			}

			// put a newline -
			buffer.append("\n");
		}
		
		// end comment line -
		buffer.append("// ---------------------------------------------------------------------------------- // \n");
		buffer.append("\n");
	}
	
	
	public static void processSecretionBlock(StringBuffer buffer, Document doc) throws Exception
	{
		// Method attributes -
		XPathFactory  _xpFactory = XPathFactory.newInstance();
		XPath _xpath = _xpFactory.newXPath();
		
		buffer.append("// ---------------------------- TRANSLOCATION --------------------------------------- // \n");
		
		// Formulate and execute xpath string -
		String strXPath = "//export_block/export/@symbol";
		NodeList exportList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_EXPORTS = exportList.getLength();
		for (int export_index=0;export_index<NUMBER_OF_EXPORTS;export_index++)
		{
			// Ok - get the symbol of the gene we are looking at -
			Node tmpNode = exportList.item(export_index);
			String strSymbol = tmpNode.getNodeValue();
			
			// populate the buffer -
			buffer.append("EXPORT_BINDING_");
			buffer.append(strSymbol);
			buffer.append(",");
			buffer.append(strSymbol);
			buffer.append("+EXTEXPORT,");
			buffer.append(strSymbol);
			buffer.append("_EXTEXPORT,-inf,inf;\n");
			buffer.append("EXPORT_TRANSPORT_");
			buffer.append(strSymbol);
			buffer.append(",");
			buffer.append(strSymbol);
			buffer.append("_EXTEXPORT,");
			buffer.append(strSymbol);
			buffer.append("xt+EXTEXPORT,0,inf;\n");
			
		}
		
		buffer.append("// ---------------------------------------------------------------------------------- // \n");
		buffer.append("\n");
	}
	
	public static void processRegulatedExpressionBlock(StringBuffer buffer,Document doc) throws Exception
	{
		// Method attributes -
		XPathFactory  _xpFactory = XPathFactory.newInstance();
		XPath _xpath = _xpFactory.newXPath();
		
		// Ok, so we'll need to process a series of blocks -
		String strXPath = "//regulated_gene/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_REGULATED_GENES = nodeList.getLength();
		for (int gene_index=0;gene_index<NUMBER_OF_REGULATED_GENES;gene_index++)
		{
			// Ok - get the symbol of the gene we are looking at -
			Node tmpNode = nodeList.item(gene_index);
			String strGeneSybmol = tmpNode.getNodeValue();
			
			buffer.append("// ----- REGULATED TRANSCRIPTION/TRANSLATION ");
			buffer.append(strGeneSybmol);
			buffer.append(" ------------------------------------- // \n");

			
			// From here, I need to run a xquery to get the list of activators -
			String strListOfActivatorsXPath = "//regulated_gene[@symbol='"+strGeneSybmol+"']/listOfActivators/activator/@symbol";
			NodeList listActivators = (NodeList)_xpath.evaluate(strListOfActivatorsXPath,doc,XPathConstants.NODESET);
			
			// Ok, so how many activators do we have?
			int NUMBER_OF_ACTIVATORS = listActivators.getLength();
			for (int activator_index=0;activator_index<NUMBER_OF_ACTIVATORS;activator_index++)
			{
				Node tmpANode = listActivators.item(activator_index);
				String strASybmol = tmpANode.getNodeValue();
				
				// populate a single TF -
				populateSingleRegulatedTranscriptionBuffer(buffer,strGeneSybmol,strASybmol);
				
			}
			
			// Ok, so we need to process the activator complex block -
			String strActivatorComplexNameXPath = "//regulated_gene[@symbol='"+strGeneSybmol+"']/listOfActivators/activator_complex/@symbol";
			NodeList listAComplexNames = (NodeList)_xpath.evaluate(strActivatorComplexNameXPath,doc,XPathConstants.NODESET);
			int NUMBER_ACTIVATOR_COMPLEXES = listAComplexNames.getLength();
			for (int activator_complex_block_index=0;activator_complex_block_index<NUMBER_ACTIVATOR_COMPLEXES;activator_complex_block_index++)
			{
			
				Node tmpACCNode = listAComplexNames.item(activator_complex_block_index);
				String strACCSybmol = tmpACCNode.getNodeValue();
				
				// Ok, when I get here I have an activator complex - we know the name. Need to grab the components of this complex -
				String strACComponentsXPath = "//regulated_gene[@symbol='"+strGeneSybmol+"']/listOfActivators/activator_complex[@symbol='"+strACCSybmol+"']/activator/@symbol";
				
				// Get the list of activator components -
				NodeList listActivatorComplexComponents = (NodeList)_xpath.evaluate(strACComponentsXPath,doc,XPathConstants.NODESET);
				int NUMBER_OF_COMPONENTS = listActivatorComplexComponents.getLength();
				ArrayList<String> tmpList = new ArrayList<String>();
				for (int tmp_index=0;tmp_index<NUMBER_OF_COMPONENTS;tmp_index++)
				{
					// Ok, so we need to get the name of the components -
					Node tmpComponentNode = listActivatorComplexComponents.item(tmp_index);
					String strComponentSymbol = tmpComponentNode.getNodeValue();
					
					// Add the symbol to the list -
					tmpList.add(strComponentSymbol);
				}
				
				// Ok, so we should have the list of the components and the final name - write the block 
				processNuclearComplexAssembly(buffer,tmpList,strACCSybmol);
			}
			
			
			// I need to run a xpath query to process the repressors -
			String strListOfRepressorsXPath = "//regulated_gene[@symbol='"+strGeneSybmol+"']/listOfRepressors/repressor/@symbol";
			NodeList listRepressors = (NodeList)_xpath.evaluate(strListOfRepressorsXPath,doc,XPathConstants.NODESET);
			
			// Ok, so how many repressord do we have?
			int NUMBER_OF_REPRESSORS = listRepressors.getLength();
			for (int repressor_index=0;repressor_index<NUMBER_OF_REPRESSORS;repressor_index++)
			{
				// Get the node -
				Node tmpRNode = listRepressors.item(repressor_index);
				String strRSybmol = tmpRNode.getNodeValue();
				
				// Process the repressor -
				populateSingleRegulatedTranscriptionBufferRepressor(buffer,strGeneSybmol,strRSybmol);
			}
			
			
			buffer.append("// ---------------------------------------------------------------------------------- // \n");
			buffer.append("\n");
		}
	}
	
	private static void processNuclearComplexAssembly(StringBuffer buffer, ArrayList<String> tmpList, String strACCSybmol) {
		// Ok, so we need to process this list -
		buffer.append("\n");
		buffer.append("// Assemble the ");
		buffer.append(strACCSybmol);
		buffer.append(" complex \n");
		
		// Get how many steps -
		
		buffer.append("// ASSEMBLE_BINDING_");
		buffer.append(strACCSybmol);
		buffer.append("\n");
		
		
	}

	private static void populateSingleRegulatedTranscriptionBufferRepressor(StringBuffer buffer, String strGeneSymbol, String strRSybmol) {
		
		// Repressor block -
		buffer.append("\n");
		buffer.append("// Repression of  ");
		buffer.append(strGeneSymbol);
		buffer.append(" expression by ");
		buffer.append(strRSybmol);
		buffer.append("\n");
		buffer.append("REPRESSOR_BINDING_");
		buffer.append(strRSybmol);
		buffer.append(",g_");
		buffer.append(strGeneSymbol);
		buffer.append("+");
		buffer.append(strRSybmol);
		buffer.append("_N,g_");
		buffer.append(strGeneSymbol);
		buffer.append("_");
		buffer.append(strRSybmol);
		buffer.append("_N,-inf,inf;\n");
		
		// Add the transport line for the repressor -
		buffer.append("\n");
		buffer.append("// Transport into/from the nucleus for ");
		buffer.append(strRSybmol);
		buffer.append("\n");
		buffer.append("NCULEAR_TRANSPORT_BINDING_");
		buffer.append(strRSybmol);
		buffer.append(",");
		buffer.append(strRSybmol);
		buffer.append("+NIMPORT,");
		buffer.append(strRSybmol);
		buffer.append("_NIMPORT,-inf,inf;\n");
		
		buffer.append("NCULEAR_TRANSPORT_");
		buffer.append(strRSybmol);
		buffer.append(",");
		buffer.append(strRSybmol);
		buffer.append("_NIMPORT,");
		buffer.append(strRSybmol);
		buffer.append("_N+NIMPORT,0,inf;\n");
		
		// Transport of mRNA out of the nucleus -
		buffer.append("EXONCULEAR_TRANSPORT_BINDING_");
		buffer.append(strRSybmol);
		buffer.append(",");
		buffer.append(strRSybmol);
		buffer.append("_N+NEXPORT,");
		buffer.append(strRSybmol);
		buffer.append("_N_NEXPORT,-inf,inf;\n");
		
		buffer.append("EXONCULEAR_TRANSPORT_");
		buffer.append(strRSybmol);
		buffer.append(",");
		buffer.append(strRSybmol);
		buffer.append("_N_NEXPORT,");
		buffer.append(strRSybmol);
		buffer.append("+NEXPORT,0,inf;\n");
		
	}

	public static void processInfrastructureBlock(StringBuffer buffer, Document doc) throws Exception
	{
		// Method attributes -
		XPathFactory  _xpFactory = XPathFactory.newInstance();
		XPath _xpath = _xpFactory.newXPath();
		
		// Formulate the xpath string -
		String strXPath = "//infrastructure/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_INFRASTRUCTURE_COMPONENTS = nodeList.getLength();
		
		// add a comment line -
		buffer.append("// ----------- INFRASTRUCTURE BLOCK ------------------------------------------------- // \n");
		
		for (int index=0;index<NUMBER_OF_INFRASTRUCTURE_COMPONENTS;index++)
		{
			Node tmpNode = nodeList.item(index);
			String strGeneSybmol = tmpNode.getNodeValue();
			
			buffer.append("SYNTHESIS_");
			buffer.append(strGeneSybmol);
			buffer.append(",[],");
			buffer.append(strGeneSybmol);
			buffer.append(",-inf,inf;\n");
		}
		
		buffer.append("// ---------------------------------------------------------------------------------- // \n");
		buffer.append("\n");
	}
	
	public static void processBasalBlock(StringBuffer buffer,Document doc) throws Exception
	{
		// Method attributes -
		XPathFactory  _xpFactory = XPathFactory.newInstance();
		XPath _xpath = _xpFactory.newXPath();
		ArrayList<String> localVector = new ArrayList<String>();
		
		// Formulate the xpath string -
		String strXPath = "//basal_gene/@symbol";
		NodeList nodeList = (NodeList)_xpath.evaluate(strXPath,doc,XPathConstants.NODESET);
		int NUMBER_OF_BASAL_GENES = nodeList.getLength();
		for (int gene_index=0;gene_index<NUMBER_OF_BASAL_GENES;gene_index++)
		{
			// Get the if of the this receptor -
			Node tmpNode = nodeList.item(gene_index);
			String strGeneSybmol = tmpNode.getNodeValue();
			
			buffer.append("// ----- BASAL TRANSCRIPTION/TRANSLATION ");
			buffer.append(strGeneSybmol);
			buffer.append(" ---------------------------------------- // \n");
			
			// Ok, bitches .. so now I will process the gene symbols -
			populateBasalTranscriptionBuffer(buffer,strGeneSybmol);
			populateTranslationBuffer(buffer,strGeneSybmol);
			
			buffer.append("// ------------------------------------------------------------------------------------ // \n");
			buffer.append("\n");
		}
		
	}
	
	// Add reactions to buffer for regulated expression of gene x with a single transcription factor y-
	private static void populateSingleRegulatedTranscriptionBuffer(StringBuffer buffer,String strGSymbol,String strTFSymbol) throws Exception
	{
		
		
		// RNAP binding step -
		buffer.append("\n");
		buffer.append("// Regulated expression for gene g_");
		buffer.append(strGSymbol);
		buffer.append("\n");
		
		buffer.append("TF_BINDING_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("g_");
		buffer.append(strGSymbol);
		buffer.append("+");
		buffer.append(strTFSymbol);
		buffer.append("_N,g_");
		buffer.append(strGSymbol);
		buffer.append("_");
		buffer.append(strTFSymbol);
		buffer.append("_N,-inf,inf;\n");
		
		buffer.append("RNAP_BINDING_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("g_");
		buffer.append(strGSymbol);
		buffer.append("_");
		buffer.append(strTFSymbol);
		buffer.append("_N+RNAP,g_");
		buffer.append(strGSymbol);
		buffer.append("_");
		buffer.append(strTFSymbol);
		buffer.append("_N_RNAP,-inf,inf;\n");
		
		// Generation of mRNA -
		buffer.append("TRANSCRIPTION_");
		buffer.append(strGSymbol);
		buffer.append(",g_");
		buffer.append(strGSymbol);
		buffer.append("_");
		buffer.append(strTFSymbol);
		buffer.append("_N_RNAP,mRNA_");
		buffer.append(strGSymbol);
		buffer.append("_N+RNAP+g_");
		buffer.append(strGSymbol);
		buffer.append("+");
		buffer.append(strTFSymbol);
		buffer.append("_N,0,inf;\n");
		
		// Transport of mRNA out of the nucleus -
		buffer.append("\n");
		buffer.append("// Transport into/from the nucleus for ");
		buffer.append(strTFSymbol);
		buffer.append("\n");
		buffer.append("NCULEAR_TRANSPORT_BINDING_");
		buffer.append(strTFSymbol);
		buffer.append(",");
		buffer.append(strTFSymbol);
		buffer.append("+NIMPORT,");
		buffer.append(strTFSymbol);
		buffer.append("_NIMPORT,-inf,inf;\n");
		
		buffer.append("NCULEAR_TRANSPORT_");
		buffer.append(strTFSymbol);
		buffer.append(",");
		buffer.append(strTFSymbol);
		buffer.append("_NIMPORT,");
		buffer.append(strTFSymbol);
		buffer.append("_N+NIMPORT,0,inf;\n");
		
		// Transport of mRNA out of the nucleus -
		buffer.append("EXONCULEAR_TRANSPORT_BINDING_");
		buffer.append(strTFSymbol);
		buffer.append(",");
		buffer.append(strTFSymbol);
		buffer.append("_N+NEXPORT,");
		buffer.append(strTFSymbol);
		buffer.append("_N_NEXPORT,-inf,inf;\n");
		
		buffer.append("EXONCULEAR_TRANSPORT_");
		buffer.append(strTFSymbol);
		buffer.append(",");
		buffer.append(strTFSymbol);
		buffer.append("_N_NEXPORT,");
		buffer.append(strTFSymbol);
		buffer.append("+NEXPORT,0,inf;\n");
	}
	
	// Add reactions to buffer for basal expression of gene "symbol"
	private static void populateBasalTranscriptionBuffer(StringBuffer buffer,String strGSymbol) throws Exception
	{
		// RNAP binding step -
		buffer.append("// Basal expression for gene g_");
		buffer.append(strGSymbol);
		buffer.append("\n");
		
		buffer.append("RNAP_BINDING_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("g_");
		buffer.append(strGSymbol);
		buffer.append("+RNAP,g_");
		buffer.append(strGSymbol);
		buffer.append("_RNAP,-inf,inf;\n");
		
		// Generation of mRNA -
		buffer.append("TRANSCRIPTION_");
		buffer.append(strGSymbol);
		buffer.append(",g_");
		buffer.append(strGSymbol);
		buffer.append("_RNAP,mRNA_");
		buffer.append(strGSymbol);
		buffer.append("_N+RNAP+g_");
		buffer.append(strGSymbol);
		buffer.append(",0,inf;\n");
		
		// Transport of mRNA out of the nucleus -
		buffer.append("EXONCULEAR_TRANSPORT_BINDING_");
		buffer.append(strGSymbol);
		buffer.append(",mRNA_");
		buffer.append(strGSymbol);
		buffer.append("_N+NEXPORT,mRNA_");
		buffer.append(strGSymbol);
		buffer.append("_N_NEXPORT,-inf,inf;\n");
		
		buffer.append("EXONCULEAR_TRANSPORT_");
		buffer.append(strGSymbol);
		buffer.append(",mRNA_");
		buffer.append(strGSymbol);
		buffer.append("_N_NEXPORT,mRNA_");
		buffer.append(strGSymbol);
		buffer.append("+NEXPORT,0,inf;\n");
		
		// Add a new line -
		buffer.append("\n");
	}
	
	// Add reactions to the buffer for basal translation of mRNA "symbol"
	private static void populateTranslationBuffer(StringBuffer buffer,String strGSymbol) throws Exception
	{
		// Ribosome binding -
		buffer.append("// Translation of gene product mRNA_");
		buffer.append(strGSymbol);
		buffer.append("\n");
		buffer.append("RIBOSOME_BINDING_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("mRNA_");
		buffer.append(strGSymbol);
		buffer.append("+RIBOSOME,RM_");
		buffer.append(strGSymbol);
		buffer.append(",-inf,inf;\n");
		
		// Scanning step -
		buffer.append("RIBOSOME_SCANNING_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("RM_");
		buffer.append(strGSymbol);
		buffer.append(",AR_");
		buffer.append(strGSymbol);
		buffer.append(",0,inf;\n");
		
		// Production of protein -
		buffer.append("ELONGATION_TERMINATON_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("AR_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("RIBOSOME+");
		buffer.append(strGSymbol);
		buffer.append(",0,inf;\n");
		
		// add a new line -
		buffer.append("DEGRADATION_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("[],0,inf;\n");
		
		// add mRNA degradation -
		buffer.append("DEGRADATION_MRNA_");
		buffer.append(strGSymbol);
		buffer.append(",mRNA_");
		buffer.append(strGSymbol);
		buffer.append(",");
		buffer.append("[],0,inf;\n");
	}
	
	
	
}
