package org.varnerlab.ccmlparser;

// Import statements
import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

/**
 *  Helper class that contains a set of simple static io routines.
 *  @author J.Varner
 */
public class GIOL extends Object {
    
   /** 
    *  Public static method holding logic for reading a file. 
    *  @returns StringBuffer
    *  @throws Exception
    */
    public static StringBuffer read(String path) throws Exception {
        // method attributes
        StringBuffer buffer=new StringBuffer();
        
        // Create new buffered reader and load file
        File iFile=new File(path);
        BufferedReader reader=new BufferedReader(new FileReader(iFile));
        String s="";
        while ((s=reader.readLine())!=null){
            // Check for comments
            int idxCPP=s.indexOf("//");
            int idxCOpen=s.indexOf("/*");
            int idxCClose=s.indexOf("*/");
            boolean idxStar=s.startsWith("*");
            int idxMatlab=s.indexOf("%");
            int whiteSpace=s.length();
            
            // If any of the comments indexes are not -1 then we have a comment line
            if (idxCPP==-1 && idxCOpen==-1 && idxCClose==-1 && !idxStar && idxMatlab==-1){
                if (whiteSpace!=0){
                    buffer.append(s);
                    //buffer.append("\n");
                }
            }
        }
        
        // Close reader
        reader.close();
        
        // return buffer to caller
        return(buffer);
    }
    
    public static StringBuffer readNewLine(String path) throws Exception {
        // method attributes
        StringBuffer buffer=new StringBuffer();
        
        // Create new buffered reader and load file
        File iFile=new File(path);
        BufferedReader reader=new BufferedReader(new FileReader(iFile));
        String s="";
        while ((s=reader.readLine())!=null){
            // Check for comments
            int idxCPP=s.indexOf("//");
            int idxCOpen=s.indexOf("/*");
            int idxCClose=s.indexOf("*/");
            boolean idxStar=s.startsWith("*");
            int idxMatlab=s.indexOf("%");
            int whiteSpace=s.length();
            
            // If any of the comments indexes are not -1 then we have a comment line
            if (idxCPP==-1 && idxCOpen==-1 && idxCClose==-1 && !idxStar && idxMatlab==-1){
                if (whiteSpace!=0){
                    buffer.append(s);
                    buffer.append("\n");
                }
            }
        }
        
        // Close reader
        reader.close();
        
        // return buffer to caller
        return(buffer);
    }
    
    
    /**
     *  Public static merthod that write StringBuffer to disk. Take two ars, the path and the buffer.
     *  @param String Path
     *  @param StringBuffer My Payload (JT Rules!)
     *  @throws Exception
     */
    public static void write(String path,StringBuffer buffer) throws Exception {
        // Create writer
        File oFile=new File(path);
        BufferedWriter writer=new BufferedWriter(new FileWriter(oFile));
        
        // Write buffer to file system and close writer
        writer.write(buffer.toString());
        writer.close();
    }
    
    public static void write(String path,String buffer) throws Exception {
        // Create writer
        File oFile=new File(path);
        BufferedWriter writer=new BufferedWriter(new FileWriter(oFile));
        
        // Write buffer to file system and close writer
        writer.write(buffer);
        writer.close();
    }
}