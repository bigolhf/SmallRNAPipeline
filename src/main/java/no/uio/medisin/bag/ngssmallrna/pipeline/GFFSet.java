/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author simonray
 */
public class GFFSet {
    static Logger logger = LogManager.getLogger();
    private final ArrayList<GFFEntry> GFFEntries;
    
    public GFFSet(){
        GFFEntries = new ArrayList<>();
    }
    /**
     * add GFFEntry to the set
     * 
     * @param gffEntry 
     */
    public void addEntry(GFFEntry gffEntry){
        GFFEntries.add(gffEntry);
    }
    
    
    /**
     * read specified GFF file
     * 
     * @param filename
     * @return number of lines read
     * 
     */
    public int readGFF(String filename){
        String line = null;
        int lineCount = 0;
        try{
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            while ((line = br.readLine()) != null) {
                this.addEntry(new GFFEntry(line));
                lineCount ++;
            }
        }
        catch(IOException exIO){
            logger.error("error parsing GFF file " + filename);
            logger.error("exception thrown on line " + lineCount);
            logger.error(line);
            logger.error(exIO);
            return lineCount;
        }
        return lineCount;
    }
}
