/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
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
    
    
    

    /**
     * 
     * return entry within the specified region
     * 
     * @param start
     * @param stop
     * @param strand
     * @param chr
     * @param bleed
     * @return 
     */
    public GFFEntry findMatch(int start, int stop, String strand, String chr, int bleed){
        Iterator itGF = GFFEntries.iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();
            if(gffEntry.getStrand().equals(strand)
                && gffEntry.getAttr().equals(chr)
                && Math.abs(gffEntry.getStart()-start) < bleed 
                && Math.abs(gffEntry.getStop() - stop) < bleed
                    ){
                return gffEntry;
            }
        }
        return null;
    }
    
    
    
    
    /**
     * Does the specified region contain a feature?
     * 
     * @param start
     * @param stop
     * @param strand
     * @param chr
     * @param bleed
     * @return 
     */
    public Boolean doesRegionContainFeature(int start, int stop, String strand, String chr, int bleed){
        Iterator itGF = GFFEntries.iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();
            if(gffEntry.getStrand().equals(strand)
                && gffEntry.getSrc().equals(chr)
                && Math.abs(gffEntry.getStart()-start) < bleed 
                && Math.abs(gffEntry.getStop() - stop) < bleed
                    ){
                return true;
            }
        }
        return false;
        
    }
    
    
    
    /**
     * write out the features in GFF format
     * 
     * @param bwFT
     * @throws IOException 
     */
    public void writeFeaturesAsGFF3(BufferedWriter bwFT) throws IOException{
        Iterator itGF = GFFEntries.iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();            
            bwFT.write(gffEntry.toGFF3String() + "\n");
        }
    }
}
