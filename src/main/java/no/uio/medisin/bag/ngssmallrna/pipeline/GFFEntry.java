/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * stores a single entry from a GFF line
 * @author simonray
 */
public class GFFEntry {
    
    static                  Logger logger    = LogManager.getLogger();
    
    private static final    int GFF_SEQID  = 0;
    private static final    int GFF_SRC    = 1;
    private static final    int GFF_TYPE   = 2;
    private static final    int GFF_START  = 3;
    private static final    int GFF_STOP   = 4;
    private static final    int GFF_SCORE  = 5;
    private static final    int GFF_STRAND = 6;
    private static final    int GFF_PHASE  = 7;
    private static final    int GFF_ATTR   = 8;
        
    static String   CR = System.getProperty("line.separator");
    
    
    private                 String  seqID;
    private                 String  featureID;
    private                 String  src;
    private                 String  type;
    private                 int     start;
    private                 int     stop;
    private                 float   score;
    private                 Strand  strand;
    private                 int     phase;
    private                 String  attr;
    
    
    public GFFEntry(String line){
        
        try{
            seqID   = line.split("\t")[GFF_SEQID];
            src     = line.split("\t")[GFF_SRC];
            type    = line.split("\t")[GFF_TYPE];
        }
        catch(Exception ex){
            logger.error("exception while parsing seqID/src/type values in GFF entry" + line);
            logger.error(ex);
            //throw new Exception("exception while parsing seqID/src/type values in GFF entry" + line);
        }
        try{
            start   = Integer.parseInt(line.split("\t")[GFF_START]);
            stop    = Integer.parseInt(line.split("\t")[GFF_STOP]);
        }
        catch(Exception ex){
            logger.error("exception while parsing start/stop values in GFF entry" + line);
            logger.error(ex);      
            start = -1;
            stop  = -1;
            //throw new Exception("exception while parsing seqID/src/type values in GFF entry" + line);
        }
        try{
            score   = Float.parseFloat(line.split("\t")[GFF_SCORE]);
        }
        catch(NumberFormatException exNF){
            score = 0;
        }
        
        strand  = StrandString.guessStrand(line.split("\t")[GFF_STRAND]);
        
        try{
            phase   = Integer.parseInt(line.split("\t")[GFF_PHASE]);
        }
        catch(NumberFormatException exNF){
            phase = 0;
        }
        
        if(line.split("\t")[GFF_ATTR].isEmpty() == false){
            attr    = line.split("\t")[GFF_ATTR];
            String attribs[] = attr.split(";");
            for(String a:attribs){
                if(a.toUpperCase().contains("ID=")){
                    featureID=a.split("=")[1];
                }
            }
        }

        
        
        
    }
    
    
    /**
     * create a new GFFEntry from a limited information set
     * 
     * @param id
     * @param s
     * @param c
     * @param b
     * @param e 
     */
    public GFFEntry(String id, String s, String c, int b, int e){
        start = b;
        stop = e;
        src = c;
        strand = StrandString.guessStrand(s);
        seqID = id;
        attr = id;
    }

    

    public String toGFF3String(){
        /*
            chr1            chromosome
            source          n/a here
            miRNA           feature type (n/a)
            start pos
            end pos
            score           n/a here               
            strand          (+/-)
            frame           n/a here
            attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p;Derives_from=MI0022705

        */
        String gff3String = this.src + "\t"
                + "." + "\t"
                + "smallRNA" + "\t"
                + start + "\t"
                + stop + "\t"
                + "." + "\t"
                + strand + "\t"
                + "." + "\t"
                + "ID=" + seqID + ";" + "Alias=" + seqID + ";" + "Name=" + seqID;
        return gff3String;
    }
        
    

    /**
     * add attribute to attribute string
     * 
     * @param attrKey
     * @param attrVal 
     */
    public void addAttr(String attrKey, String attrVal){
        attr = attr.concat(";" + attrKey + "=" + attrVal );            
    }   
    
    

    /**
     * @return the seqID
     */
    public String getSeqID() {
        return seqID;
    }

    /**
     * @param seqID the seqID to set
     */
    public void setSeqID(String seqID) {
        this.seqID = seqID;
    }

    /**
     * @return the src
     */
    public String getSrc() {
        return src;
    }

    /**
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @return the stop
     */
    public int getStop() {
        return stop;
    }

    /**
     * @return the score
     */
    public float getScore() {
        return score;
    }

    /**
     * @return the strand
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * @return the phase
     */
    public int getPhase() {
        return phase;
    }

    /**
     * @return the attr
     */
    public String getAttr() {
        return attr;
    }
    
    
    
    /**
     * add a sequence by appending to the attribute string
     * 
     * @param seq 
     */
    public void setSequence(String seq){
        attr = attr.concat(";seq=" + seq);            
    }
    
    
    
    /**
     * extract sequence from the attribute string
     * 
     * @return 
     */
    public String getSequence(){
        if(attr.contains("seq=")){
            int startPos = attr.indexOf("seq=")+4;
            int stopPos = attr.indexOf(";", startPos);
            return attr.substring(startPos, startPos);
        }
        return "";
    }
    
    
    
    
    /**
     * write the entry as in FASTA format as used by MiRBase
     * for consistency with MiRBAse the header line must have the format
     * >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
     * 
     * @return 
     */
    public String toMirbaseFastAString(){
        if(attr.contains("seq=")){
            int startPos = attr.indexOf("seq=")+4;
            int stopPos = attr.indexOf(";", startPos);
            return ">" + this.seqID + " " + this.seqID.split("-")[1] + " " + this.seqID.split("-")[0] + " " + this.seqID.split("-")[1] 
                    + CR + attr.substring(startPos, stopPos);
        }
        return "";
    }

    /**
     * @return the featureID
     */
    public String getFeatureID() {
        return featureID;
    }

    /**
     * @param featureID the featureID to set
     */
    public void setFeatureID(String featureID) {
        this.featureID = featureID;
    }
}
