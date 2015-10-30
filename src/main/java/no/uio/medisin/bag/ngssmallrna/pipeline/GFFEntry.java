/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
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
    
    private                 String  seqID;
    private                 String  src;
    private                 String  type;
    private                 int     start;
    private                 int     stop;
    private                 float   score;
    private                 String  strand;
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
        }
        try{
            start   = Integer.parseInt(line.split("\t")[GFF_START]);
            stop    = Integer.parseInt(line.split("\t")[GFF_STOP]);
        }
        catch(Exception ex){
            logger.error("exception while parsing start/stop values in GFF entry" + line);
            logger.error(ex);       
            
        }
        if(line.split("\t")[GFF_SCORE].isEmpty() == false)
            score   = Float.parseFloat(line.split("\t")[GFF_SCORE]);
        
        strand  = line.split("\t")[GFF_STRAND];
        
        if(line.split("\t")[GFF_PHASE].isEmpty() == false)
            phase   = Integer.parseInt(line.split("\t")[GFF_PHASE]);
        
        if(line.split("\t")[GFF_ATTR].isEmpty() == false)
            attr    = line.split("\t")[GFF_ATTR];
        
        
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
    public String getStrand() {
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
}
