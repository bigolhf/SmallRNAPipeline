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
 * Undefined fields are replaced with the "." character, as described in the original GFF spec.
 *
 * Column 1: "seqid"
 * Column 2: "source"
 * Column 3: "type"
 * Columns 4 & 5: "start" and "end"
 * Column 6: "score"
 * Column 7: "strand"
 * Column 8: "phase"
 * Column 9: "attributes"
 * A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. 
 * These tags have predefined meanings:
 * 
 *   ID            - Indicates the ID of the feature. must be unique within the scope of the GFF file. 
 *   Name          - Display name for the feature. not necessarily unique
 *   Alias         - A secondary name for the feature. e.g. locus names and accession numbers. not necessarily unique
 *   Parent        - Indicates the parent of the feature. 
 *   Target        - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. 
 *   Gap           - The alignment of the feature to the target if the two are not collinear 
 *   Derives_from  - Used to disambiguate the relationship between one feature and another 
 *   Note          - A free text note.
 *   Dbxref        - A database cross reference. 
 *   Ontology_term - A cross reference to an ontology term. 
 *   Is_circular   - A flag to indicate whether a feature is circular. 
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
        
    static String   CR = System.getProperty("line.separator");
    
    
    private                 String  refSeqID;
    private                 String  featureID;
    private                 String  src;
    private                 String  type;
    private                 int     start;
    private                 int     stop;
    private                 float   score;
    private                 Strand  strand;
    private                 int     phase;
    private                 String  attrString;
    
    
    public GFFEntry(String line){
        
        try{
            refSeqID    = line.split("\t")[GFF_SEQID];
            src         = line.split("\t")[GFF_SRC];
            type        = line.split("\t")[GFF_TYPE];
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
            attrString    = line.split("\t")[GFF_ATTR];
            String attribs[] = attrString.split(";");
            for(String a:attribs){
                if(a.toUpperCase().contains("ID=")){
                    featureID=a.split("=")[1];
                }
            }
        }

        
        
        
    }
    
    
    
    
    
    
    /**
     * create a new GFF Entry.
     * We require all fields to be specified because there are too many possible
     * subsets of parameters that might be specified
     * 
     * @param nRefSeqID
     * @param nSource
     * @param nType
     * @param nStart
     * @param nStop
     * @param nScore
     * @param nStrand
     * @param nPhase
     * @param nAttrString 
     * 
     */
    public GFFEntry(String nRefSeqID, String nSource, String nType, int nStart, int nStop, String nScore, String nStrand, String nPhase, String nAttrString){
        
        refSeqID = nRefSeqID;
        src = nSource;
        type = nType;
        start = nStart;
        stop = nStop;
        try{
            score   = Float.parseFloat(nScore);
        }
        catch(NumberFormatException exNF){
            score = 0;
        }
        strand = StrandString.guessStrand(nStrand);
        try{
            phase   = Integer.parseInt(nPhase);
        }
        catch(NumberFormatException exNF){
            phase = 0;
        }
        attrString = nAttrString;
        if(nAttrString.isEmpty() == false){
            String attribs[] = attrString.split(";");
            for(String a:attribs){
                if(a.toUpperCase().contains("ID=")){
                    featureID=a.split("=")[1];
                }
            }
        }
    }

    
    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param qStart
     * @param qStop
     * @param qStrand
     * @param qChr
     * @param featureName
     * @param bleed
     * @return 
     */
    public Boolean doesRegionOverlap(int qStart, int qStop, Strand qStrand, String qChr, String featureName, int bleed){
        return 
          this.getType().toUpperCase().equals(featureName.toUpperCase())  
            && this.getStrand() == qStrand
            && this.getSeqID().equals(qChr)
            && Math.abs(this.getStart()- qStart) < bleed 
            && Math.abs(this.getStop() - qStop) < bleed;
    }
    
    

    
    
    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param qStart
     * @param qStop
     * @param qStrand
     * @param qChr
     * @param featureName
     * @return 
     */
    public Boolean doesFeatureContainRegion(int qStart, int qStop, Strand qStrand, String qChr, String featureName){
        return 
          this.getType().toUpperCase().equals(featureName.toUpperCase())  
            && this.getStrand() == qStrand
            && this.getSeqID().equals(qChr)
            && this.getStart() <= qStart 
            && this.getStop() >= qStop;
    }




    
    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param queryGFFEntry
     * @return 
     */
    public Boolean doesFeatureContainRegion(GFFEntry queryGFFEntry){
        return 
           this.getStrand() == queryGFFEntry.getStrand()
            && this.getSeqID().equals(queryGFFEntry.getSeqID())
            && this.getStart() <= queryGFFEntry.getStart()
            && this.getStop() >= queryGFFEntry.getStop();
    }




    
    /**
     * does the specified region overlap this gffEntry regardless of featureType?
     * 
     * @param qStart
     * @param qStop
     * @param qStrand
     * @param qChr
     * @param bleed
     * @return 
     */
    public Boolean doesRegionOverlap(int qStart, int qStop, Strand qStrand, String qChr, int bleed){
        return 
            this.getStrand() == qStrand
            && this.getSrc().equals(qChr)
            && Math.abs(this.getStart()- qStart) < bleed 
            && Math.abs(this.getStop() - qStop) < bleed;
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
        String gff3String = this.refSeqID + "\t"
                + this.src + "\t"
                + this.type + "\t"
                + start + "\t"
                + stop + "\t"
                + "." + "\t"
                + strand + "\t"
                + "." + "\t"
                + "ID=" + featureID + ";"  + "Seq=" + this.getSequence();
        return gff3String;
    }
        
    

    /**
     * add attribute to attribute string
     * 
     * @param attrKey
     * @param attrVal 
     */
    public void addAttr(String attrKey, String attrVal){
        attrString = attrString.concat(";" + attrKey + "=" + attrVal );            
    }   
    
    
    
    /**
     * return the specified Attribute value
     * 
     * @param attrKey
     * @return 
     */
    public String getAttrValue(String attrKey){
        String attrs[] = attrString.split(";");
        for (String attr: attrs){
            if(attr.contains(attrKey)){
                return attr.split("=")[1].trim();
            }
        }
        return null;
    }
    
    /**
     * @return the seqID
     */
    public String getSeqID() {
        return refSeqID;
    }

    /**
     * @param seqID the seqID to set
     */
    public void setSeqID(String seqID) {
        this.refSeqID = seqID;
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
        return attrString;
    }
    
    
    
    /**
     * add a sequence by appending to the attribute string
     * 
     * @param seq 
     */
    public void setSequence(String seq){
        attrString = attrString.concat(";seq=" + seq);            
    }
    
    
    
    /**
     * extract sequence from the attribute string
     * 
     * @return 
     */
    public String getSequence(){
        if(attrString.contains("seq=")){
            int startPos = attrString.indexOf("seq=")+4;
            int stopPos = attrString.indexOf(";", startPos);
            if(stopPos==-1)
                stopPos = attrString.length() - 1;
            return attrString.substring(startPos, stopPos);
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
        if(attrString.contains("seq=")){            
            int startPos = attrString.indexOf("seq=")+4;
            int stopPos = attrString.indexOf(";", startPos);
            if (stopPos==-1)
                stopPos=attrString.length()-1;
            return ">" + featureID + "|" + this.refSeqID  + ":" + start + "-" + stop + "(" + strand + ")" + CR + attrString.substring(startPos, stopPos);
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
