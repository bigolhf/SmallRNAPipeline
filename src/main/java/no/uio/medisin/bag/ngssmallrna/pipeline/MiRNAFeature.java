/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import no.uio.medisin.bag.jmirpara.SimpleSeq;
/**
 * stores an miRBase entry
 * This needs to be made comparable so the entries stored in a list can be sorted
 * to speed things up
 * 
 * @author sr
 */
public class MiRNAFeature {

    static Logger                       logger                      = LogManager.getLogger();
    
    private String mimatID;
    private String note;
    private String name;
    private String parent;
    private String chromosome;
    private int    startPos;
    private int    endPos;
    private String strand;
    private String sequence;
    private String isomiRString;
    
    
    //(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t")
    
    private static final int nameCol    = 0;
    private static final int startCol   = 1;
    private static final int cigarCol   = 2;
    private static final int mdCol      = 3;
    private static final int seqCol     = 4;
    
    

    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for predicted miRNAs)
     * 
     * @param n String      name
     * @param c String      chromosome
     * @param s int         start position
     * @param e int         end position
     * @param t String      Strand
     * 
     */
    public MiRNAFeature(String n, String c, int s, int e, String t){
        this(n, c, s, e, t, "", "", "");
    }
    
    
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for miRBase entry miRNAs)
     * 
     * @param n String      name
     * @param c String      chromosome
     * @param s int         start position
     * @param e int         end position
     * @param t String      Strand
     * @param m String      MIMAT (miRBase) ID
     * @param p String      MI (miRBase) Parent ID
     * @param seq :     String sequence
     * 
     */
    public MiRNAFeature(String n, String c, int s, int e, String t, String m, String p, String seq){
        name = n;
        chromosome = c;
        startPos = s;
        endPos = e;
        strand = t;
        parent = p;
        mimatID = m;
        sequence = seq;
        isomiRString = "";
        
    }
    
    
    public MiRNAFeature(MiRNAFeature m){
        
        name            = m.name;
        chromosome      = m.chromosome;
        startPos        = m.startPos;
        endPos          = m.endPos;
        strand          = m.strand;
        parent          = m.parent;
        mimatID         = m.mimatID;
        sequence        = m.sequence;
        isomiRString    = m.isomiRString;
        
    }
    /**
     * checks whether Chromosome strings are the same, while attempting
     * to allow for the presence or absence of a variation on the 'Chr' 
     * prefix
     * 
     * @param queryChr
     * @return 
     */
    public Boolean chromosomeMatch(String queryChr){
        
        return MiRNAFeature.removeChromosomePrefix(chromosome).equals(MiRNAFeature.removeChromosomePrefix(queryChr));
        
    }
    
    
    
    
    /**
     * attempt to remove any prefix of the form 'Chr' from the chromosome string
     * 
     * @param chrString
     * @return 
     */
    public static String removeChromosomePrefix(String chrString){
        
        if(chrString.contains("chr")){
            chrString = chrString.replace("chr", "");
        }
        else{
            if(chrString.contains("CHR")){
                chrString = chrString.replace("CHR", "");
            }
            else{
                if(chrString.contains("Chr")){
                    chrString = chrString.replace("Chr", "");
                }
            }
            
        }
        return chrString;
        
    }
        
    
    /**
     * add information to define an isomiR for this entry
     * 
     * @param name
     * @param start
     * @param cigar
     * @param md 
     * @param seq
     * 
     */
    public void addIsomiR(String name, int start, String cigar, String md, String seq){

        isomiRString = isomiRString.concat(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t");
        
    }
    
    /**
     * delete the isomiR string
     * 
     */
    public void removeIsomiRs(){
        isomiRString = "";
    }

    
    /**
     * write out isomiRs.
     * 
     * only report reads that are above a baseline, defined in terms of the fraction 
     * of the total number of reads for the miRNA. e.g. if there are 100 reads, and 
     * baseline is 5, then only isomiRs with more than 5 reads will be reported
     * 
     * @param baselinePercent : int     only report isomiRs with reads that
     * @param minCounts       : int     total counts for isomiR must be greater
     *                                  than this value
     * @return 
     */
    public String reportIsomiRs(int baselinePercent, int minCounts){
        String reportStr = this.getName() + ":\t[" + this.getChromosome() + "]\t" 
                + this.getStartPos() + "\t" + this.getEndPos() + this.getSequence() + "\n";
        String [] isomiRs = isomiRString.split("\t");
        
        int totalCounts = this.getTotalCounts();
        if (totalCounts < minCounts) return "";
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        
        String isomiRStr = "";
        for(String isomiR: isomiRs){
            String[] values = isomiR.split(";");
            if(Double.parseDouble(isomiR.split(";")[nameCol].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                isomiRStr = isomiRStr.concat("name: " + values[0] + "\t"
                            + "start: " + values[startCol] + "\t"
                            + "cigar: " + values[cigarCol] + "\t"
                            + "MD: " + values[mdCol] + "\t" 
                            + "SQ: " + values[seqCol] + "\n"
                );
            }
        }
        
        if (isomiRStr.equals("")) return "";
        return reportStr.concat(isomiRStr);
    }
    
    
    
    
    /**
     * report the isomiR in a manner that is visually appealing
     * 
     * @param baselinePercent
     * @param minCounts
     * @return 
     */
    public String prettyReportIsomiRs(int baselinePercent, int minCounts){
        
        String reportStr = this.getName() + "|" + this.getMimatID() + " :\tchr" + this.getChromosome() + "\t" 
                + this.getStartPos() + " --> " + this.getEndPos() + " (" + this.getStrand() + ") : " + this.getSequence() + "\n";
        
        int totalCounts = this.getTotalCounts();
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        logger.debug(reportStr);
        if (totalCounts < minCounts) return "";
        
    
        String [] isomiRs = isomiRString.split("\t");
        
        int longestName = this.getName().length();
        int longestSeq = this.getSequence().length();
        int longestCounts = 0;
        int longestMD = 0;
        int minStart = this.getIsomiRMinStart();
        int maxStop = this.getIsomiRMaxStop();
        
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[nameCol].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                if(isomiR.split(";")[nameCol].length() > longestName)
                    longestName = isomiR.split(";")[nameCol].length();

                if(isomiR.split(";")[nameCol].split("-")[1].length() > longestCounts)
                    longestCounts = isomiR.split(";")[nameCol].split("-")[1].length();

                if(isomiR.split(";")[seqCol].length() > longestSeq)
                    longestSeq = isomiR.split(";")[seqCol].length();
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[seqCol].equals(this.getSequence().replace("U", "T")))
                        longestName++;
                }
                else{
                    if(isomiR.split(";")[seqCol].equals(SimpleSeq.complement(this.getSequence()).replace("U", "T")))
                        longestName++;                    
                }
                
                    

                if(isomiR.split(";")[mdCol].length() > longestMD)
                    longestMD = isomiR.split(";")[mdCol].length();
                
                
            }          

        }    
        
        int leftMargin = 10;
        int ColMargin = 5;
        
        String isoString = "";
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[nameCol].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){

                isoString           = isoString.concat(StringUtils.repeat(" ", leftMargin));
                String isoName      = isomiR.split(";")[nameCol];
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[seqCol].equals(this.getSequence().replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                else{
                    if(isomiR.split(";")[seqCol].equals(SimpleSeq.complement(this.getSequence()).replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                isoString           = isoString.concat(StringUtils.repeat(" ", longestName - isoName.length()) + isoName);
                
                int isoStart        = Integer.parseInt(isomiR.split(";")[startCol]);
                String isoSeq       = isomiR.split(";")[seqCol];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (isoStart - this.getIsomiRMinStart())));
                isoString           = isoString.concat(isoSeq);
                isoString           = isoString.concat(StringUtils.repeat(" ", (this.getIsomiRMaxStop()-(isoStart + isoSeq.length())) + ColMargin));
                
                String isoMD        = isomiR.split(";")[mdCol];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin) + isoMD + StringUtils.repeat(" ", (longestMD - isoMD.length()) + ColMargin));
                
                String isoCounts    = isomiR.split(";")[nameCol].split("-")[1];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (longestCounts - isoCounts.length())) + isoCounts + StringUtils.repeat(" ", ColMargin) + "\n");
            
            }
            
        }
        
        if (isoString.equals("")) return "";
        
        return reportStr.concat(isoString + "\n\n\n");

    }
    
    
    
    
        /**
     * an isomiR can be classified as 5´ modification , 3´ modification or polymorphic
     * 
     * @param baselinePercent
     * @param minCounts 
     * 
     * @return ArrayList    : list of isomiR points
     */
    public ArrayList characterizeIsomiRs(int baselinePercent){
        
        ArrayList isomiRPts = new ArrayList<>();
        // we can identify 5´ modification from start position
        int totalCounts = this.getTotalCounts();
        String [] isomiRStrings = isomiRString.split("\t");
        String isoString = "";
        logger.debug(this.name);

        for(String isomiR: isomiRStrings){

            
            if(Double.parseDouble(isomiR.split(";")[nameCol].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                logger.debug(isomiR.split(";")[nameCol] +  this.getStrand() +  "\n" + " : " + Double.parseDouble(isomiR.split(";")[nameCol].split("-")[1]) / (double) totalCounts);


                //logger.info(isomiR);
                int isoStart        = Integer.parseInt(isomiR.split(";")[startCol]);
                String isoSeq       = isomiR.split(";")[seqCol];
                int isoEnd          = isoStart + isoSeq.length() - 1;

                int noOf5pSteps = 0;
                int noOf3pSteps = 0;
                int noOfPolySteps = 0;                

                int wStart = 0;
                int wEnd = 0;
                int iStart = 0;
                int iEnd = 0;

                
                if(this.getStrand().equals("+")){
                    logger.debug(this.sequence + "\n" + isomiR.split(";")[seqCol]);
                    if(isoStart != this.startPos){                   
                        noOf5pSteps = isoStart - this.startPos;                   
                    }

                    if(isoEnd != this.endPos){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug(noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug(isoStart + ", " + isoEnd + ", " + (this.startPos- isoStart) + ", " + sequence.length());

                    if(this.startPos >= isoStart ){
                        wStart = 0;
                        iStart = this.startPos- isoStart;
                    }
                    else{
                        wStart = isoStart - this.startPos;
                        iStart = 0;
                    }

                    if(this.endPos < isoEnd){
                        wEnd = sequence.length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.endPos);
                    }
                    else{
                        wEnd = sequence.length() - (this.endPos - isoEnd);
                        iEnd = isoSeq.length();
                    }

                    if(sequence.substring(wStart, wEnd).replace("U", "T").equals(isoSeq.substring(iStart, iEnd))==false){
                        for(int b=0; b<wEnd-wStart; b++){
                            if(sequence.charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }
                }
                else{
                    logger.debug(this.sequence + "\t" + startPos + "\t" + endPos + "\n");
                    logger.debug(SimpleSeq.complement(isomiR.split(";")[seqCol]) + "\t" + isoStart + "\t" + isoEnd + "\n");
                    
                    if(isoStart != this.startPos){                   
                        noOf5pSteps = isoStart - this.startPos;                   
                    }

                    if(isoEnd != this.endPos){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug(noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug(isoStart + ", " + isoEnd + ", " + (this.startPos- isoStart) + ", " + sequence.length());

                    if(this.startPos >= isoStart ){
                        wStart = 0;
                        iStart = this.startPos- isoStart;
                    }
                    else{
                        wStart = isoStart - this.startPos;
                        iStart = 0;
                    }

                    if(this.endPos < isoEnd){
                        wEnd = sequence.length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.endPos);
                    }
                    else{
                        wEnd = sequence.length() - (this.endPos - isoEnd);
                        iEnd = isoSeq.length();
                    }

                    if(sequence.substring(wStart, wEnd).replace("U", "T").equals(SimpleSeq.complement(isoSeq.substring(iStart, iEnd)))==false){
                        String complementWTSeq = SimpleSeq.complement(sequence);
                        for(int b=0; b<wEnd-wStart; b++){
                            if(complementWTSeq.charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }                 
                    
                }


                String isoCounts    = isomiR.split(";")[nameCol].split("-")[1];
                
                HashMap isomiRPt = new HashMap();
                isomiRPt.put("5p", noOf5pSteps);
                isomiRPt.put("3p", noOf3pSteps);
                isomiRPt.put("poly", noOfPolySteps);
                isomiRPt.put("fraction", Double.parseDouble(isoCounts)/(double)totalCounts);
                
                isomiRPts.add(isomiRPt);

            }

        }
        
        return isomiRPts;

    }
    

    
    
    
    /**
     * sum counts for this isomiR
     * 
     * @return 
     */
    public int getTotalCounts(){
        
        int totalCounts = 0;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){            
            totalCounts += Integer.parseInt(isomiR.split(";")[nameCol].split("-")[1]);
        }
        return totalCounts;
    }
    
    
    /**
     * find the smallest start position within the isomiRs
     * 
     * @return 
     */
    public int getIsomiRMinStart(){
        
        int isomiRMinStart = 1000000000;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[startCol]) < isomiRMinStart) 
                isomiRMinStart = Integer.parseInt(isomiR.split(";")[startCol]);
        }
        return isomiRMinStart;
        
    }
    
    
    
    
    /**
     * find the smallest start position within the isomiRs
     * 
     * @return 
     */
    public int getIsomiRMaxStop(){
        
        int isomiRMaxStop = -1;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[startCol]) + isomiR.split(";")[seqCol].length() > isomiRMaxStop) 
                isomiRMaxStop = Integer.parseInt(isomiR.split(";")[startCol]) + isomiR.split(";")[seqCol].length();
        }
        return isomiRMaxStop;
        
    }
    
    
    
    
    /**
     * 
     * @param miRFeat
     * @return
     */
    @Deprecated
    public int compareTo(MiRNAFeature miRFeat) {

        int thisChr = -1;
        int queryChr = -1;
        if(this.getChromosome().contains("chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("chr", ""));
        }
        
        if(this.getChromosome().contains("CHR")){
            thisChr = Integer.parseInt(this.getChromosome().replace("CHR", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("CHR", ""));
        }
        
        if(this.getChromosome().contains("Chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("Chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("Chr", ""));
        }
        
     return thisChr - queryChr;

     }
   
    
    
    
    /**
     * base equality on the mimatID, which should be unique
     * 
     * @param qObject
     * @return 
     */
    @Override
    public boolean equals(Object qObject){
        if (qObject != null && qObject instanceof MiRNAFeature)
        {
            return (this.mimatID.equals(((MiRNAFeature) qObject).mimatID));
        }
        return false;

    }
    
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((mimatID == null) ? 0 : mimatID.hashCode());
        return result;
    }    
    
    
    /**
     * @return the mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    /**
     * @param mimatID the mimatID to set
     */
    public void setMimatID(String mimatID) {
        this.mimatID = mimatID;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the parent
     */
    public String getParent() {
        return parent;
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(String parent) {
        this.parent = parent;
    }

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @param startPos the startPos to set
     */
    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @param endPos the endPos to set
     */
    public void setEndPos(int endPos) {
        this.endPos = endPos;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(String strand) {
        this.strand = strand;
    }

    /**
     * @return the isomiRString
     */
    public String getIsomiRString() {
        return isomiRString;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }
    
}
