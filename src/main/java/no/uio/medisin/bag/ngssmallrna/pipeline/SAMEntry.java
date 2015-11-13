/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * parse a SAM line
 * 
 * the line can be two types. A read entry or a data line (starts with a @)
 * 
 *
    1   QNAME	   Query template NAME
    2   FLAG	   bitwise FLAG
    3   RNAME	   Reference sequence NAME
    4   POS	   1-based leftmost mapping POSition
    5   MAPQ	   MAPping Quality
    6   CIGAR	   CIGAR string
    7   RNEXT	   Ref. name of the mate/next read
    8   PNEXT	   Position of the mate/next read
    9   TLEN	   observed Template LENgth
    10  SEQ	   segment SEQuence
    11  QUAL	   ASCII of Phred-scaled base QUALity+33
 *                       
 *
 *
 * @author simonray
 */
public final class SAMEntry {
    
    static Logger                       logger                      = LogManager.getLogger();
    
    
    private static final int            QNAME   = 0;
    private static final int            FLAG    = 1;
    private static final int            RNAME   = 2;
    private static final int            POS     = 3;
    private static final int            MAPQ    = 4;
    private static final int            CIGAR   = 5;
    private static final int            RNEXT   = 6;
    private static final int            PNEXT   = 7;
    private static final int            TLEN    = 8;
    private static final int            SEQ     = 9;
    private static final int            QUAL    = 10;
    private static final int            MDSTR   = 11;
    
    private static final short          UNMAPPED  = 0x04;

    private String                      qName;
    private short                       flags;
    private String                      rName;
    private int                         startPos;
    private int                         endPos;
    private String                      strand;
    private int                         mapQ;
    private String                      cigar;
    private String                      rNext;
    private int                         pNext;
    private int                         tLen;
    private String                      seq;
    private String                      qual;
    private String                      mdString;
    private Boolean                     header;
    private Boolean                     mapped;
    
    private Boolean                     fullInfo = true;
    public SAMEntry(String samLine){
        parseLine(samLine);
    }
    
    /**
     * is line a header line?
     * 
     * @param line
     * @return 
     */
    public static final Boolean isHeaderLine(String line){
        return line.startsWith("@");
    }
    
    
    /**
     * extract information from input line
     * at the moment, this only extracts the information I am interested in
     * 
     * @param samLine String 
     * 
     */
    public void parseLine(String samLine){

        try{
            if(samLine.startsWith("@")) {
                header = true;
                return;
               
            }
 
            header = false;
            flags = Short.parseShort(samLine.split("\t")[FLAG]);
            if ((flags & UNMAPPED) == UNMAPPED) {
                mapped = false;
                return;                
            }
            
            mapped = true;


            qName = samLine.split("\t")[QNAME];
            
            if(samLine.split("\t")[FLAG].equals("16") || samLine.split("\t")[FLAG].equals("0")){
                if (samLine.split("\t")[FLAG].equals("16")) {
                    strand = "-";
                }
                else{
                    strand = "+";
                }

                mapQ = Integer.parseInt(samLine.split("\t")[MAPQ]);
                startPos = Integer.parseInt(samLine.split("\t")[POS]);
                String cigarStr = samLine.split("\t")[CIGAR].replace("M", "").trim();
                endPos = getStartPos() + Integer.parseInt(cigarStr) - 1;
                rName= samLine.split("\t")[RNAME].trim();
                mdString = samLine.split("\t")[MDSTR];
                cigar = samLine.split("\t")[CIGAR];
                rNext = samLine.split("\t")[RNEXT];
                pNext = Integer.parseInt(samLine.split("\t")[PNEXT]);
                tLen = Integer.parseInt(samLine.split("\t")[TLEN]);
                if(fullInfo){
                    seq = samLine.split("\t")[SEQ];
                    qual = samLine.split("\t")[QUAL];
                }
            }
        }
        catch(Exception ex){
            logger.error("error parsing samLine " + samLine);
            logger.error(ex);
        }
        
    }

    /**
     * @return the qName
     */
    public String getqName() {
        return qName;
    }

    /**
     * @param qName the qName to set
     */
    public void setqName(String qName) {
        this.qName = qName;
    }

    /**
     * @return the flags
     */
    public short getFlags() {
        return flags;
    }

    /**
     * @return the rName
     */
    public String getrName() {
        return rName;
    }

    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @return the mapQ
     */
    public int getMapQ() {
        return mapQ;
    }

    /**
     * @return the rNext
     */
    public String getrNext() {
        return rNext;
    }

    /**
     * @return the tLen
     */
    public int gettLen() {
        return tLen;
    }

    /**
     * @return the seq
     */
    public String getSeq() {
        return seq;
    }

    /**
     * @return the qual
     */
    public String getQual() {
        return qual;
    }

    /**
     * @return the mdString
     */
    public String getMdString() {
        return mdString;
    }

    /**
     * @return the cigar
     */
    public String getCigar() {
        return cigar;
    }

    /**
     * @return the pNext
     */
    public int getpNext() {
        return pNext;
    }

    /**
     * @return the header
     */
    public Boolean isHeaderLine() {
        return header;
    }

    /**
     * @return the mapped
     */
    public Boolean isMappedRead() {
        return mapped;
    }

}
