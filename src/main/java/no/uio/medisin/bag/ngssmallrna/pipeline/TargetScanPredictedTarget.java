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
public class TargetScanPredictedTarget {
    
    static Logger                       logger                          = LogManager.getLogger();
    
    public final static int             MIRFAMILY_COL       = 0;
    public final static int             GENEID_COL          = 1;
    public final static int             GENESYMB_COL        = 2;
    public final static int             TRANSID_COL         = 3;
    public final static int             SPECIESID_COL       = 4;
    public final static int             UTTSTART_COL        = 5;
    public final static int             UTREND_COL          = 6;
    public final static int             MSASTART_COL        = 7;
    public final static int             MSASTOP_COL         = 8;
    public final static int             SEEDMATCH_COL       = 9;
    public final static int             PCT_COL             = 10;

    private String                      miRFamily;	
    private String                      GeneID;	
    private String                      GeneSymbol;	
    private String                      TranscriptID;	
    private int                         SpeciesID;	
    private int                         UTRstart;	
    private int                         UTRend;	
    private int                         MSAstart;	
    private int                         MSAend;	
    private String                      Seedmatch;
    private String                      PCT;

    
    
    public TargetScanPredictedTarget(String targetLine){
        
        miRFamily       = targetLine.split("\t")[MIRFAMILY_COL].trim();
        GeneID          = targetLine.split("\t")[GENEID_COL].trim();
        GeneSymbol      = targetLine.split("\t")[GENESYMB_COL].trim();
        TranscriptID    = targetLine.split("\t")[TRANSID_COL].trim();
        SpeciesID       = Integer.parseInt(targetLine.split("\t")[SPECIESID_COL].trim());
        UTRstart        = Integer.parseInt(targetLine.split("\t")[UTTSTART_COL].trim());
        UTRend          = Integer.parseInt(targetLine.split("\t")[UTREND_COL].trim());
        MSAstart        = Integer.parseInt(targetLine.split("\t")[MSASTART_COL].trim());
        MSAend          = Integer.parseInt(targetLine.split("\t")[MSASTOP_COL].trim());
        Seedmatch       = targetLine.split("\t")[SEEDMATCH_COL].trim();
        PCT             = targetLine.split("\t")[PCT_COL].trim();
        
    }
    
    
    
    /**
     * @return the miRFamily
     */
    public String getMiRFamily() {
        return miRFamily;
    }

    /**
     * @param miRFamily the miRFamily to set
     */
    public void setMiRFamily(String miRFamily) {
        this.miRFamily = miRFamily;
    }

    /**
     * @return the GeneID
     */
    public String getGeneID() {
        return GeneID;
    }

    /**
     * @param GeneID the GeneID to set
     */
    public void setGeneID(String GeneID) {
        this.GeneID = GeneID;
    }

    /**
     * @return the GeneSymbol
     */
    public String getGeneSymbol() {
        return GeneSymbol;
    }

    /**
     * @param GeneSymbol the GeneSymbol to set
     */
    public void setGeneSymbol(String GeneSymbol) {
        this.GeneSymbol = GeneSymbol;
    }

    /**
     * @return the TranscriptID
     */
    public String getTranscriptID() {
        return TranscriptID;
    }

    /**
     * @param TranscriptID the TranscriptID to set
     */
    public void setTranscriptID(String TranscriptID) {
        this.TranscriptID = TranscriptID;
    }

    /**
     * @return the SpeciesID
     */
    public int getSpeciesID() {
        return SpeciesID;
    }

    /**
     * @param SpeciesID the SpeciesID to set
     */
    public void setSpeciesID(int SpeciesID) {
        this.SpeciesID = SpeciesID;
    }

    /**
     * @return the UTRstart
     */
    public int getUTRstart() {
        return UTRstart;
    }

    /**
     * @param UTRstart the UTRstart to set
     */
    public void setUTRstart(int UTRstart) {
        this.UTRstart = UTRstart;
    }

    /**
     * @return the UTRend
     */
    public int getUTRend() {
        return UTRend;
    }

    /**
     * @param UTRend the UTRend to set
     */
    public void setUTRend(int UTRend) {
        this.UTRend = UTRend;
    }

    /**
     * @return the MSAstart
     */
    public int getMSAstart() {
        return MSAstart;
    }

    /**
     * @param MSAstart the MSAstart to set
     */
    public void setMSAstart(int MSAstart) {
        this.MSAstart = MSAstart;
    }

    /**
     * @return the MSAend
     */
    public int getMSAend() {
        return MSAend;
    }

    /**
     * @param MSAend the MSAend to set
     */
    public void setMSAend(int MSAend) {
        this.MSAend = MSAend;
    }

    /**
     * @return the Seedmatch
     */
    public String getSeedmatch() {
        return Seedmatch;
    }

    /**
     * @param Seedmatch the Seedmatch to set
     */
    public void setSeedmatch(String Seedmatch) {
        this.Seedmatch = Seedmatch;
    }

    /**
     * @return the PCT
     */
    public String getPCT() {
        return PCT;
    }

    /**
     * @param PCT the PCT to set
     */
    public void setPCT(String PCT) {
        this.PCT = PCT;
    }
    
    
}
