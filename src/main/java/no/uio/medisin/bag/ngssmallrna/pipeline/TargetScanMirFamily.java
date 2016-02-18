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
public class TargetScanMirFamily {
    
    static Logger                       logger                          = LogManager.getLogger();
    
    
    public final static int             MIRFAMILY_COL       = 0;
    public final static int             SEED_COL            = 1;
    public final static int             SPECIESID_COL       = 2;
    public final static int             MIRBASEID_COL       = 3;
    public final static int             SEQ_COL             = 4;
    public final static int             CONSERVATION_COL    = 5;
    public final static int             MIRBASEACC_COL      = 6;
    
    
    
    private String                      miRFamily;
    private String                      seedM8;
    private String                      speciesID;
    private String                      miRBaseID;
    private String                      matureSequence;
    private String                      familyConservation;
    private String                      miRBaseAccession;

    
    public TargetScanMirFamily(String familyLine){
        
        miRFamily           = familyLine.split("\t")[MIRFAMILY_COL].trim();
        seedM8              = familyLine.split("\t")[SEED_COL].trim();
        speciesID           = familyLine.split("\t")[SPECIESID_COL].trim();
        miRBaseID           = familyLine.split("\t")[MIRBASEID_COL].trim();
        matureSequence      = familyLine.split("\t")[SEQ_COL].trim();
        familyConservation  = familyLine.split("\t")[CONSERVATION_COL].trim();
        if (familyLine.split("\t").length == 7)
            miRBaseAccession    = familyLine.split("\t")[MIRBASEACC_COL].trim();
        
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
     * @return the seedM8
     */
    public String getSeedm8() {
        return seedM8;
    }

    /**
     * @param seedM8 the Seedm8 to set
     */
    public void setSeedm8(String seedM8) {
        this.seedM8 = seedM8;
    }

    /**
     * @return the SpeciesID
     */
    public String getSpeciesID() {
        return speciesID;
    }

    /**
     * @param SpeciesID the SpeciesID to set
     */
    public void setSpeciesID(String SpeciesID) {
        this.speciesID = SpeciesID;
    }

    /**
     * @return the miRBaseID
     */
    public String getMiRBaseID() {
        return miRBaseID;
    }

    /**
     * @param miRBaseID the miRBaseID to set
     */
    public void setMiRBaseID(String miRBaseID) {
        this.miRBaseID = miRBaseID;
    }

    /**
     * @return the matureSequence
     */
    public String getMatureSequence() {
        return matureSequence;
    }

    /**
     * @param matureSequence the matureSequence to set
     */
    public void setMatureSequence(String matureSequence) {
        this.matureSequence = matureSequence;
    }

    /**
     * @return the familyConservation
     */
    public String getFamilyConservation() {
        return familyConservation;
    }

    /**
     * @param familyConservation the familyConservation to set
     */
    public void setFamilyConservation(String familyConservation) {
        this.familyConservation = familyConservation;
    }

    /**
     * @return the miRBaseAccession
     */
    public String getMiRBaseAccession() {
        return miRBaseAccession;
    }

    /**
     * @param miRBaseAccession the miRBaseAccession to set
     */
    public void setMiRBaseAccession(String miRBaseAccession) {
        this.miRBaseAccession = miRBaseAccession;
    }
}
