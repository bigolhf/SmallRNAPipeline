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
public class TargetScanConservedFamilyEntry {
    
    static Logger                       logger                          = LogManager.getLogger();
    
    private String                      miRFamily;
    private String                      Seedm8;
    private String                      SpeciesID;
    private String                      miRBaseID;
    private String                      matureSequence;
    private String                      familyConservation;
    private String                      miRBaseAccession;

    
    
    
    
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
     * @return the Seedm8
     */
    public String getSeedm8() {
        return Seedm8;
    }

    /**
     * @param Seedm8 the Seedm8 to set
     */
    public void setSeedm8(String Seedm8) {
        this.Seedm8 = Seedm8;
    }

    /**
     * @return the SpeciesID
     */
    public String getSpeciesID() {
        return SpeciesID;
    }

    /**
     * @param SpeciesID the SpeciesID to set
     */
    public void setSpeciesID(String SpeciesID) {
        this.SpeciesID = SpeciesID;
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
