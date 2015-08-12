/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * this stores the dispersion pattern for a set of isomiRs associated with a
 * miRNA.
 * 
 * This is kept separately from the miRNAFeature class to stop the class getting
 * too big and the overhead associated with creating an instance and in many
 * 
 * 
 * @author simon rayner
 */
public class IsomiRSet {
    private String mimatID;
    private String runID;
    
    
    private ArrayList<HashMap> isomiRPts;
    
    private ArrayList<IsomirPoint> isomiRPts2;

    
    public IsomiRSet(String mID, String rID, ArrayList<HashMap> isoPts){
        mimatID = mID;
        runID = rID;
        isomiRPts = isoPts;
    }
    
    
    
    /**
     * report the isomiR dispersion for this miRNA and experiment
     * 
     * @return String
     * 
     */
    public String tabReportIsomiRSet(){
        String reportStr = "";
        
        reportStr = reportStr.concat(mimatID + "\t" + runID + "\t" + isomiRPts.size() + "\n" ); 
        reportStr = reportStr.concat("5p" + "\t" + "3p" + "poly" + "\t" + "fraction" + "\n");
        for(HashMap hmIsomiR: isomiRPts){
            reportStr = reportStr.concat(hmIsomiR.get("5p") + "\t" + hmIsomiR.get("3p") 
                    + "\t" + hmIsomiR.get("poly") + "\t" + hmIsomiR.get("fraction") +"\n" );
        }
        return reportStr;
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
     * @return the runID
     */
    public String getRunID() {
        return runID;
    }

    /**
     * @param runID the runID to set
     */
    public void setRunID(String runID) {
        this.runID = runID;
    }

    /**
     * @return the isomiRPts
     */
    public ArrayList<HashMap> getIsomiRPts() {
        return isomiRPts;
    }

    /**
     * @param isomiRPts the isomiRPts to set
     */
    public void setIsomiRPts(ArrayList<HashMap> isomiRPts) {
        this.isomiRPts = isomiRPts;
    }

    /**
     * @return the isomiRPts2
     */
    public ArrayList<IsomirPoint> getIsomiRPts2() {
        return isomiRPts2;
    }

    /**
     * @param isomiRPts2 the isomiRPts2 to set
     */
    public void setIsomiRPts2(ArrayList<IsomirPoint> isomiRPts2) {
        this.isomiRPts2 = isomiRPts2;
    }
    
    
}
