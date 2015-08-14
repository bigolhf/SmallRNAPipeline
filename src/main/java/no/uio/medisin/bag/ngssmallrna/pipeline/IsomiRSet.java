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
    
    private double totalPairwiseDist    = 0.0;
    private double total3pDist          = 0.0;
    private double total5pDist          = 0.0;
    private double totalPolyDist        = 0.0;
    
    
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
        
        reportStr = reportStr.concat(mimatID + ":" + runID + "\t" + isomiRPts.size() + "\t" + totalPairwiseDist 
          + "\t" + total5pDist + "\t" + total3pDist + "\t" + totalPolyDist + "\n" ); 
        reportStr = reportStr.concat("5p" + "\t" + "3p" + "poly" + "\t" + "fraction" + "\t" + "TPD" + "\n");
        for(HashMap hmIsomiR: isomiRPts){
            reportStr = reportStr.concat(hmIsomiR.get("5p") + "\t" + hmIsomiR.get("3p") 
                    + "\t" + hmIsomiR.get("poly") + "\t" + hmIsomiR.get("fraction") + "\n" );
        }
        return reportStr;
    }
    
    
    
    /**
     * calculate sum of pairwise distances between all isomiRs, scaled by the fraction that each isomiR is present
     * A rather rough way of estimating dispersion, there has to be a better way
     * 
     * 
     */
    public void calcDistParameters(){
        
        totalPairwiseDist = 0;
        for(int i=0; i<isomiRPts.size(); i++){
            
            HashMap hmIsomiRi = isomiRPts.get(i);
            for(int j=i+1; j<isomiRPts.size(); j++){
                
                HashMap hmIsomiRj = isomiRPts.get(j);
                //(double) ((Integer) marks.get(i)).intValue();
                ((Integer) 1).doubleValue();
                double ff=((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"));

                total5pDist   += Math.abs(((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("5p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                total3pDist   += Math.abs(((Integer)hmIsomiRi.get("3p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("3p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                totalPolyDist += Math.abs(((Integer)hmIsomiRi.get("poly")).doubleValue() * ((Double)hmIsomiRi.get("fraction")) - ((Integer)hmIsomiRj.get("poly")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                totalPairwiseDist += Math.sqrt(
                  Math.pow(((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("5p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                + Math.pow(((Integer)hmIsomiRi.get("3p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("3p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                + Math.pow(((Integer)hmIsomiRi.get("poly")).doubleValue() * ((Double)hmIsomiRi.get("fraction")) - ((Integer)hmIsomiRj.get("poly")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                );
                
            }
            

        }

        int ii=0;
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

    /**
     * @return the totPairwiseDist
     */
    public double getTotPairwiseDist() {
        return totalPairwiseDist;
    }

    /**
     * @return the total3pDist
     */
    public double getTotal3pDist() {
        return total3pDist;
    }

    /**
     * @return the total5pDist
     */
    public double getTotal5pDist() {
        return total5pDist;
    }

    /**
     * @return the totalPolyDist
     */
    public double getTotalPolyDist() {
        return totalPolyDist;
    }
    
    
}
