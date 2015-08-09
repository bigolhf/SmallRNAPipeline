/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 *
 * @author simonray
 */
public class IsomirPoint {
    private int fivePrimeSteps;
    int threePrimeSteps;
    int polySteps;
    int fraction;
            
    
    
    /**
     * report this isomiRPt
     * 
     * @return 
     */
    public String tabReportIsomirPt(){
        String reportStr = "";
        
        reportStr = reportStr.concat(String.valueOf(threePrimeSteps));
        
        return reportStr;
    }
}
