/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.IOException;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  Adapter Trimming Step
 *  Performs adapter trimming via a call to Trimmomatic.
 * 
 *   Input is a raw FASTQ file
 *   Output is a trimmed FASTQ file
 * 
 * 
 * @author sr
 */

public class TrimAdapters extends NGSStep{
    
    static Logger logger = LogManager.getLogger();
    
    private StepInputData stepInputData;
    private StepResultData stepResultData;
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public TrimAdapters(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
            
        }
        catch(IOException exIO){
            
        }
    }
    
    public void performTrimming(){
        logger.info("trimming");
    }
            
    
}
