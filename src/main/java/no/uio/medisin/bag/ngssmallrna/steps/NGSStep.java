/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author sr
 */
abstract public class NGSStep {
    
    protected StepInputData             stepInputData   = null;
    protected StepResultData            stepResultData  = null;
    
    protected               String      inFolder        = null;
    protected               String      outFolder       = null;
    
    protected static final  String      FILESEPARATOR   = System.getProperty("file.separator");

/*    
    abstract void       verifyInputData() throws IOException, NullPointerException;
    abstract void       verifyOutputData();
    abstract void       parseConfigurationData(HashMap configData) throws Exception;
    abstract HashMap    generateExampleConfigurationData();
    abstract void       execute() throws IOException;
*/    
    
    /**
     * set paths for input and output data folders
     * 
     */
    final void setPaths(){
        
        stepInputData.getInputFolder();
        String projectFolder = stepInputData.getProjectRoot() + FILESEPARATOR + stepInputData.getProjectID();
        projectFolder = projectFolder.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();
        
        inFolder = projectFolder + FILESEPARATOR + stepInputData.getInputFolder();
        outFolder = projectFolder + FILESEPARATOR + stepInputData.getOutputFolder();
        
    }
    
    /**
     * strip out duplicate folder delimiter from a file path
     * 
     * @param path
     * @return 
     */
    final String cleanPath(String path){
        return path.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
    }
    
}
