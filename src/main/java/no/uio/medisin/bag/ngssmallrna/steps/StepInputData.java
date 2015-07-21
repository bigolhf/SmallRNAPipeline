/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;
/**
 *  Specifies the information needed to locate the input data associated with a step
 * ProcessBuilder
 * htsjdk
 * @author sr
 */
public class StepInputData {
    
    static Logger logger = LogManager.getLogger();
    private HashMap stepParams;
    private String projectID;
    private String projectRoot;
    private ArrayList<SampleDataEntry> sampleData;
    
    public StepInputData(HashMap params, String pid, String pRoot, ArrayList<SampleDataEntry>sdata){
        
        stepParams  = params;
        projectID   = pid;
        projectRoot = pRoot;
        sampleData  = sdata;
        
    }
    
    
    /**
     * check the input data is good
     * @throws IOException 
     */
    public void verifyInputData() throws IOException
    {
        logger.info("verify input data");
        // check the project root exists
        // check the project folder exists within root
        
        for (Object param: getStepParams().keySet()){
            logger.info((String) param);
        }
        
        
        if (!(new File(projectRoot)).exists()){
            throw new IOException();
        }
        logger.info("Project root <" + projectRoot + "> exists");
        
        String projectFolder = projectRoot + System.getProperty("file.separator") + projectID;
        if (!(new File(projectFolder)).exists()){
            throw new IOException();
        }
        logger.info("Project folder <" + projectFolder + ">  exists");
        

        
    }
    /**
     * @return the projectID
     */
    public String getProjectID() {
        return projectID;
    }

    /**
     * @param projectID the projectID to set
     */
    public void setProjectID(String projectID) {
        this.projectID = projectID;
    }

    /**
     * @return the projectRoot
     */
    public String getProjectRoot() {
        return projectRoot;
    }

    /**
     * @param projectRoot the projectRoot to set
     */
    public void setProjectRoot(String projectRoot) {
        this.projectRoot = projectRoot;
    }

    /**
     * @return the stepParams
     */
    public HashMap getStepParams() {
        return stepParams;
    }

    /**
     * @return the sampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return sampleData;
    }
    
    
}
