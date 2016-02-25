/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Adapter Trimming Step
 *  perform QC of FASTQ file list.
 * 
 *   Input is a FASTQ file
 *   Output is a quality report
 * 
 * 
 * @author sr
 */

public class StepFastqQC extends NGSStep implements NGSBase{
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final String          STEP_ID_STRING          = "FastqQC";
    private static final String         ID_SOFTWARE             = "unzipSoftware";    
    private static final String         ID_THREADS              = "noOfThreads";
    
    private static final String         INFILE_EXTENSION         = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION        = ".fastq";
    
    private int                         noOfSampleReads         = 100000;
    private String                      unzipSoftware           = "";

    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepFastqQC(StepInputData sid){
       stepInputData = sid;
    }
    
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        
        
        
        this.setUnzipSoftware((String) configData.get(ID_SOFTWARE));
 
        logger.info("passed");
    }
    
    
    
    /**
     * unzip specified file types
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /*
            sample unzip command       
            pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz     
        */
        logger.info(STEP_ID_STRING + ": execute step");        
        
        String fastqFile1in = "";
        String fastqFile1out = "";
        String fastqFile2in = "";
        String fastqFile2out = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();                
                
                fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
                fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
                
                if (sampleData.getFastqFile2() == null) continue;
                                
            }
            catch(Exception ex){
                logger.error("error executing pigz unzip command");
                throw new IOException("error executing pigz unzip command");
            }
        }
        
        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if(new File(this.getUnzipSoftware()).exists() == false){
            logger.error("unzip software not found at location < " + this.getUnzipSoftware() +">");
            throw new IOException("unzip software not found at location < " + this.getUnzipSoftware() +">");
        }
                
                
                            
        // check the data files
        String fastqFile1in = "";
        String fastqFile1out = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            
            fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
            fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            if (new File(fastqFile1in).exists()==false && new File(fastqFile1out).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";
            fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION));
            fastqFile2out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2());

            
            if ((new File(fastqFile2in)).exists()==false && (new File(fastqFile2in)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
            }
            if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
            
        }

    }
    
    
    
    @Override
    public void verifyOutputData(){
        logger.info("no output verification required");
    }
    
    
    
    /**
     * generate sample configuration data so the user can see what can be 
     * specified 
     * 
     * @return 
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");
        
        HashMap<String, Object> configData = new HashMap();
        
        configData.put(ID_SOFTWARE, "/usr/local/pigz");
        configData.put(ID_THREADS, 4);
        
        return configData;
        
    }

    
    
    
    
    
    /**
     * @return the unzipSoftware
     */
    public String getUnzipSoftware() {
        return unzipSoftware;
    }

    /**
     * @param unzipSoftware the unzipSoftware to set
     */
    public void setUnzipSoftware(String unzipSoftware) {
        this.unzipSoftware = unzipSoftware;
    }
    
    
    
    
    
}
