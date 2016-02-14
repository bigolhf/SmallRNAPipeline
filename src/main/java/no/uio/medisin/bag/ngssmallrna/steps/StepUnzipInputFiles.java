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
import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Adapter Trimming Step
 *  Unzip FASTQ files using pigz.
 * 
 *   Input is a zipped FASTQ file
 *   Output is a unzipped FASTQ file
 * 
 * 
 * @author sr
 */

public class StepUnzipInputFiles extends NGSStep implements NGSBase{
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final String          STEP_ID_STRING          = "UnzipInputFiles";
    private static final String         ID_SOFTWARE             = "unzipSoftware";    
    private static final String         ID_THREADS              = "noOfThreads";
    
    private static final String         INFILE_EXTENSION         = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION        = ".fastq";
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";

    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepUnzipInputFiles(StepInputData sid){
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
        
        
        try{
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }
        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
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
        String fastqFile1 = "";
        String fastqFile2 = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FILESEPARATOR + stepInputData.getProjectID();                
                String outputFolder = pathToData + FILESEPARATOR + outFolder;
                
                
                fastqFile1 = outputFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
                if(new File(this.cleanPath(fastqFile1)).exists()){
                                
                    ArrayList<String> cmd1 = new ArrayList<>();
                    cmd1.add(this.getUnzipSoftware());
                    cmd1.add("-d");
                    cmd1.add(fastqFile1);
                    cmd1.add("-p " + this.getNoOfThreads());

                    String cmdUnzip = StringUtils.join(cmd1, " ");
                    cmdUnzip = cmdUnzip.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
                    logger.info(STEP_ID_STRING + " fastq1 unzip command:\t" + cmdUnzip);

                    Runtime rt1 = Runtime.getRuntime();
                    Process proc1 = rt1.exec(cmdUnzip);
                    BufferedReader brStdin1  = new BufferedReader(new InputStreamReader(proc1.getInputStream()));
                    BufferedReader brStdErr1 = new BufferedReader(new InputStreamReader(proc1.getErrorStream()));

                        String line1 = null;
                        logger.info("<OUTPUT>");
                        while ( (line1 = brStdin1.readLine()) != null)
                            logger.info(line1);
                        logger.info("</OUTPUT>");

                        logger.info("<ERROR>");
                        while ( (line1 = brStdErr1.readLine()) != null)
                            logger.info(line1);
                        logger.info("</ERROR>");                

                    int exitVal = proc1.waitFor();            
                    logger.info("Process1 exitValue: " + exitVal);   

                    brStdin1.close();
                    brStdErr1.close();
                }
                else{               
                    logger.info("fastq file 1 <" + fastqFile1 + "> exists. Skipping");
                }

                fastqFile2 = outputFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION);
                if(new File(this.cleanPath(fastqFile2)).exists()){
            
                    ArrayList<String> cmd2 = new ArrayList<>();
                    cmd2.add(this.getUnzipSoftware());
                    cmd2.add("-d");
                    cmd2.add(fastqFile2);
                    cmd2.add("-p " + this.getNoOfThreads());

                    String cmdUnzip2 = StringUtils.join(cmd2, " ");
                    cmdUnzip2 = cmdUnzip2.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
                    logger.info(STEP_ID_STRING + " fastq2 unzip command:\t" + cmdUnzip2);

                    Runtime rt2 = Runtime.getRuntime();
                    Process proc2 = rt2.exec(cmdUnzip2);
                    BufferedReader brStdin2  = new BufferedReader(new InputStreamReader(proc2.getInputStream()));
                    BufferedReader brStdErr2 = new BufferedReader(new InputStreamReader(proc2.getErrorStream()));

                        String line2 = null;
                        logger.info("<OUTPUT>");
                        while ( (line2 = brStdin2.readLine()) != null)
                            logger.info(line2);
                        logger.info("</OUTPUT>");

                        logger.info("<ERROR>");
                        while ( (line2 = brStdErr2.readLine()) != null)
                            logger.info(line2);
                        logger.info("</ERROR>");                

                    int exitVal2 = proc2.waitFor();            
                    logger.info("Process2 exitValue: " + exitVal2);   

                    brStdin2.close();
                    brStdErr2.close();
                }
                else{               
                    logger.info("fastq file 2 <" + fastqFile2 + "> exists. Skipping");
                }
                                
            }
            catch(InterruptedException ex){
                logger.error("error executing pigz unzip command");
                throw new IOException("error executing pigz unzip command");
            }
        }
        
        
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
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            
            //Fastq 1
            if (fastqFile1==null) throw new IOException("no Fastq1 file specified");
            
            if (new File(this.cleanPath(fastqFile1)).exists()==false){
                logger.info(STEP_ID_STRING + ": fastq File1 <" + fastqFile1 + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq File1 <" + fastqFile1 + "> does not exist");
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
            String fastqFile2 = inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION);

            
            if ((new File(fastqFile2)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq File2 <" + fastqFile2 + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
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
     * @return the noOfThreads
     */
    public int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * @param noOfThreads the noOfThreads to set
     */
    public void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads;
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
