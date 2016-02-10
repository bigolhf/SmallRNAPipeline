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

public class StepUnzipInputFiles extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    public static final String          STEP_ID_STRING          = "UNZIP:";
    private static final String         ID_SOFTWARE             = "unzip_software:";    
    private static final String         ID_THREADS              = "no_of_threads:";
    
    private static final String         INFILE_EXTENSION         = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION        = ".fastq";
    //private static final String         inFolder                = "fastq_files";
    //private static final String         outFolder               = "fastq_files";
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";

    

    
    
    
    //private StepInputData               stepInputData;
    //private StepResultData              stepResultData;
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepUnzipInputFiles(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
        }
        catch(IOException exIO){
            
        }
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
        if(configData.get(ID_THREADS)==null) {
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        
        
        try{
            Integer.parseInt((String) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }
        
        if (Integer.parseInt((String) configData.get(ID_THREADS)) <= 0){
            throw new IllegalArgumentException(ID_THREADS + " <" + (String) configData.get(ID_THREADS) + "> must be positive");
        }
        setNoOfThreads(Integer.parseInt((String) configData.get(ID_THREADS)));

        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        /*
        
        unzipFastqParams.put("unzipSoftware",           this.getUnzipSoftware());
        unzipFastqParams.put("trimNoOfThreads",         this.getTrimNoOfThreads());
               
        pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz     

        */
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();                
                String outputFolder = pathToData + FileSeparator + outFolder;
                
                
                String fastqFile1 = outputFolder + FileSeparator + sampleData.getFastqFile1().replace(INFILE_EXTENSION, OUTFILE_EXTENSION);
                if(new File(fastqFile1).exists()){
                    logger.info("fastq file 1 <" + fastqFile1 + "> exists. Skipping");
                    continue;
                }
                                
                ArrayList<String> cmd1 = new ArrayList<>();
                cmd1.add(this.getUnzipSoftware());
                cmd1.add("-d");
                cmd1.add(fastqFile1);
                cmd1.add("-p " + this.getNoOfThreads());

                String cmdUnzip = StringUtils.join(cmd1, " ");
                cmdUnzip = cmdUnzip.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info(STEP_ID_STRING + " command:\t" + cmdUnzip);

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
                
                String fastqFile2 = outputFolder + FileSeparator + sampleData.getFastqFile2().replace(INFILE_EXTENSION, OUTFILE_EXTENSION);
                ArrayList<String> cmd2 = new ArrayList<>();
                cmd2.add(this.getUnzipSoftware());
                cmd2.add("-d");
                cmd2.add(fastqFile2);
                cmd2.add("-p " + this.getNoOfThreads());

                String cmdUnzip2 = StringUtils.join(cmd1, " ");
                cmdUnzip2 = cmdUnzip2.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info(STEP_ID_STRING + " command:\t" + cmdUnzip);

                Runtime rt2 = Runtime.getRuntime();
                Process proc2 = rt2.exec(cmdUnzip);
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
            catch(InterruptedException ex){
                logger.error("error executing pigz unzip command\n" + ex.toString());
            }
        }
        
        
    }
    
    
    
            
    /**
     * check the input data and parameters before we begin
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException, NullPointerException{
        
        logger.info("verify input data");
        
        // does unzip software exist?
        Validate.notNull((String) this.getUnzipSoftware());
        
        // is no of threads a positive integer?
        if (this.getNoOfThreads() <= 0)
        {
            logger.error("number of threads <" + this.getNoOfThreads() + "> must be positive");    
            return;            
        }
                    
        // check the data files
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = (String)sampleData.getFastqFile1();
            String fastqFile2 = (String)sampleData.getFastqFile2();
            
            //Fastq 1
            if (fastqFile1==null) throw new IOException("no Fastq1 file specified");
            
            if ((new File(fastqFile1)).exists()==false){
                throw new IOException("unzipFastqFiles: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException("unzipFastqFiles: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (fastqFile2==null) continue;
            
            if ((new File(fastqFile2)).exists()==false){
                throw new IOException("unzipFastqFiles: fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException("unzipFastqFiles: incorrect file extension for fastq file 2 <" 
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
        
        HashMap configData = new HashMap();
        HashMap paramData = new HashMap();
        
        paramData.put(ID_THREADS, 4);
        configData.put(STEP_ID_STRING, paramData);

        
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
