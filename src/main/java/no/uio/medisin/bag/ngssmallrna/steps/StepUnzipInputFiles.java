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
import java.nio.file.Files;
import java.util.ArrayList;
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
    
    private static final String         infileExtension     = ".fastq.gz";
    private static final String         outfileExtension    = ".fastq";
    private static final String         inFolder            = "fastq_files";
    private static final String         outFolder           = "fastq_files";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

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
    
    @Override
    public void execute() throws IOException{
        /*
        
        unzipFastqParams.put("unzipSoftware",           this.getUnzipSoftware());
        unzipFastqParams.put("trimNoOfThreads",         this.getTrimNoOfThreads());
               
        */
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
                
                String outputFolder = pathToData + FileSeparator + outFolder;
                String fastqFile1 = outputFolder + FileSeparator + sampleData.getFastqFile1().replace(infileExtension, outfileExtension);
                String fastqFile2 = outputFolder + FileSeparator + sampleData.getFastqFile2().replace(infileExtension, outfileExtension);
                
                // need to add Fastq2 command
                if(new File(fastqFile1).exists()){
                    logger.info("fastq file 1 <" + fastqFile1 + "> exists. Skipping");
                    continue;
                }
                
                String inputFile = outputFolder + FileSeparator + sampleData.getFastqFile1();
                
                
                
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add((String) stepInputData.getStepParams().get("unzipSoftware"));
                cmd.add("-d");
                cmd.add(inputFile);
                cmd.add("-p " + stepInputData.getStepParams().get("trimNoOfThreads"));
                cmd.add(pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getFastqFile1());

                    /*
                    pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz 
                    */      

                String cmdUnzip = StringUtils.join(cmd, " ");
                cmdUnzip = cmdUnzip.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("pigz unzip command:\t" + cmdUnzip);

                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdUnzip);
                BufferedReader brStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                BufferedReader brStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
                
                    String line = null;
                    logger.info("<OUTPUT>");
                    while ( (line = brStdin.readLine()) != null)
                        logger.info(line);
                    logger.info("</OUTPUT>");

                    logger.info("<ERROR>");
                    while ( (line = brStdErr.readLine()) != null)
                        logger.info(line);
                    logger.info("</ERROR>");
                
                
                int exitVal = proc.waitFor();            
                logger.info("Process exitValue: " + exitVal);   
                
                brStdin.close();
                brStdErr.close();
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
        Validate.notNull((String) stepInputData.getStepParams().get("unzipSoftware"));
        
        // is no of threads a positive integer?
        try{
            Integer.parseInt((String) stepInputData.getStepParams().get("trimNoOfThreads"));
        }
        catch(NumberFormatException exNm){
            logger.error("number of threads <" + (String) stepInputData.getStepParams().get("trimNoOfThreads") + "> is not an integer");
            return;
        }
        
        if (Integer.parseInt((String) stepInputData.getStepParams().get("trimNoOfThreads")) <= 0)
        {
            logger.error("number of threads <" + (String) stepInputData.getStepParams().get("trimNoOfThreads") + "> must be positive");    
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
                  + sampleData.getFastqFile1() + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(infileExtension.toUpperCase())==false)
            {
                throw new IOException("unzipFastqFiles: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + infileExtension + "> as extension");
            }
            
            
            if (fastqFile2==null) continue;
            
            if ((new File(fastqFile2)).exists()==false){
                throw new IOException("unzipFastqFiles: fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(infileExtension.toUpperCase())==false)
            {
                throw new IOException("unzipFastqFiles: incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + infileExtension + "> as extension");
            }
                        
            
        }

    }
    
    
    
    @Override
    public void outputResultData(){
        
    }
}
