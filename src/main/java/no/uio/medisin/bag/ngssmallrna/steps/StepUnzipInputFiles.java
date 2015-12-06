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
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
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
    public void execute(){
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
                String outputFile = outputFolder + FileSeparator + sampleData.getDataFile().replace(infileExtension, outfileExtension);
                if(new File(outputFile).exists()){
                    logger.info("Output file <" + outputFile + "> exists. Skipping");
                    continue;
                }
                
                String inputFile = outputFolder + FileSeparator + sampleData.getDataFile();
                
                
                
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add((String) stepInputData.getStepParams().get("unzipSoftware"));
                cmd.add(inputFile);
                cmd.add("-threads " + stepInputData.getStepParams().get("trimNoOfThreads"));
                cmd.add(pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile());

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
            catch(IOException|InterruptedException ex){
                logger.error("error executing pigz unzip command\n" + ex.toString());
            }
        }
        
        
    }
    
    
    
            
    @Override
    public void verifyInputData(){
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            if (sampleData.getDataFile().toUpperCase().endsWith(infileExtension.toUpperCase())==false)
            {
                throw new IllegalArgumentException("AdapterTrimming: incorrect file extension for input file <" 
                  + sampleData.getDataFile() + ">. " 
                  + "should have <" + infileExtension + "> as extension");
            }
            
            if (sampleData.getDataFile().toUpperCase().endsWith(outfileExtension.toUpperCase())==true)
            {
                logger.warn("AdapterTrimming: input file has output file extension (.trim.fastq)");
                logger.warn("this file has already been trimmed");
            }
            
        }
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
    }
    
    @Override
    public void outputResultData(){
        
    }
}
