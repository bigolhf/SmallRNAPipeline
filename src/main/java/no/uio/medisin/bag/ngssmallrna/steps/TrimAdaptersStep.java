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
 *  Performs adapter trimming via a call to Trimmomatic.
 * 
 *   Input is a raw FASTQ file
 *   Output is a trimmed FASTQ file
 * 
 * 
 * @author sr
 */

public class TrimAdaptersStep extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private static final String         infileExtension     = ".fastq";
    private static final String         outfileExtension    = ".trim.fastq";
    private static final String         inFolder            = "fastq_files";
    private static final String         outFolder           = "adapter_trimmed";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public TrimAdaptersStep(StepInputData sid){
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
        
            trimAdapterParams.put("trimAdapterSoftware", this.getAdapterTrimmingSoftware());
            trimAdapterParams.put("trimAdapterFile", this.getTrimAdapterFile());
            trimAdapterParams.put("trimNoOfMismatches", this.getTrimNoOfMismatches());
            trimAdapterParams.put("trimMinAlignScore", this.getTrimMinAlignScore());
            trimAdapterParams.put("trimNoOfThreads", this.getTrimNoOfThreads());

        */
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add("java -jar");
                cmd.add((String) stepInputData.getStepParams().get("trimAdapterSoftware"));
                cmd.add("SE");
                cmd.add("-phred64");
    //          cmd.add("-trimlog " + pathToData + FileSeparator + sampleData.getDataFile() + ".trimlog");  // this will create huge logfiles. Disabled for now
                cmd.add("-threads " + stepInputData.getStepParams().get("trimNoOfThreads"));
                cmd.add(pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile());

                String outputFolder = pathToData + FileSeparator + outFolder;
                outputFolder = outputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
                Boolean f = new File(outputFolder).mkdir();       
                if (f) logger.info("created output folder <" + outputFolder + "> for results" );

                cmd.add(outputFolder + FileSeparator + sampleData.getDataFile().replace(infileExtension, outfileExtension));
                cmd.add("ILLUMINACLIP:" + stepInputData.getStepParams().get("trimAdapterFile") 
                    + ":" + stepInputData.getStepParams().get("trimNoOfMismatches")
                    + ":30"
                    + ":" + stepInputData.getStepParams().get("trimMinAlignScore")
                );
//                if (stepInputData.getStepParams().get("trimMinAvgReadQuality") != null){
//                    cmd.add("AVGQUAL:" + stepInputData.getStepParams().get("trimMinAvgReadQuality"));
//                }

                    /*
                    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar 
                    SE 
                    -phred64	
                    <input file fastq file>
                    <trimmed output fastq file> 	
                    ILLUMINACLIP:Trimmomatic-0.33/adapters/TruSeqE-SE.fa:2:30:7 2:30:7 (Mismatches:N/A:Alignment Score)
                    AVGQUAL (drop read if average quality is below this threshold)
                    -trimlog (write log)
                    -threads
                    */      

                String cmdTrimAdapters = StringUtils.join(cmd, " ");
                cmdTrimAdapters = cmdTrimAdapters.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Adapter Trim command:\t" + cmdTrimAdapters);

                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdTrimAdapters);
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
                logger.error("error executing AdapterTrimming command\n" + ex.toString());
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
