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

public class DStepPairedTrimAdapters extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private static final String         infileExtension     = ".fastq";
    private static final String         outfileExtension    = ".trim.fastq";
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public DStepPairedTrimAdapters(StepInputData sid){
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
        
        this.setPaths();
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add("java -jar");
                cmd.add((String) stepInputData.getStepParams().get("trimAdapterSoftware"));
                cmd.add("SE");
                cmd.add("-phred64");
    //          cmd.add("-trimlog " + pathToData + FileSeparator + sampleData.getDataFile() + ".trimlog");  // this will create huge logfiles. Disabled for now
                cmd.add("-threads " + stepInputData.getStepParams().get("trimNoOfThreads"));
                cmd.add(inFolder + FileSeparator + sampleData.getFastqFile1());

                Boolean f = new File(outFolder).mkdir();       
                if (f) logger.info("created output folder <" + outFolder + "> for results" );

                cmd.add(outFolder + FileSeparator + sampleData.getFastqFile1().replace(infileExtension, outfileExtension));
                cmd.add("ILLUMINACLIP:" + stepInputData.getStepParams().get("trimAdapterFile") 
                    + ":" + stepInputData.getStepParams().get("trimNoOfMismatches")
                    + ":30"
                    + ":" + stepInputData.getStepParams().get("trimMinAlignScore")
                );
//                if (stepInputData.getStepParams().get("trimMinAvgReadQuality") != null){
//                    cmd.add("AVGQUAL:" + stepInputData.getStepParams().get("trimMinAvgReadQuality"));
//                }

                /*

                java -jar /data/projects/simonray/software/Trimmomatic-0.33/trimmomatic-0.33.jar 
                PE 
                <input file fastq file 1>
                <input file fastq file 2>
                <output file fastq file 1 paired reads>
                <output file fastq file 1 unpaired reads>
                <output file fastq file 2 paired reads>
                <output file fastq file 2 unpaired reads>
                ILLUMINACLIP:<path to adapter sequence file>:<MISMATCHES>:<MATCH_SCORE>:10 LEADING:3 TRAILING:3 MINLEN:30
                AVGQUAL (drop read if average quality is below this threshold)
                -trimlog (write log)
                -threads 16 

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
    public void verifyInputData() throws IOException{
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = (String)sampleData.getFastqFile1();
            String fastqFile2 = (String)sampleData.getFastqFile2();
            
            if (fastqFile1==null) throw new IOException("no Fastq1 file specified");
            
            if ((new File(fastqFile1)).exists()==false){
                throw new IOException("AdapterTrimming: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(infileExtension.toUpperCase())==false){
                throw new IOException("AdapterTrimming: incorrect file extension for input file <" 
                  + fastqFile1 + ">. " 
                  + "should have <" + infileExtension + "> as extension");
            }
            
            //Fastq 2
            if (fastqFile2==null) continue;
            
            if ((new File(fastqFile2)).exists()==false){
                throw new IOException("AdapterTrimming: fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(infileExtension.toUpperCase())==false)
            {
                throw new IOException("AdapterTrimming: incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + infileExtension + "> as extension");
            }
        }

    }
    
    
    
    @Override
    public void verifyOutputData(){
        
    }
}
