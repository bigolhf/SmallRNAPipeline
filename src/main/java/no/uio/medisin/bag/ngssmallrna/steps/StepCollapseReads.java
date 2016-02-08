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
 *  Collapse Read Step
 *  1. Convert FASTQ files in input list to an equivalent set of FASTA files.
 *  2. Count up duplicate reads and store this information in the header line 
 * 
 * @author sr
 */

public class StepCollapseReads extends NGSStep{
    
    static  Logger                      logger = LogManager.getLogger();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private static final String         infileExtension         = ".trim.fastq";
    private static final String         faOutputExtension       = ".trim.fasta";
    private static final String         clpfaOutputExtension    = ".trim.clp.fasta";
    
    
    
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepCollapseReads(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
        }
        catch(IOException exIO){
            
        }
    }
    
    
    
    
    
    
    @Override
    public void execute(){
        this.setPaths();
        
        /*
        
            collapseReadsParams.put("fastqTofasta", this.getFastq2fastaSoftware());
            collapseReadsParams.put("collapseFasta", this.getCollapseFastaSoftware());

        */
        String cmdFQ2FA = "";
        String cmdClp = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                Boolean f = new File(outFolder).mkdir(); 
                if (f) logger.info("created output folder <" + outFolder + "> for results" );

                
 
                /*
                    fastq_to_fasta -i 1000.fastq -o 1000.fasta -Q33
                */                      
                ArrayList<String> cmdQ2A = new ArrayList<>();
                cmdQ2A.add((String) stepInputData.getStepParams().get("fastqTofasta"));

                String fastqInputFile = inFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", infileExtension);
                fastqInputFile = fastqInputFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();                
                if (new File(fastqInputFile).exists() == false)
                {
                    throw new IllegalArgumentException("CollapseReads: fastq2fasta input file <" + fastqInputFile + ">. " + "does not exist");
                }
                cmdQ2A.add("-i");
                cmdQ2A.add(fastqInputFile);

                cmdQ2A.add("-o");

                String fastaOutputFile = outFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", faOutputExtension);
                cmdQ2A.add(fastaOutputFile);
                
                cmdQ2A.add("-Q33");


                cmdFQ2FA = StringUtils.join(cmdQ2A, " ");
                cmdFQ2FA = cmdFQ2FA.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Fastq to Fasta command:\t" + cmdFQ2FA);
                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdFQ2FA);
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
                    System.out.println("Process exitValue: " + exitVal);            

                brStdin.close();
                brStdErr.close();
            
            
            
                /*
                    fastx_collapse -i 1000.fastq -o 1000.fasta -Q33
                */                      
                ArrayList<String> cmd2 = new ArrayList<>();
                cmd2.add((String) stepInputData.getStepParams().get("collapseFasta"));

                String fastaInputFile = fastaOutputFile;
                if (new File(fastaInputFile).exists() == false)
                {
                    throw new IllegalArgumentException("CollapseReads: collapse reads input file <" + fastaInputFile + ">. " + "does not exist");
                }
                cmd2.add("-i");
                cmd2.add(fastaInputFile);

                cmd2.add("-o");

                String clpOutputFile = outFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", clpfaOutputExtension);
                cmd2.add(clpOutputFile);
                

                cmdClp = StringUtils.join(cmd2, " ");
                cmdClp = cmdClp.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Fastq to Fasta command:\t" + cmdClp);

                Runtime rtClp = Runtime.getRuntime();
                Process procClp = rtClp.exec(cmdClp);
                BufferedReader brStdinClp  = new BufferedReader(new InputStreamReader(procClp.getInputStream()));
                BufferedReader brStdErrClp = new BufferedReader(new InputStreamReader(procClp.getErrorStream()));
                
                line = null;
                logger.info("<OUTPUT>");
                while ( (line = brStdinClp.readLine()) != null)
                    logger.info(line);
                logger.info("</OUTPUT>");
                
                logger.info("<ERROR>");
                while ( (line = brStdErrClp.readLine()) != null)
                    logger.info(line);
                logger.info("</ERROR>");
                
                
                int exitValclp = procClp.waitFor();            
                System.out.println("Process exitValue: " + exitValclp);            
            
                brStdinClp.close();
                brStdErrClp.close();
            
            }
            catch(IOException|InterruptedException exIE){
                logger.error("error executing CollapseReads command\n");
                logger.error("CMD1 is " + cmdFQ2FA);
                logger.error("CMD2 is " + cmdClp);
                logger.error(exIE.toString());
            }
        }
        
        
    }
    

    
    /**
     * check the folder/files exist
     */
    @Override
    public void verifyInputData(){
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
        }
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
    }
    
    @Override
    public void outputResultData(){
        
    }
}
