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
 *  Collapse Read Step
 *  1. Convert FASTQ files in input list to an equivalent set of FASTA files.
 *  2. Count up duplicate reads and store this information in the header line 
 * 
 * @author sr
 */

public class StepCollapseReads extends NGSStep{
    
    static  Logger                      logger = LogManager.getLogger();
    static  String                      FileSeparator = System.getProperty("file.separator");

    public static final String          STEP_ID_STRING          = "COLLAPSE_READS:";
    private static final String         ID_Q2A_SOFTWARE         = "/usr/local/bin/fastq_to_fasta";   
    private static final String         ID_COLLAPSE_SOFTWARE    = "/usr/local/bin/fastx_collapser";   
    
    
    private static final String         INFILE_EXTENSION        = ".trim.fastq";
    private static final String         FASTA_OUTFILE_EXTENSION = ".trim.fasta";
    private static final String         CLP_OUTFILE_EXTENSION   = ".trim.clp.fasta";
    
    private String                      fastq2fasta_software    = "";
    private String                      collapseFastaSoftware   = "";
    
    
    

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
    
    
    
    
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
    
        if(configData.get(ID_Q2A_SOFTWARE)==null) {
            throw new NullPointerException("<" + configData.get(ID_Q2A_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_COLLAPSE_SOFTWARE)==null) {
            throw new NullPointerException("<" + configData.get(ID_COLLAPSE_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        
        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        this.setPaths();
        
        /*
        
            collapseReadsParams.put("fastqTofasta", this.getFastq2fastaSoftware());
            collapseReadsParams.put("collapseFasta", this.getCollapseFastaSoftware());

        */
        String cmdFQ2FA = "";
        String cmdClp = "";
        
        String fastqInputFile = "";
        String fastaOutputFile = "";
        String clpOutputFile = ""; 
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
                cmdQ2A.add(this.getFastq2fasta_software());

                fastqInputFile = inFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
                fastqInputFile = fastqInputFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();                
                cmdQ2A.add("-i");
                cmdQ2A.add(fastqInputFile);

                cmdQ2A.add("-o");

                fastaOutputFile = outFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", FASTA_OUTFILE_EXTENSION);
                cmdQ2A.add(fastaOutputFile);
                
                cmdQ2A.add("-Q33");


                cmdFQ2FA = StringUtils.join(cmdQ2A, " ");
                cmdFQ2FA = cmdFQ2FA.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Fastq2Fasta command:\t" + cmdFQ2FA);
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
            }
            catch(IOException|InterruptedException exIE){
                logger.error("error executing Fastq2Fasta command\n");
                logger.error("CMD is " + cmdFQ2FA);
                throw new IOException("error executing Fastq2Fasta command\n" + cmdFQ2FA);
            }
            
            
            
                /*
                    fastx_collapse -i 1000.fastq -o 1000.fasta -Q33
                */                      
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                ArrayList<String> cmd2 = new ArrayList<>();
                cmd2.add(this.getCollapseFastaSoftware());

                String fastaInputFile = fastaOutputFile;
                cmd2.add("-i");
                cmd2.add(fastaInputFile);

                cmd2.add("-o");

                clpOutputFile = outFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", CLP_OUTFILE_EXTENSION);
                cmd2.add(clpOutputFile);
                

                cmdClp = StringUtils.join(cmd2, " ");
                cmdClp = cmdClp.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Fastq to Fasta command:\t" + cmdClp);

                Runtime rtClp = Runtime.getRuntime();
                Process procClp = rtClp.exec(cmdClp);
                BufferedReader brStdinClp  = new BufferedReader(new InputStreamReader(procClp.getInputStream()));
                BufferedReader brStdErrClp = new BufferedReader(new InputStreamReader(procClp.getErrorStream()));
                
                String line = null;
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
                logger.error("CMD is " + cmdClp);
                throw new IOException("error executing CollapseReads command\n" + cmdClp);
             }
        }
        
        
    }
    

    
    /**
     * check data makes sense 
     */
    @Override
    public void verifyInputData() throws IOException, NullPointerException{
        
        logger.info("verify input data");
        
        // does software exist?
        Validate.notNull(this.getFastq2fasta_software());
        Validate.notNull(this.getCollapseFastaSoftware());
        
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
        
        paramData.put(ID_Q2A_SOFTWARE, "/usr/local/bin/fastq_to_fasta");
        paramData.put(ID_COLLAPSE_SOFTWARE, "/usr/local/bin/fastx_collapser");
        configData.put(STEP_ID_STRING, paramData);
        
        return configData;
    }   
    
    
    
    
    
    /**
     * @return the fastq2fasta_software
     */
    public String getFastq2fasta_software() {
        return fastq2fasta_software;
    }

    /**
     * @param fastq2fasta_software the fastq2fasta_software to set
     */
    public void setFastq2fasta_software(String fastq2fasta_software) {
        this.fastq2fasta_software = fastq2fasta_software;
    }

    /**
     * @return the collapseFastaSoftware
     */
    public String getCollapseFastaSoftware() {
        return collapseFastaSoftware;
    }

    /**
     * @param collapseFastaSoftware the collapseFastaSoftware to set
     */
    public void setCollapseFastaSoftware(String collapseFastaSoftware) {
        this.collapseFastaSoftware = collapseFastaSoftware;
    }
}
