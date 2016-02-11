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
import java.util.Collection;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  perform clean up
 *  this primarily involves zipping up files to save space
 *  because we don't know which steps were run or where the files are located
 *  we have to find files by searching the project folder
 *  thus, this will only work for linux command line, there is no attempt
 *  to make it work for Windows
 * 
 *  
 * 
 *   Input is a zipped FASTQ file
 *   Output is a unzipped FASTQ file
 * 
 * need to add the ability to define which types of files to compress. e.g., fastq, fasta, sam
 * 
 * @author sr
 */

public class StepCleanUp extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    public  static final String     STEP_ID_STRING                  = "CleanUp";
    
    
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private static final String         INFILE_EXTENSION     = ".fastq.gz";
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
    public StepCleanUp(StepInputData sid){
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
        String projectRoot = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        File directory = new File(projectRoot);
        ArrayList<String> fileTypes = (ArrayList<String>)stepInputData.getStepParams().get("fileTypes");
        for(String fileType: fileTypes){
//            Collection<File> c = FileUtils.listFiles(directory, new WildcardFileFilter(fileType), null);
            String fileTypesArray[] = new String[fileTypes.size()];
            Collection<File> c = FileUtils.listFiles(directory, fileTypes.toArray(fileTypesArray), true);
            Iterator itFL = c.iterator();
            for(File f: c){
                logger.info(f.toString());
                try{
    
                    ArrayList<String> cmd = new ArrayList<>();
                    cmd.add((String) stepInputData.getStepParams().get("zipSoftware"));
                    cmd.add(f.toString());
                    cmd.add("-p " + stepInputData.getStepParams().get("trimNoOfThreads"));

                        /*
                        pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz 
                        */      

                    String cmdZip = StringUtils.join(cmd, " ");
                    cmdZip = cmdZip.replace(FileSeparator + FileSeparator, FileSeparator);
                    logger.info("pigz zip command:\t" + cmdZip);

                    Runtime rt = Runtime.getRuntime();
                    Process proc = rt.exec(cmdZip);
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
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
                
                String outputFolder = pathToData + FileSeparator + outFolder;
                String outputFile = outputFolder + FileSeparator + sampleData.getFastqFile1().replace(INFILE_EXTENSION, outfileExtension);
                if(new File(outputFile).exists()){
                    logger.info("Output file <" + outputFile + "> exists. Skipping");
                    continue;
                }
                
                String inputFile = outputFolder + FileSeparator + sampleData.getFastqFile1();
                
                
                
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add((String) stepInputData.getStepParams().get("unzipSoftware"));
                cmd.add(inputFile);
                cmd.add("-threads " + stepInputData.getStepParams().get("trimNoOfThreads"));
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
            catch(IOException|InterruptedException ex){
                logger.error("error executing pigz unzip command\n" + ex.toString());
            }
        }
        
        
    }
    
    
    
            
    @Override
    public void verifyInputData(){
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
    }
    
    @Override
    public void verifyOutputData(){
        
    }
}
