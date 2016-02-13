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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;
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

public class StepCleanUp extends NGSStep implements NGSBase{
    
    static Logger                       logger = LogManager.getLogger();
    public  static final String         STEP_ID_STRING          = "CleanUp";
    private static final String         ID_SOFTWARE             = "unzipSoftware";    
    private static final String         ID_THREADS              = "noOfThreads";
    private static final String         ID_FILE_TYPES           = "fileTypes";
    
    
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";
    private ArrayList<String>           fileTypes;
    
    
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepCleanUp(StepInputData sid){
        stepInputData = sid;
    }
    




    /**
     * This parses out the hashMap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        if(configData.get(ID_SOFTWARE)==null) {
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_FILE_TYPES)==null) {
            throw new NullPointerException("<" + configData.get(ID_FILE_TYPES) + "> : Missing Definition in Configuration File");
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
        
        try{
            this.setFileTypes((ArrayList<String> )configData.get(ID_FILE_TYPES));
        }
        catch(Exception ex){
            throw new IOException("couldn't cast " + configData.get(ID_FILE_TYPES) + "value to ArrayList");
        }

        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        String projectRoot = stepInputData.getProjectRoot() + FILESEPARATOR + stepInputData.getProjectID();
        File directory = new File(projectRoot);
        String cmdZip = "";
        for(String fileType: fileTypes){
            String fileTypesArray[] = new String[fileTypes.size()];
            Collection<File> c = FileUtils.listFiles(directory, fileTypes.toArray(fileTypesArray), true);
            Iterator itFL = c.iterator();
            for(File f: c){
                logger.info(f.toString());
                try{
    
                    ArrayList<String> cmd = new ArrayList<>();
                    cmd.add(this.getUnzipSoftware());
                    cmd.add(f.toString());
                    cmd.add("-p " + this.getNoOfThreads());

                    /*
                    pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz 
                    */      

                    cmdZip = StringUtils.join(cmd, " ");
                    cmdZip = cmdZip.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
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
                    logger.error("error executing pigz unzip command\n" + cmdZip);
                    logger.error(ex.toString());
                    throw new IOException(STEP_ID_STRING + ": error executing pigz unzip command" + cmdZip);
                }
             
            }            
        }
        
        
        
    }
    
    
    
            
    @Override
    public void verifyInputData(){
        
        logger.info("verify input data");        
        this.setPaths();
        
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
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

        configData.put(ID_SOFTWARE, "/usr/local/pigz");
        configData.put(ID_THREADS, 4);
        configData.put(ID_FILE_TYPES, new ArrayList<>(Arrays.asList("fastq", "fasta", "sam")));

        
        return configData;
    }
    
    
    

    @Override
    public void verifyOutputData(){
        
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

    /**
     * @return the fileTypes
     */
    public ArrayList<String> getFileTypes() {
        return fileTypes;
    }

    /**
     * @param fileTypes the fileTypes to set
     */
    public void setFileTypes(ArrayList<String> fileTypes) {
        this.fileTypes = fileTypes;
    }
}
