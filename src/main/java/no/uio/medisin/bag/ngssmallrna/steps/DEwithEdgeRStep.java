/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.security.SecureRandom;
import java.util.Random;

import no.uio.medisin.bag.ngssmallrna.pipeline.IsomiRSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.MirBaseSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   1. parse SAM file to extract and process the miRNA reads to determine isomiR content
 *   2. merge the counts from each sample to generate a single count file 
 *      that can be used for differential expression analysis
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class DEwithEdgeRStep extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "mirna_isomir_analysis";
    private static final String         deAnalysisOutFolder         = "de_analysis";
    private              String         pathToDEAnalysisOutputFolder= "";
    private              String         rScriptFilename             = "";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         isomirSummaryExtension      = ".trim.clp.gen.iso_summary.tsv";
    private static final String         isomirPrettyExtension       = ".trim.clp.gen.iso_pretty.tsv";
    private static final String         miRCountsExtension          = ".trim.clp.gen.mircounts.tsv";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    
    private MirBaseSet                  miRBaseMiRNAList            = new MirBaseSet();
    private List<MiRNAFeature>          miRNAList                   = new ArrayList<>();
    private List<MiRNAFeature>          miRNAHitList;
    private ArrayList<IsomiRSet>        isomiRList;
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public DEwithEdgeRStep(StepInputData sid){
        stepInputData = sid;
    }
    
    @Override
    public void execute(){
        /*
            diffExpressionAnalysisParams.put("pvalue", this.getDiffExpressionPVal());
            diffExpressionAnalysisParams.put("host", this.getBowtieMappingReferenceGenome());
            diffExpressionAnalysisParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());
        */
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }
    
        
        
        /*
            1. read in all sample count files and merge
            2. output merged count file
            3. generate R script to perform DE using EdgeR
            4. process EdgeR output file 
        */
        try{
            miRBaseMiRNAList.loadMiRBaseData((String) stepInputData.getStepParams().get("host"), (String) stepInputData.getStepParams().get("miRBaseHostGFFFile"));
        }
        catch(IOException exIO){
            logger.error("error reading miRBase reference file <" + (String) stepInputData.getStepParams().get("miRBaseHostGFFFile") + ">\n" + exIO.toString());
        }
        
        
        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        String miRNAInputFolder = pathToData + FileSeparator + inFolder;
        pathToDEAnalysisOutputFolder = pathToData + FileSeparator + deAnalysisOutFolder;
        
        pathToDEAnalysisOutputFolder = pathToDEAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(pathToDEAnalysisOutputFolder).mkdir();       
        if (fA) logger.info("created output folder <" + pathToDEAnalysisOutputFolder + "> for results" );
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();

        
        logger.info("Merging Count Files");
        String[] countStrings = new String[miRBaseMiRNAList.getNumberOfEntries()];
        Arrays.fill(countStrings, "");
        itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String  miRCountsFile  = miRNAInputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", miRCountsExtension);
            miRCountsFile = miRCountsFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();
            try{
                int m=0;
                BufferedReader brmiRCounts  = new BufferedReader(new FileReader(new File(miRCountsFile)));
                    String countLine = "";
                    while((countLine = brmiRCounts.readLine()) != null){
                        countStrings[m] = countStrings[m].concat("\t" + countLine.split("\t")[1].trim() );
                        m++;
                    }
                brmiRCounts.close();
            }
            catch(IOException ex){
                logger.error("error reading count files for merging <" + miRCountsFile + "> \n" + ex.toString());
            }
        }
        
        logger.info("Writing merged count files");
        String mergedCountsFile      = pathToDEAnalysisOutputFolder + FileSeparator + stepInputData.getProjectID() + ".merged.mirna_counts.tsv";    
        try{
            BufferedWriter bwMc = new BufferedWriter(new FileWriter(new File(mergedCountsFile)));
            int m=0;
            for(MiRNAFeature miR: miRBaseMiRNAList.getMiRBaseMiRNAList()){
                bwMc.write(miR.getMimatID() + ":" + miR.getName() + countStrings[m] + "\n");
                m++;
            }
            
            bwMc.close();
        }
        catch(IOException exIO){
            logger.info("error writing merged counts File <" + mergedCountsFile + ">\n" + exIO);        
        }
        
    }
    
    
    /**
     * construct the R script to perform DE analysis for this dataset
     * 
     * need to generate the input files based on the specified groups
     * 
     */
    private void buildRScript(){

        BigInteger big = new BigInteger(130, new Random());
        rScriptFilename = pathToDEAnalysisOutputFolder + FileSeparator + new BigInteger(130, new SecureRandom()).toString(32) + ".R";
        rScriptFilename = rScriptFilename.replace(FileSeparator + FileSeparator, FileSeparator);
        
        ArrayList<String> cmd = new ArrayList<>();

        cmd.add("Rscript");
        cmd.add(rScriptFilename);
         
    }
    
    
    
    private void executeRScript(){
        

        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        ArrayList<String> cmd = new ArrayList<>();

        cmd.add("Rscript");
        cmd.add(rScriptFilename);


        String cmdRunRScript = StringUtils.join(cmd, " ");
        cmdRunRScript = cmdRunRScript.replace(FileSeparator + FileSeparator, FileSeparator);
        logger.info("Rscript command:\t" + cmdRunRScript);

        try{
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(cmdRunRScript);
            BufferedReader brAStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

                String line = null;
                logger.info("<OUTPUT>");
                while ((line = brAStdin.readLine()) != null)
                    logger.info(line);
                logger.info("</OUTPUT>");

                logger.info("<ERROR>");
                while ( (line = brAStdErr.readLine()) != null){
                    logger.info(line);
                }
                // need to parse the output from Bowtie to get the mapping summary
                logger.info("</ERROR>");

                int exitVal = proc.waitFor();            
                logger.info("Process exitValue: " + exitVal);   

            brAStdin.close();
            brAStdErr.close();        
        }
        catch(IOException exIO){
            logger.info("error executing RScript command\n" + exIO);
        }
        catch(InterruptedException exIE){
            logger.info("error executing RScript command\n" + exIE);            
        }
    }
    
    
    /**
     * Verify Input Data for parsing SAM file for miRNAs
     * 
     */        
    @Override
    public void verifyInputData(){
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            /*
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
            */
            
        }
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
    }
    
    
    
    
    @Override
    public void outputResultData(){
        
    }
}
