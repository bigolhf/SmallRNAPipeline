/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  parse SAM file to extract and process the miRNA reads
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class ParseSAMForMiRNAsStep extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         abundReadsOutFolder         = "bowtie_abundant_mapped";
    private static final String         genomeReadsOutFolder        = "bowtie_genome_mapped";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         fastqAbundantAlnExtension   = ".trim.clp.abun.fasta";
    private static final String         fastqAbundantUnAlnExtension = ".trim.clp.notabun.fasta";
    private static final String         samAbundantAlnExtension     = ".trim.clp.abun.sam";
    private static final String         fastqGenomeAlnExtension     = ".trim.clp.gen.sam";
    private static final String         fastqGenomeUnAlnExtension   = ".trim.clp.unmap.fasta";
    private static final String         samGenomeAlnExtension       = ".trim.clp.gen.sam";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public ParseSAMForMiRNAsStep(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
        }
        catch(IOException exIO){
            
        }
    }
    
    @Override
    public void execute(){
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();                
                String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
                String line = null;
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                    while((line=brSAM.readLine())!= null){
                        
                        if(line.split("\t")[1].equals("16")){

                            List<String> output = new ArrayList<>();
                            Matcher match = Pattern.compile("[0-9]+|[a-z]+|[A-Z]+").matcher(line.split("\t")[12].split(":")[2]);
                            String outputString = "";
                            while (match.find()) {
                                output.add(match.group());
                                outputString = outputString.concat(match.group() + "|");
                            }
                            logger.info(outputString);
                            
                        }
                    }
                brSAM.close();
                
                
            }
            catch(IOException ex){
                logger.error("error executing AdapterTrimming command\n" + ex.toString());
            }
        }
        
        
    }
    
    
            
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
