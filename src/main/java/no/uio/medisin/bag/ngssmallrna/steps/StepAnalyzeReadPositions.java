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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MappedRead;
import no.uio.medisin.bag.ngssmallrna.pipeline.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * this step performs additional analyses on a sample set of mapped reads in a SAM
 * file. It can look at features such as
 * 
 *      1. distribution of start and end positions
 *      2. distribution of mapped read length
 * 
 * 
 * @author sr
 */
public class StepAnalyzeReadPositions extends NGSStep{
    
    
    public static final int QNAME = 1;
    public static final int FLAG = 2;
    public static final int RNAME = 3;
    public static final int POS = 4;
    public static final int MAPQ = 5;
    public static final int CIGAR = 6;
    public static final int RNEXT = 7;
    public static final int PNEXT = 8;
    public static final int TLEN = 9;
    public static final int SEQ = 10;
    public static final int QUAL = 11;
    
    static Logger                       logger                      = LogManager.getLogger();
    static String                       FileSeparator               = System.getProperty("file.separator");

    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         posAnalysisFolder           = "position_analysis";

    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         positionsExtension          = ".trim.clp.gen.pos.tsv";
    
    
    private StepInputData               stepInputData;
    
    private GFFSet                      gffSet                      = new GFFSet();    
    private ArrayList<MappedRead>       mappedReads                 = new ArrayList<>();
    private ArrayList<MappedRead>       filteredReads               = new ArrayList<>(); // generated from mappedReads
   /*
    int[][] startPositions = new int[169100][9];
    int[][] readLengths = new int[169100][9];
    
    int[][] mrnaStartPositions = new int[4665][9];
    int[][] mrnaStopPositions  = new int[4665][9];
    
    int[] chrFeatureCount = new int[10];
    */
    public StepAnalyzeReadPositions()
    {
       
    }
    
    
    
    
    /**
     * parse out SAM file to retrieve start and stop information for each read
     * 
        analyzeSAMStartPositionsParams.put("bleed", this.getSamParseStartPosBleed());
        analyzeSAMStartPositionsParams.put("feature_types", this.getSamParseFeatureTypes());
        analyzeSAMStartPositionsParams.put("host", this.getBowtieMappingReferenceGenome());
        analyzeSAMStartPositionsParams.put("genomeReferenceGFFFile", this.getGenomeAnnotationGFF());
        analyzeSAMStartPositionsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
        analyzeSAMStartPositionsParams.put("bowtieReferenceGenome", this.getBowtieMappingReferenceGenome());
     *  
     * 
     * @param filename
     * @throws IOException 
     */
    @Override
    public void execute() 
    {
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }

        /**
            read genes.gtf or genes.gff3 file
        */
        String annotationFile = "";
        String pathToAnnotation = stepInputData.getStepParams().get("bowtieMapGenomeRootFolder") 
            + FileSeparator + stepInputData.getStepParams().get("bowtieReferenceGenome") + "/Annotation/Genes";
        if(new File(pathToAnnotation + "genes.gtf").exists())
            annotationFile = pathToAnnotation + "genes.gtf";
        else        
            if(new File(pathToAnnotation + "genes.gtf").exists())
                annotationFile = pathToAnnotation + "genes.gtf";
        
        try{
            if(annotationFile == null)
                throw new IOException("no annotation file was found for reference genome " 
                        + stepInputData.getStepParams().get("bowtieReferenceGenome"));
            gffSet.readGFF(annotationFile);
        }
        catch(IOException exIO){
            logger.error("Exception trying to read Annotation file ");
            logger.error(exIO);
        }
        
        /**
         *  read and parse the SAM files
         */
        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        String samLine = null;
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            String positionFile  = posAnalysisFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", positionsExtension);
            String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
            
            logger.info("sam input file is " + samInputFile);
            logger.info("results will be written to " + positionFile);
            
            samLine = null;
        
            int lineCount = 0;
            try{
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                    while ((samLine = brSAM.readLine()) != null) {

                        logger.debug(samLine);
                        SAMEntry e;
                        SAMEntry samEntry = new SAMEntry(samLine);
                        mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(), 
                                samEntry.getqName(), samEntry.getStrand()));

                    }
                brSAM.close();
                logger.debug("read " + mappedReads.size() + " mapped entries");

            }
            catch(IOException smIO){
                logger.error("error parsing SAM file " + samInputFile);
                logger.error(smIO);
            }
        
            
            /**
             * investigate the mapped reads
             * condense overlapping reads into single feature, i.e., those which have start/stop positions
             * that vary by <= bleed
             * 
             * if a read is going to be used to identify a feature, then it is probably sufficient for it to arise
             * in only one sample, because it won´t be identified in a differential expression analysis anyway
             * 
             */
            try{
                Collections.sort(mappedReads);
                
                int bleed = (int) stepInputData.getStepParams().get("bleed");
                int separation = (int) stepInputData.getStepParams().get("separation");
                
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;

                MappedRead mappedRead = ((MappedRead)mappedReads.get(0));
                String currentChr       = mappedRead.getChr();
                int currentStart        = mappedRead.getStartPos();
                int currentStop         = mappedRead.getEndPos();
                String currentStrand    = mappedRead.getStrand();
                
                Iterator itMR = mappedReads.iterator();
                while(itMR.hasNext()){
                    mappedRead = (MappedRead)itMR.next();
                    if ( mappedRead.getStartPos() - currentStart > separation){
                        // time to start a new feature
                    }
                }
            }
            catch(Exception ex){
                logger.error(ex);
            }
           
        
        }
    }
    
    public static void main(String args[]){
        StepAnalyzeReadPositions readPos = new StepAnalyzeReadPositions();
        try{
            //readPos.readGFF("/home/sr/research/sweden/msy2.mrna.gff");
            readPos.execute();
//            readPos.execute("/home/sr/research/sweden/2726_S10.trim.gen.PB1.sam");        
        }
        catch(Exception exIO){
            logger.info("error reading SAM file");
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
