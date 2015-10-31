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
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MappedRead;
import no.uio.medisin.bag.ngssmallrna.pipeline.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
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
    
    
   
    private static final int            COVERAGE_SPAN = 200;
    
    static Logger                       logger                      = LogManager.getLogger();
    static String                       FileSeparator               = System.getProperty("file.separator");

    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         posAnalysisOutFolder           = "position_analysis";

    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         positionsExtension          = ".trim.clp.gen.pos.tsv";
    private static final String         featuresExtension           = ".trim.clp.gen.features.tsv";
    
    
    private StepInputData               stepInputData;
    
    private GFFSet                      gffSet                      = new GFFSet();    
    private ArrayList<MappedRead>       mappedReads                 = new ArrayList<>();
    private ArrayList<MappedRead>       filteredReads               = new ArrayList<>(); // generated from mappedReads
    
    int[] coverage5 = new int[COVERAGE_SPAN]; 
    int[] coverage3 = new int[COVERAGE_SPAN];
    
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
        analyzeSAMStartPositionsParams.put("longest_feature", this.getSamParseLongestFeature());                    
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
            if(new File(pathToAnnotation + "genes.gff3").exists())
                annotationFile = pathToAnnotation + "genes.gtf3";
        
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
        String posAnalysisOutputFolder = pathToData + FileSeparator + posAnalysisOutFolder;
        posAnalysisOutputFolder = posAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(posAnalysisOutputFolder).mkdir();       
        if (fA) logger.info("created output folder <" + posAnalysisOutputFolder + "> for results" );
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        String samLine = null;
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            String positionFile  = posAnalysisOutFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", positionsExtension);
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
                                samEntry.getqName(), samEntry.getStrand(), Integer.parseInt(samEntry.getrName().split("-")[1])));

                    }
                brSAM.close();
                logger.debug("read " + mappedReads.size() + " mapped entries");

            }
            catch(IOException smIO){
                logger.error("error parsing SAM file " + samInputFile);
                logger.error(smIO);
            }
        
            
            /**
             * 
             * condense overlapping reads into single features, i.e., those which have start/stop positions
             * that vary by <= ´bleed´. Start a new feature when the separation between current read and feature is
             * great than ´separation´.
             * 
             * if a read is going to be used to identify a feature, then it is probably sufficient for it to arise
             * in only one sample, because it won´t be identified in a differential expression analysis anyway
             * 
             * need to accommodate the possibility of two features in the same position but on different strands
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
                String currentStrand    = mappedRead.getStrand();
                
                /*
                    need to accommodate the possibility of two features in the
                    same position but on different strands
                */
                int currentStart5       = 100000000;
                int currentStop5        = -1;
                int currentStart3       = 100000000;
                int currentStop3        = -1;
                
                /*
                    using an int array might not be the best way to proceed, but try for now
                */
                int coverage5Start;
                int coverage3Start;

                        
                if(currentStrand.equals(GFFEntry.PLUSSTRAND)){
                    currentStart5       = mappedRead.getStartPos();
                    currentStop5        = mappedRead.getEndPos();
                    coverage5Start      = currentStart5 - COVERAGE_SPAN/2;
                    coverage3Start      = coverage5Start + COVERAGE_SPAN;
                }else{
                    currentStart3       = mappedRead.getStartPos();
                    currentStop3        = mappedRead.getEndPos();                    
                    coverage3Start      = currentStart3 + COVERAGE_SPAN;
                    coverage5Start      = coverage3Start - COVERAGE_SPAN/2;
                }
                
                
                
                String featureOutFile = pathToData + FileSeparator + posAnalysisOutFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", featuresExtension);
                int featureCount = 0;
                BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureOutFile)));
                Iterator itMR = mappedReads.iterator();
                while(itMR.hasNext()){
                    mappedRead = (MappedRead)itMR.next();
                    Boolean startNewFeature = false;
                    if (mappedRead.getChr().equals(currentChr)==false){
                        startNewFeature = true;
                    }else{
                        switch(currentStrand){
                            case GFFEntry.PLUSSTRAND:
                                if ( mappedRead.getStartPos() - currentStart5 > separation){
                                    startNewFeature = true;
                                }else{
                                    if(mappedRead.getStartPos() < currentStart5){
                                        currentStart5 = mappedRead.getStartPos();
                                    }
                                    if(mappedRead.getEndPos() > currentStop5){
                                        currentStop5 = mappedRead.getEndPos();
                                    }
                                    for(int i=currentStart5-coverage5Start; i<mappedRead.getEndPos()-coverage5Start; i++){
                                        coverage5[i] += mappedRead.getCount();
                                    }
                                }
                                
                            case GFFEntry.NEGSTRAND:
                                if ( mappedRead.getStartPos() - currentStart5 < separation){
                                    startNewFeature = true;
                                }else{
                                    if(mappedRead.getStartPos() > currentStart5){
                                        currentStart5 = mappedRead.getStartPos();
                                    }
                                    if(mappedRead.getEndPos() < currentStop5){
                                        currentStop5 = mappedRead.getEndPos();
                                    }
                                    for(int i=currentStart3-coverage3Start; i>mappedRead.getEndPos()-coverage3Start; i--){
                                        coverage3[i] += mappedRead.getCount();
                                    }
                                }
                            
                            case GFFEntry.UNKSTRAND:
                                logger.error("no Strand information for read, cannot process ");
                                logger.error(mappedRead.toString());
                                
                                
                        }
                        int longestFeature = (int) stepInputData.getStepParams().get("longest_feature");
                        if(startNewFeature){
                            bwFT.write(featureCount + "\t" + currentStart5 + "\t" + currentStop5 + "\t" 
                                + (currentStop5 - currentStart5 + 1) + "\t" + this.countCoverage5(currentStart5, currentStop5)
                                + "\t" + this.countDispersion5(currentStart5, currentStop5) + "\n");
                            // have to write out both strands
                            // 5' strand
                            // for now, add all features so I can understand 
                            // how they vary in terms of characteristics
                            // length, average reads, dispersion
                            if(currentStop5 - currentStart5 < longestFeature){
                                
                            }
                            featureCount++;
                        }
                        
                             
                    }
                    
                }
            }
            catch(IOException exIO){
                logger.error(exIO);
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
     * calculate average count across the mapped region for + Strand
     * 
     * @param start
     * @param stop
     * @return 
     */
    private double countCoverage5(int start, int stop){
        int countTotal = 0;
        for(int i=start; i<stop; i++){
            countTotal += coverage5[i];
        }
        return (double)countTotal/(double)(stop-start+1);
    }
    
    
    /**
     * calculate average count across the mapped region for - Strand
     * 
     * @param start
     * @param stop
     * @return 
     */
    private double countCoverage3(int start, int stop){
        int countTotal = 0;
        for(int i=start; i<stop; i++){
            countTotal += coverage3[i];
        }
        return (double)countTotal/(double)(stop-start+1);
    }
    
    
    /**
     * estimate dispersion of counts across the feature for + strand
     * 
     * @param start
     * @param stop
     * @return 
     */
    private double countDispersion5(int start, int stop){
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for(int i=start; i<stop; i++){
            stats.addValue(coverage5[i]);
        }
        return stats.getStandardDeviation();
    }
    
    
    
    /**
     * estimate dispersion of counts across the feature for - strand
     * 
     * @param start
     * @param stop
     * @return 
     */
    private double countDisperation3(int start, int stop){
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for(int i=start; i<stop; i++){
            stats.addValue(coverage3[i]);
        }
        return stats.getStandardDeviation();
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
