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
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.jmirpara.SimpleSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.GenomeSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.MappedRead;
import no.uio.medisin.bag.ngssmallrna.pipeline.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.Strand;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * this step performs additional analyses on a sample set of mapped reads in a
 * SAM file. It can look at features such as
 *
 * 1. distribution of start and end positions 2. distribution of mapped read
 * length
 *
 *
 * @author sr
 */
public class StepAnalyzeSAMforStartPositions extends NGSStep implements NGSBase{
    
    public  static final String     STEP_ID_STRING      = "AnalyzeStartPositions";
    private static final String     ID_BLEED            = "bleed";
    private static final String     ID_BASELINE         = "baselinePercent";
    private static final String     ID_ISOMIRS          = "analyzeIsomirs";
    private static final String     ID_MIRBASE_VERSION  = "mirbaseVersion";
    private static final String     ID_REF_GENOME       = "host";
    private static final String     ID_SHORTEST_FEATURE = "shortestFeature";
    private static final String     ID_LONGEST_FEATURE  = "longestFeature";
    private static final String     ID_MIN_COUNTS       = "minCounts";
    private static final String     ID_SEPARATION       = "featureSeparation";
    
    private static final String     INFILE_EXTENSION    = ".trim.clp.gen.sam";
    private static final String     POS_FILE_EXT        = ".trim.clp.gen.pos.tsv";
    private static final String     FEAT_FILE_EXT       = ".trim.clp.gen.features.tsv";

    private static final int        COVERAGE_SPAN       = 200;
    
    static Logger logger = LogManager.getLogger();
    
    
    private int                     locationBleed       = 2;
    private int                     baselinePercent     = 5;
    private Boolean                 analyzeIsomirs      = false;
    private int                     miRBaseRelease      = 20;
    private String                  ReferenceGenome     = "";
    private int                     shortestRead        = 0;
    private int                     longestRead         = 0;
    private int                     minCounts           = 0;
    private int                     separation          = 0;

    
    private GenomeSeq               genomeFasta;
    private GFFSet                  gffSet              = new GFFSet();    
    private ArrayList<MappedRead>   mappedReads         = new ArrayList<>();
    
    int[]                           coverage5           = new int[COVERAGE_SPAN];    
    int[]                           coverage3           = new int[COVERAGE_SPAN];
    
    private GFFSet                  featureSet          = new GFFSet(); // stores the identified features
    

    public StepAnalyzeSAMforStartPositions(StepInputData sid) {
        stepInputData = sid;
    }

    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_BLEED)==null) {
            throw new NullPointerException("<" + ID_BLEED + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_SHORTEST_FEATURE)==null) {
            throw new NullPointerException("<" + ID_SHORTEST_FEATURE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_LONGEST_FEATURE)==null) {
            throw new NullPointerException("<" + ID_LONGEST_FEATURE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIN_COUNTS)==null) {
            throw new NullPointerException("<" + ID_MIN_COUNTS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_SEPARATION)==null) {
            throw new NullPointerException("<" + ID_SEPARATION + "> : Missing Definition in Configuration File");
        }
        

      
        
        try{
            Integer.parseInt((String) configData.get(ID_BASELINE));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_BASELINE + " <" + configData.get(ID_BASELINE) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_BASELINE)) <= 0){
            throw new IllegalArgumentException(ID_BASELINE + " <" + (String) configData.get(ID_BASELINE) + "> must be an integer between 0 and 100");
        }
        this.setBaselinePercent(Integer.parseInt((String) configData.get(ID_BASELINE)));

        
        try{
            Integer.parseInt((String) configData.get(ID_SHORTEST_FEATURE));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_SHORTEST_FEATURE + " <" + ID_SHORTEST_FEATURE + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_SHORTEST_FEATURE)) <= 0 ){
            throw new IllegalArgumentException(ID_SHORTEST_FEATURE + " <" + (String) configData.get(ID_SHORTEST_FEATURE) + "> must be > 0 ");
        }
        this.setShortestRead(Integer.parseInt((String) configData.get(ID_SHORTEST_FEATURE)));
        

        try{
            Integer.parseInt((String) configData.get(ID_LONGEST_FEATURE));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_LONGEST_FEATURE + " <" + ID_LONGEST_FEATURE + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_LONGEST_FEATURE)) <= 0 ){
            throw new IllegalArgumentException(ID_LONGEST_FEATURE + " <" + (String) configData.get(ID_LONGEST_FEATURE) + "> must be > 0 ");
        }
        this.setLongestRead(Integer.parseInt((String) configData.get(ID_BLEED)));
        

        try{
            Integer.parseInt((String) configData.get(ID_MIN_COUNTS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MIN_COUNTS + " <" + ID_MIN_COUNTS + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MIN_COUNTS)) <= 0 ){
            throw new IllegalArgumentException(ID_MIN_COUNTS + " <" + (String) configData.get(ID_MIN_COUNTS) + "> must be > 0 ");
        }
        this.setMinCounts(Integer.parseInt((String) configData.get(ID_MIN_COUNTS)));
        

        try{
            Integer.parseInt((String) configData.get(ID_BLEED));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_BLEED + " <" + ID_BLEED + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_BLEED)) <= 0 ){
            throw new IllegalArgumentException(ID_BLEED + " <" + (String) configData.get(ID_BLEED) + "> must be > 0 ");
        }
        this.setLocationBleed(Integer.parseInt((String) configData.get(ID_BLEED)));
        

        try{
            Integer.parseInt((String) configData.get(ID_SEPARATION));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_SEPARATION + " <" + ID_SEPARATION + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_SEPARATION)) <= 0 ){
            throw new IllegalArgumentException(ID_SEPARATION + " <" + (String) configData.get(ID_SEPARATION) + "> must be > 0 ");
        }
        this.setLocationBleed(Integer.parseInt((String) configData.get(ID_BLEED)));
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + (String) configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        try{
            Boolean.parseBoolean((String) configData.get(ID_ISOMIRS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_BLEED + " <" + (String) configData.get(ID_BLEED) + "> cannot be cast as Boolean");
        }        
        this.setAnalyzeIsomirs(Boolean.parseBoolean((String) configData.get(ID_ISOMIRS)));
        

        logger.info("passed");
    }
    
    
    
    
    
    
    /**
     * parse out SAM file to retrieve start and stop information for each read
     *
     *
     * @param filename
     * @throws IOException
     */
    @Override
    public void execute() throws IOException{
        
        stepInputData.verifyInputData();            
        this.setPaths();
    
        
        /**
         * read genome fasta
         */
        String hostCode = this.getReferenceGenome();
        genomeFasta = new GenomeSeq(hostCode);
        String pathToFasta = stepInputData.getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + hostCode + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH;
        String genomeFastaFile = pathToFasta + FILESEPARATOR + "genome.fa";
        genomeFastaFile = genomeFastaFile.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();
        try{
            logger.info("reading genome file + <" + genomeFastaFile + ">");
            this.genomeFasta.readFastaGenome(genomeFastaFile);
            logger.info("finished ");
            logger.info("read " + genomeFasta.getNoOfBases() + " bases");
            logger.info("spanning " + genomeFasta.getNoOfChr() + " chromosomes");
        }
        catch(IOException exIO){
            logger.error("exception reading Genome reference file + <" + genomeFastaFile + ">");
            logger.error(exIO.toString());
            throw new IOException(STEP_ID_STRING + ": exception reading Genome reference file + <" + genomeFastaFile + ">");
        }
        /**
         * read genes.gtf or genes.gff3 file
         */
        String annotationFile = "";
        String pathToAnnotation = stepInputData.getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + hostCode + ReferenceDataLocations.ID_GENE_ANNOTATION;
        File f = new File(pathToAnnotation + FILESEPARATOR + "genes.gtf");
        if (new File(pathToAnnotation + FILESEPARATOR + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gtf";
        } else if (new File(pathToAnnotation + FILESEPARATOR + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gff";
        }
        
        try {
            if (annotationFile == null) {
                throw new IOException("no annotation file was found for reference genome "
                        + hostCode);
            }
            gffSet.readGFF(annotationFile);
        } catch (IOException exIO) {
            logger.error("Exception trying to read Annotation file ");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": exception reading Genome reference file + <" + genomeFastaFile + ">");
        }

        /**
         * read and parse the SAM files
         */
        
        int shortestFeature = this.getShortestRead();
        int longestFeature  = this.getLongestRead();
        
//        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
//        pathToData = pathToData.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        
//        String outputFolder = pathToData + FileSeparator + stepInputData.getOutputFolder();
//        String posAnalysisOutputFolder = pathToData + FileSeparator + posAnalysisOutFolder;
//        outputFolder = outputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(outFolder).mkdir();        
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        Iterator itSD = this.stepInputData.getSampleData().iterator();
                int featureCount = 0;
        String samLine = null;
        while (itSD.hasNext()) {
            SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
            
            String positionFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", POS_FILE_EXT);
            String featureOutFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", FEAT_FILE_EXT);
            featureOutFile = featureOutFile.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();

            String samInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            
            logger.info("sam input file is " + samInputFile);
            logger.info("results will be written to " + positionFile);
            
            samLine = null;
            
            int lineCount = 0;
            try {
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                while ((samLine = brSAM.readLine()) != null) {
                    
                    SAMEntry e;
                    SAMEntry samEntry = new SAMEntry(samLine);
                    if (samEntry.isHeaderLine() == false && samEntry.isMappedRead()) {
                        mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(),
                                samEntry.getrName(), samEntry.getStrand(), Integer.parseInt(samEntry.getqName().split("-")[1])));
                    }
                    
                }
                brSAM.close();
                logger.debug("read " + mappedReads.size() + " mapped entries");
                
            } catch (IOException smIO) {
                logger.error("error parsing SAM file " + samInputFile);
                logger.error(samLine);
                logger.error(smIO);
                throw new IOException(STEP_ID_STRING + ": error parsing SAM file <" + samInputFile + ">");
            }

            /**
             *
             * condense overlapping reads into single features, i.e., those
             * which have start/stop positions that vary by <= ´bleed´. Start a
             * new feature when the separation between current read and feature
             * is great than ´separation´.
             *
             * if a read is going to be used to identify a feature, then it is
             * probably sufficient for it to arise in only one sample, because
             * it won´t be identified in a differential expression analysis
             * anyway
             *
             * need to accommodate the possibility of two features in the same
             * position but on different strands therefore we sort by strand and
             * position
             */
            try {
                Collections.sort(mappedReads);

                int bleed = this.getLocationBleed();
                
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;
                
                MappedRead mappedRead = ((MappedRead) mappedReads.get(0));
                //logger.debug(mappedRead.toString());
                String currentChr = mappedRead.getChr();
                Strand currentStrand = mappedRead.getStrand();

                /*
                 need to accommodate the possibility of two features in the
                 same position but on different strands
                 */
                int currentStart5 = 100000000;
                int currentStop5 = -1;
                int currentStart3 = -1;
                int currentStop3 = 100000000;

                /*
                 using an int array might not be the best way to proceed, but try for now
                 */
                int coverage5Start = -1;
                int coverage3Start = -1;
                
                if (currentStrand.equals(Strand.PLUS)) {
                    currentStart5 = mappedRead.getStartPos();
                    currentStop5 = mappedRead.getEndPos();
                    coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                } else {
                    currentStart3 = mappedRead.getStartPos();
                    currentStop3 = mappedRead.getEndPos();                    
                    coverage3Start = currentStart3 - COVERAGE_SPAN / 2;
                }
                
                //String featureOutFile = pathToData + FileSeparator + posAnalysisOutFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", featuresExtension);
//                positionFile
                //featureOutFile = featureOutFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();x
                
                BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureOutFile)));
                bwFT.write("ID\tChr\t \tStart\tStop\tLength\tCoverage\tDispersion\n");

                Iterator itMR = mappedReads.iterator();
                Boolean startNewFeature3 = false;
                Boolean startNewFeature5 = false;
                while (itMR.hasNext()) {
                    mappedRead = (MappedRead) itMR.next();
                    //logger.info(mappedRead.toString());
                    if (mappedRead.getChr().equals(currentChr) == false || mappedRead.getStrand().equals(currentStrand) == false) {
                        if (mappedRead.getChr().equals(currentChr) == false) {
                            //logger.debug("different chr. start new feature");
                        }
                        if (mappedRead.getStrand().equals(currentStrand) == false) {
                            logger.debug("different strand. start new feature");                            
                        }
                        currentStrand = mappedRead.getStrand();
                        currentChr = mappedRead.getChr();
                        if (currentStrand == Strand.PLUS) {
                            startNewFeature3 = false;
                            startNewFeature5 = false;
                            currentStart5 = mappedRead.getStartPos();
                            currentStop5 = mappedRead.getEndPos();
                            coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                            addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                        } else {
                            startNewFeature3 = false;
                            startNewFeature5 = false;
                            currentStart3 = mappedRead.getStartPos();
                            currentStop3 = mappedRead.getEndPos();
                            coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                            addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                        }
                        
                    } else {
                        switch (currentStrand) {
                            case PLUS:
                                if (mappedRead.getStartPos() - currentStart5 > separation) {
                                    //logger.debug("read too far away. start new 5' feature");                                    
                                    startNewFeature5 = true;
                                } else {
                                    if (mappedRead.getStartPos() < currentStart5) {
                                        //logger.debug("extend 5´ start from " + currentStart5 + " to " + mappedRead.getStartPos());
                                        currentStart5 = mappedRead.getStartPos();
                                        if( currentStart5 <= coverage5Start)
                                            startNewFeature5 = true;
                                    }
                                    if (mappedRead.getEndPos() > currentStop5) {
                                        //logger.debug("extend 5´ end from " + currentStop5 + " to " + mappedRead.getEndPos());
                                        currentStop5 = mappedRead.getEndPos();
                                        if( currentStop5 >= coverage5Start + COVERAGE_SPAN )
                                            startNewFeature5 = true;
                                    }
                                    try {
                                        //logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                        this.addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                                        //logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                    } catch (ArrayIndexOutOfBoundsException exAB) {
                                        logger.error("5' Array out of bounds");
                                    }
                                }
                                break;
                            
                                
                            case MINUS:
                                if ((mappedRead.getStartPos() - currentStart3 > separation) || (currentStop3 - mappedRead.getEndPos() > separation)) {
                                    //logger.debug("read too far away. start new 3' feature");
                                    startNewFeature3 = true;
                                } else {
                                    //logger.debug("add read to existing feature");
                                    if (mappedRead.getStartPos() > currentStart3) {
                                        //logger.debug("extend 3´ start from " + currentStart3 + " to " + mappedRead.getStartPos());
                                        currentStart3 = mappedRead.getStartPos();
                                        if(currentStart3 <= coverage3Start)
                                            startNewFeature3 = true;
                                    }
                                    if (mappedRead.getEndPos() < currentStop3) {
                                        //logger.debug("extend 3´ end from " + currentStop3 + " to " + mappedRead.getEndPos());
                                        currentStop3 = mappedRead.getEndPos();
                                        if(currentStop3 >= coverage3Start + COVERAGE_SPAN) 
                                            startNewFeature3 = true;
                                    }
                                    try {
                                        //logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                        this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                                        //logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                    } catch (ArrayIndexOutOfBoundsException exAB) {
                                        logger.error("3' Array out of bounds");
                                        logger.error(mappedRead.toString());
                                        logger.error(exAB);
                                    }
                                }
                                break;
                            
                                
                            case UNKNOWN:
                                logger.error("no Strand information for read, cannot process ");
                                logger.error(mappedRead.toString());
                                break;
                            
                        }
                        
                    }
                    
                    if (startNewFeature5) {
                        if(currentStop5 - currentStart5 + 1 <= longestFeature 
                                && currentStop5 - currentStart5 + 1>= shortestFeature 
                                && countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start) >= minCounts){
                            bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + currentStart5 + "\t" + currentStop5 + "\t"
                                    + (currentStop5 - currentStart5 + 1) + "\t" + this.countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start)
                                    + "\t" + this.countDispersion5(currentStart5 - coverage5Start, currentStop5 - coverage5Start) + "\n");
                            if(featureSet.doesRegionContainFeature(currentStart5, currentStop5, currentStrand, currentChr, bleed)==false){
                                GFFEntry newEntry = new GFFEntry(hostCode  + "-" + Integer.toString(featureCount) + ":5|" + currentChr  + ":" + currentStart5+ "-" + currentStop5, currentStrand.toString(), currentChr, currentStart5, currentStop5);
                                String featureSeq = genomeFasta.getSubSeq(currentChr, currentStrand, currentStart5, currentStop5);
                                newEntry.setSequence(featureSeq);
                                featureSet.addEntry(newEntry);
                            }                            
                        }
                        
                        currentStart5 = mappedRead.getStartPos();
                        currentStop5 = mappedRead.getEndPos();

                        coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                        clearCoverage5();
                        try {
                            this.addCounts5(mappedRead.getStartPos() - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                        } catch (ArrayIndexOutOfBoundsException exAB) {
                            logger.error("5' Array out of bounds when writing data");
                            throw new IOException(STEP_ID_STRING + ": 5' Array out of bounds when writing data");
                        }
                        
                        currentStrand = mappedRead.getStrand();
                        currentChr = mappedRead.getChr();
                        featureCount++;
                        startNewFeature5 = false;
                        
                    } else {                                         
                        if (startNewFeature3 ) {
                            if(currentStop3 - currentStart3 + 1<= longestFeature 
                                && currentStop3 - currentStart3 + 1>= shortestFeature
                                && countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3) > minCounts){
                                bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + currentStart3 + "\t" + currentStop3 + "\t"
                                        + (currentStop3 - currentStart3 + 1) + "\t" + this.countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3)
                                        + "\t" + this.countDispersion3(coverage3Start - currentStop3, coverage3Start - currentStart3) + "\n");
                                if(featureSet.doesRegionContainFeature(currentStart3, currentStop3, currentStrand, currentChr, bleed)==false){
                                    GFFEntry newEntry = new GFFEntry(hostCode + "-" + Integer.toString(featureCount) + ":3|" + currentChr  + ":" + currentStart5+ "-" + currentStop5, currentStrand.toString(), currentChr, currentStart3, currentStop3);
                                    String featureSeq = SimpleSeq.complement(genomeFasta.getSubSeq(currentChr, currentStrand, currentStart3, currentStop3));
                                    newEntry.setSequence(featureSeq);
                                    featureSet.addEntry(newEntry);                                    
                                }
                                
                            }
                            
                            currentStart3 = mappedRead.getStartPos();
                            currentStop3 = mappedRead.getEndPos();

                            coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                            clearCoverage3();
                            try {
                                this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                            } catch (ArrayIndexOutOfBoundsException exAB) {
                                logger.error("3' Array out of bounds when writing data");
                                throw new IOException(STEP_ID_STRING + ": 3' Array out of bounds when writing data");
                            }
                            currentStrand = mappedRead.getStrand();
                            currentChr = mappedRead.getChr();                           
                            featureCount++;
                            startNewFeature3 = false;
                        }
                    }
                    
                }
                bwFT.close();
            } catch (IOException exIO) {
                logger.error(exIO);
                throw new IOException(STEP_ID_STRING + ": error writing feature details file <" + featureOutFile + ">");
            }
            
        }
        
        String featureFile = outFolder + FILESEPARATOR + stepInputData.getProjectID() + ".features.tsv";
        try{
            BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureFile)));
                featureSet.writeFeaturesAsGFF3(bwFT, stepInputData.getProjectID());
            bwFT.close();
        }
        catch(IOException exIO){
            logger.error("error writing feature file <" + featureFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing feature tsv file <" + featureFile + ">");
        }

        String freqFile = outFolder + FILESEPARATOR + stepInputData.getProjectID() + ".freq.tsv";
        try{
            BufferedWriter bwFQ = new BufferedWriter(new FileWriter(new File(freqFile)));
                featureSet.writeLengthDistribution(bwFQ, 1, longestFeature);
            bwFQ.close();
        }
        catch(IOException exIO){
            logger.error("error writing frequency file <" + freqFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing frequency file <" + freqFile + ">");
        }

        
        
        String fastaFile = outFolder + FILESEPARATOR + stepInputData.getProjectID() + ".fasta";
        try{
            BufferedWriter bwFA = new BufferedWriter(new FileWriter(new File(fastaFile)));
                featureSet.writeFeaturesAsFastA(bwFA);
            bwFA.close();
        }
        catch(IOException exIO){
            logger.error("error writing fasta feature file <" + fastaFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing fasta feature file <" + fastaFile + ">");
        }
    }

    /**
     * calculate average count across the mapped region for + Strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countCoverage5(int start, int stop) {
        int countTotal = 0;
        for (int i = start; i < stop; i++) {
            countTotal += coverage5[i];
        }
        return (double) countTotal / (double) (stop - start + 1);
    }

    /**
     * add counts to 3' coverage
     *
     * @param start
     * @param stop
     * @param counts
     */
    private void addCounts3(int start, int stop, int counts) {
        
        for (int i = stop; i < start; i++) {
            coverage3[i] += counts;
        }
        
    }

    /**
     * add counts to 5' coverage
     *
     * @param start
     * @param stop
     * @param counts
     */
    private void addCounts5(int start, int stop, int counts) {
        
        for (int i = start; i < stop; i++) {
            coverage5[i] += counts;
        }
        
    }

    /**
     * calculate average count across the mapped region for - Strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countCoverage3(int start, int stop) {
        int countTotal = 0;
        for (int i = start; i < stop; i++) {
            countTotal += coverage3[i];
        }
        return (double) countTotal / (double) (stop - start + 1);
    }

    /**
     * return string showing coverage over specified 3' region
     *
     * @param start
     * @param stop
     * @return
     */
    private String printCoverage3(int start, int stop, String delimiter) {
        String coverageStr = "";
        for (int i = start; i < stop; i++) {
            coverageStr = coverageStr.concat((Integer.toString(coverage3[i])) + delimiter);
        }
        return coverageStr;
    }

    /**
     * return string showing coverage over specified 3' region
     *
     * @param start
     * @param stop
     * @return
     */
    private String printCoverage5(int start, int stop, String delimiter) {
        String coverageStr = "";
        for (int i = start; i < stop; i++) {
            coverageStr = coverageStr.concat((Integer.toString(coverage5[i])) + delimiter);
        }
        return coverageStr;
    }

    /**
     * reset the 3' coverage array
     */
    private void clearCoverage3() {
        for (int i = 0; i < coverage3.length; i++) {
            coverage3[i] = 0;
        }        
    }

    /**
     * estimate dispersion of counts across the feature for + strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countDispersion5(int start, int stop) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = start; i < stop; i++) {
            stats.addValue(coverage5[i]);
        }
        return stats.getStandardDeviation();
    }

    /**
     * reset the 5' coverage array
     */
    private void clearCoverage5() {
        for (int i = 0; i < coverage5.length; i++) {
            coverage5[i] = 0;
        }        
    }

    /**
     * estimate dispersion of counts across the feature for - strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countDispersion3(int start, int stop) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = start; i < stop; i++) {
            stats.addValue(coverage3[i]);
        }
        return stats.getStandardDeviation();
    }

    /**
     * Verify Input Data for parsing SAM file for miRNAs
     *
     */    
    @Override
    public void verifyInputData()  throws IOException, NullPointerException{
        logger.info("verify input data");
        

        String pathToFasta = stepInputData.getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH;
        String genomeFastaFile = pathToFasta + FILESEPARATOR + "genome.fa";
        genomeFastaFile = genomeFastaFile.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();
        if (new File(genomeFastaFile + FILESEPARATOR + "genes.gtf").exists()==false){
            throw new IOException("no fasta file was found for reference genome <" 
                    + this.getReferenceGenome() + "> at location <" 
                    + genomeFastaFile + ">");

        }
        
        String annotationFile = "";
        String pathToAnnotation = stepInputData.getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_GENE_ANNOTATION;
        File f = new File(pathToAnnotation + FILESEPARATOR + "genes.gtf");
        if (new File(pathToAnnotation + FILESEPARATOR + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gtf";
        } else if (new File(pathToAnnotation + FILESEPARATOR + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gff";
        }
        if (annotationFile == null) {
            throw new IOException("no fasta file was found for reference genome <" 
                    + this.getReferenceGenome() + "> at location <"
                    + annotationFile + ">");
        }
        
        
                    
        // check the data files
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = (String)sampleData.getFastqFile1();
            
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
        }
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

        paramData.put(ID_REF_GENOME, "hsa");
        paramData.put(ID_BLEED, 2);
        paramData.put(ID_BASELINE, 5);
        paramData.put(ID_SHORTEST_FEATURE, 2);
        paramData.put(ID_LONGEST_FEATURE, 2);
        paramData.put(ID_MIN_COUNTS, 10);
        
        configData.put(STEP_ID_STRING, paramData);

        return configData;
    }





    @Override
    public void verifyOutputData() {
        
    }

    /**
     * @return the locationBleed
     */
    public int getLocationBleed() {
        return locationBleed;
    }

    /**
     * @param locationBleed the locationBleed to set
     */
    public void setLocationBleed(int locationBleed) {
        this.locationBleed = locationBleed;
    }

    /**
     * @return the baselinePercent
     */
    public int getBaselinePercent() {
        return baselinePercent;
    }

    /**
     * @param baselinePercent the baselinePercent to set
     */
    public void setBaselinePercent(int baselinePercent) {
        this.baselinePercent = baselinePercent;
    }

    /**
     * @return the analyzeIsomirs
     */
    public Boolean getAnalyzeIsomirs() {
        return analyzeIsomirs;
    }

    /**
     * @param analyzeIsomirs the analyzeIsomirs to set
     */
    public void setAnalyzeIsomirs(Boolean analyzeIsomirs) {
        this.analyzeIsomirs = analyzeIsomirs;
    }

    /**
     * @return the miRBaseRelease
     */
    public int getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(int miRBaseRelease) {
        this.miRBaseRelease = miRBaseRelease;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return ReferenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.ReferenceGenome = ReferenceGenome;
    }

    /**
     * @return the shortestRead
     */
    public int getShortestRead() {
        return shortestRead;
    }

    /**
     * @param shortestRead the shortestRead to set
     */
    public void setShortestRead(int shortestRead) {
        this.shortestRead = shortestRead;
    }

    /**
     * @return the longestRead
     */
    public int getLongestRead() {
        return longestRead;
    }

    /**
     * @param longestRead the longestRead to set
     */
    public void setLongestRead(int longestRead) {
        this.longestRead = longestRead;
    }

    /**
     * @return the min_counts
     */
    public int getMinCounts() {
        return minCounts;
    }

    /**
     * @param min_counts the min_counts to set
     */
    public void setMinCounts(int min_counts) {
        this.minCounts = min_counts;
    }

    /**
     * @return the separation
     */
    public int getSeparation() {
        return separation;
    }

    /**
     * @param separation the separation to set
     */
    public void setSeparation(int separation) {
        this.separation = separation;
    }
    
}
