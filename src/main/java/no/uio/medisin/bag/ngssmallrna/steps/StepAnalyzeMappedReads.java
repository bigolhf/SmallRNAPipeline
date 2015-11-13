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
 * this step performs additional analyses on a sample set of mapped reads in a
 * SAM file. It can look at features such as
 *
 * 1. distribution of start and end positions 2. distribution of mapped read
 * length
 *
 *
 * @author sr
 */
public class StepAnalyzeMappedReads extends NGSStep {
    
    private static final int COVERAGE_SPAN = 200;
    
    static Logger logger = LogManager.getLogger();
    static String FileSeparator = System.getProperty("file.separator");
    
    private static final String inFolder = "bowtie_genome_mapped";
    private static final String posAnalysisOutFolder = "position_analysis";
    
    private static final String infileExtension = ".trim.clp.gen.sam";
    private static final String positionsExtension = ".trim.clp.gen.pos.tsv";
    private static final String featuresExtension = ".trim.clp.gen.features.tsv";
    
    private StepInputData stepInputData;
    private StepResultData stepResultData;
    
    private GFFSet gffSet = new GFFSet();    
    private ArrayList<MappedRead> mappedReads = new ArrayList<>();
    private ArrayList<MappedRead> filteredReads = new ArrayList<>(); // generated from mappedReads
    
    int[] coverage5 = new int[COVERAGE_SPAN];    
    int[] coverage3 = new int[COVERAGE_SPAN];

    /*
     int[][] startPositions = new int[169100][9];
     int[][] readLengths = new int[169100][9];
    
     int[][] mrnaStartPositions = new int[4665][9];
     int[][] mrnaStopPositions  = new int[4665][9];
    
     int[] chrFeatureCount = new int[10];
     */
    public StepAnalyzeMappedReads(StepInputData sid) {
        stepInputData = sid;
    }

    /**
     * parse out SAM file to retrieve start and stop information for each read
     *
     * analyzeSAMStartPositionsParams.put("bleed",
     * this.getSamParseStartPosBleed());
     * analyzeSAMStartPositionsParams.put("separation",
     * this.getSamParseFeatureSeparation());
     * analyzeSAMStartPositionsParams.put("feature_types",
     * this.getSamParseFeatureTypes());
     * analyzeSAMStartPositionsParams.put("longest_feature",
     * this.getSamParseLongestFeature());
     * analyzeSAMStartPositionsParams.put("host",
     * this.getBowtieMappingReferenceGenome());
     * analyzeSAMStartPositionsParams.put("genomeReferenceGFFFile",
     * this.getGenomeAnnotationGFF());
     * analyzeSAMStartPositionsParams.put("bowtieMapGenomeRootFolder",
     * this.getGenomeRootFolder());
     * analyzeSAMStartPositionsParams.put("bowtieReferenceGenome",
     * this.getBowtieMappingReferenceGenome());
     *
     *
     * @param filename
     * @throws IOException
     */
    @Override
    public void execute() {
        try {
            stepInputData.verifyInputData();            
        } catch (IOException exIO) {
            logger.info("exception parsing InputData" + exIO);
        }

        /**
         * read genes.gtf or genes.gff3 file
         */
        String annotationFile = "";
        String pathToAnnotation = stepInputData.getStepParams().get("bowtieMapGenomeRootFolder")
                + FileSeparator + stepInputData.getStepParams().get("bowtieReferenceGenome") + "/Annotation/Genes";
        File f = new File(pathToAnnotation + FileSeparator + "genes.gtf");
        if (new File(pathToAnnotation + FileSeparator + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FileSeparator + "genes.gtf";
        } else if (new File(pathToAnnotation + FileSeparator + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FileSeparator + "genes.gff";
        }
        
        try {
            if (annotationFile == null) {
                throw new IOException("no annotation file was found for reference genome "
                        + stepInputData.getStepParams().get("bowtieReferenceGenome"));
            }
            gffSet.readGFF(annotationFile);
        } catch (IOException exIO) {
            logger.error("Exception trying to read Annotation file ");
            logger.error(exIO);
        }

        /**
         * read and parse the SAM files
         */
        int longestFeature = (int) stepInputData.getStepParams().get("longest_feature");
        
        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        String posAnalysisOutputFolder = pathToData + FileSeparator + posAnalysisOutFolder;
        posAnalysisOutputFolder = posAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(posAnalysisOutputFolder).mkdir();        
        if (fA) {
            logger.info("created output folder <" + posAnalysisOutputFolder + "> for results");
        }
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        String samLine = null;
        while (itSD.hasNext()) {
            SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
            
            String positionFile = posAnalysisOutFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", positionsExtension);
            String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
            
            logger.info("sam input file is " + samInputFile);
            logger.info("results will be written to " + positionFile);
            
            samLine = null;
            
            int lineCount = 0;
            try {
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                while ((samLine = brSAM.readLine()) != null) {
                    
                    logger.debug(samLine);
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
                logger.error(smIO);
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
                /*
                 Iterator itMR0 = mappedReads.iterator();
                 int m=0;
                 while(itMR0.hasNext()){
                 MappedRead mRead = (MappedRead) itMR0.next();
                 if(m % 1000 ==0)
                 logger.debug(mRead.toString());
                 m++;
                 if (m > 5000) break;
                 }
                 */                
                int bleed = (int) stepInputData.getStepParams().get("bleed");
                int separation = (int) stepInputData.getStepParams().get("separation");
                
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;
                
                MappedRead mappedRead = ((MappedRead) mappedReads.get(0));
                logger.debug(mappedRead.toString());
                String currentChr = mappedRead.getChr();
                String currentStrand = mappedRead.getStrand();

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
                
                if (currentStrand.equals(GFFEntry.PLUSSTRAND)) {
                    currentStart5 = mappedRead.getStartPos();
                    currentStop5 = mappedRead.getEndPos();
                    coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                } else {
                    currentStart3 = mappedRead.getStartPos();
                    currentStop3 = mappedRead.getEndPos();                    
                    coverage3Start = currentStart3 - COVERAGE_SPAN / 2;
                }
                
                String featureOutFile = pathToData + FileSeparator + posAnalysisOutFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", featuresExtension);
                featureOutFile = featureOutFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();
                int featureCount = 0;
                
                BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureOutFile)));
                Iterator itMR = mappedReads.iterator();
                Boolean startNewFeature3 = false;
                Boolean startNewFeature5 = false;
                while (itMR.hasNext()) {
                    mappedRead = (MappedRead) itMR.next();
                    logger.info(mappedRead.toString());
                    if (mappedRead.getChr().equals(currentChr) == false || mappedRead.getStrand().equals(currentStrand) == false) {
                        if (mappedRead.getChr().equals(currentChr) == false) {
                            logger.debug("different chr. start new feature");
                        }
                        if (mappedRead.getStrand().equals(currentStrand) == false) {
                            logger.debug("different strand. start new feature");                            
                        }
                        currentStrand = mappedRead.getStrand();
                        currentChr = mappedRead.getChr();
                        if (currentStrand.equals(GFFEntry.PLUSSTRAND)) {
                            startNewFeature3 = false;
                            startNewFeature5 = true;
                            currentStart5 = mappedRead.getStartPos();
                            currentStop5 = mappedRead.getEndPos();
                            coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                            addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                        } else {
                            startNewFeature3 = true;
                            startNewFeature5 = false;
                            currentStart3 = mappedRead.getStartPos();
                            currentStop3 = mappedRead.getEndPos();
                            coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                            addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                        }
                        
                    } else {
                        switch (currentStrand) {
                            case GFFEntry.PLUSSTRAND:
                                if (mappedRead.getStartPos() - currentStart5 > separation) {
                                    logger.debug("read too far away. start new 5' feature");                                    
                                    startNewFeature5 = true;
                                } else {
                                    if (mappedRead.getStartPos() < currentStart5) {
                                        logger.debug("extend 5´ start from " + currentStart5 + " to " + mappedRead.getStartPos());
                                        currentStart5 = mappedRead.getStartPos();
                                    }
                                    if (mappedRead.getEndPos() > currentStop5) {
                                        logger.debug("extend 5´ end from " + currentStop5 + " to " + mappedRead.getEndPos());
                                        currentStop5 = mappedRead.getEndPos();
                                    }
                                    try {
                                        logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                        this.addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                                        logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                    } catch (ArrayIndexOutOfBoundsException exAB) {
                                        logger.error("5' Array out of bounds");
                                    }
                                }
                                break;
                            
                            case GFFEntry.NEGSTRAND:
                                if ((mappedRead.getStartPos() - currentStart3 > separation) || (currentStop3 - mappedRead.getEndPos() > separation)) {
                                    logger.debug("read too far away. start new 3' feature");
                                    startNewFeature3 = true;
                                } else {
                                    logger.debug("add read to existing feature");
                                    if (mappedRead.getStartPos() > currentStart3) {
                                        logger.debug("extend 3´ start from " + currentStart3 + " to " + mappedRead.getStartPos());
                                        currentStart3 = mappedRead.getStartPos();
                                    }
                                    if (mappedRead.getEndPos() < currentStop3) {
                                        logger.debug("extend 3´ end from " + currentStop3 + " to " + mappedRead.getEndPos());
                                        currentStop3 = mappedRead.getEndPos();
                                    }
                                    try {
                                        logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                        this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                                        logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                    } catch (ArrayIndexOutOfBoundsException exAB) {
                                        logger.error("3' Array out of bounds");
                                        logger.error(exAB);
                                    }
                                }
                                break;
                            
                            case GFFEntry.UNKSTRAND:
                                logger.error("no Strand information for read, cannot process ");
                                logger.error(mappedRead.toString());
                                break;
                            
                        }
                        
                    }
                    
                    if (startNewFeature5) {
                        logger.debug("write out 5' feature");
                        if (currentStop5 - currentStart5 < longestFeature) {
                            bwFT.write("KEEP:\t");
                        } else {
                            bwFT.write("DROP:\t");
                        }
                        logger.debug(featureCount + "\t" + currentStart5 + "\t" + currentStop5 + "\t"
                                + (currentStop5 - currentStart5 + 1));
                        logger.debug(this.countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start));
                        logger.debug(countDispersion5(currentStart5 - coverage5Start, currentStop5 - coverage5Start));
                        logger.debug(printCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start, "|"));
                        bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + currentStart5 + "\t" + currentStop5 + "\t"
                                + (currentStop5 - currentStart5 + 1) + "\t" + this.countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start)
                                + "\t" + this.countDispersion5(currentStart5 - coverage5Start, currentStop5 - coverage5Start) + "\t"
                                + this.printCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start, "|") + "\n");
                        
                        logger.debug("set currentStart5 from " + currentStart5 + " to " + mappedRead.getStartPos());
                        logger.debug("set currentStop5 from " + currentStop5 + " to " + mappedRead.getEndPos());
                        currentStart5 = mappedRead.getStartPos();
                        currentStop5 = mappedRead.getEndPos();
                        logger.debug("set coverage range from " + (currentStart5 - COVERAGE_SPAN / 2) + " to " + (currentStart5 + COVERAGE_SPAN / 2));
                        coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                        clearCoverage5();
                        try {
                            logger.debug(printCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start, "|"));
                            this.addCounts5(mappedRead.getStartPos() - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                            logger.debug(printCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start, "|"));
                        } catch (ArrayIndexOutOfBoundsException exAB) {
                            logger.error("5' Array out of bounds when writing data");
                        }
                        
                        currentStrand = mappedRead.getStrand();
                        currentChr = mappedRead.getChr();
                        featureCount++;
                        startNewFeature5 = false;
                        
                    } else {
                        
                        if (startNewFeature3) {
                            logger.debug("write out 3' feature");
                            if (currentStop3 - currentStart3 < longestFeature) {
                                logger.debug("KEEP:\n");
                                bwFT.write("KEEP:\t");
                            } else {
                                logger.debug("DROP:\n");
                                bwFT.write("DROP:\t");
                            }
                            logger.debug(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + "\t" + currentStart3 + "\t" + currentStop3 + "\t"
                                    + (currentStop3 - currentStart3 + 1));
                            logger.debug(this.countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3));
                            logger.debug(countDispersion3(coverage3Start - currentStop3, coverage3Start - currentStart3));
                            logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                            bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + "\t" + currentStart3 + "\t" + currentStop3 + "\t"
                                    + (currentStop3 - currentStart3 + 1) + "\t" + this.countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3)
                                    + "\t" + this.countDispersion3(coverage3Start - currentStop3, coverage3Start - currentStart3) + "\t"
                                    + this.printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                            
                            logger.debug("set currentStart3 from " + currentStart3 + " to " + mappedRead.getStartPos());
                            logger.debug("set currentStop3 from " + currentStop3 + " to " + mappedRead.getEndPos());
                            currentStart3 = mappedRead.getStartPos();
                            currentStop3 = mappedRead.getEndPos();
                            logger.debug("set coverage range from " + (currentStart3 + COVERAGE_SPAN / 2) + " to " + (currentStart3 - COVERAGE_SPAN / 2));
                            coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                            clearCoverage3();
                            try {
                                logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                                logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                            } catch (ArrayIndexOutOfBoundsException exAB) {
                                logger.error("3' Array out of bounds when writing data");
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
            }
            
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
    public void verifyInputData() {
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()) {
            SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
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
    public void outputResultData() {
        
    }
    
}
