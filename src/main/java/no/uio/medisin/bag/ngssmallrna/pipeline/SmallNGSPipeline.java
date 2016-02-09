/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import static no.uio.medisin.bag.ngssmallrna.pipeline.SmallNGSCmd.logger;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeIsomiRDispersion;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapSingleReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCollapseReads;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepData;
import no.uio.medisin.bag.ngssmallrna.steps.StepParseSAMForMiRNAs;
import no.uio.medisin.bag.ngssmallrna.steps.StepDEwithEdgeR;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeMappedReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapPairedReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCleanUp;
import no.uio.medisin.bag.ngssmallrna.steps.StepInputData;
import no.uio.medisin.bag.ngssmallrna.steps.StepUnzipInputFiles;
import no.uio.medisin.bag.ngssmallrna.steps.StepSingleTrimAdapters;
import org.yaml.snakeyaml.Yaml;
import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.annotate.JsonAutoDetect;
import org.codehaus.jackson.annotate.JsonMethod;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

/**
 *
 * @author sr
 */
public class SmallNGSPipeline {
    
    static  Yaml                        yaml = new Yaml();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private Boolean                     generateSampleFiles = false;
    
    private String                      configurationFile = "";
    private String                      pipelineFile = "";
    private String                      dataFile = "";
    
    private HashMap                     pipelineConfigurationDataHash;    
    
    private String                      softwareRootFolder = "";
    private String                      zipSoftware = "";
    private String                      adapterTrimmingSoftware = "";
    private String                      fastq2fastaSoftware = "";
    private String                      collapseFastaSoftware = "";
    private String                      bowtieMappingCommand = "";
    private String                      rScriptCommand = "";
    
    private String                      genomeRootFolder = "";
    private String                      mirbaseRootFolder = "";
    
    
    private String                      trimAdapterFile = "";
    private int                         trimNoOfMismatches = 2;
    private int                         trimMinAlignScore = 7;
    private int                         trimNoOfThreads = 4;
    private int                         trimMinAvgReadQuality = 30;
    
    private String                      bowtieSMappingAlignMode = "";
    private int                         bowtieSMappingNoOfMismatches = 2;
    private int                         bowtieSMappingNoOfThreads = 4;
    private String                      bowtieSMappingReferenceGenome = "";
    
    private String                      bowtiePMappingAlignMode = "";
    private int                         bowtiePMappingNoOfMismatches = 2;
    private int                         bowtiePMappingNoOfThreads = 4;
    private String                      bowtiePMappingReferenceGenome = "";
    
    private int                         samParseForMiRNAsBleed = 2;
    private int                         samParseForMiRNAsBaselinePercent = 5;
    private String                      samParseForMiRNAsMiRBaseVersion = "21";
    private Boolean                     samParseForMiRNAsAnalyzeIsomirs = false;
    
    private int                         samParseStartPosBleed = 2;
    private int                         samParseFeatureSeparation = 10;
    private int                         samParseLongestFeature = 40;
    private int                         samParseShortestFeature = 14;
    private int                         samParseMinCounts = 10;
    private ArrayList<String>           samParseFeatureTypes = new ArrayList<>();
    
    private double                      analyzeIsomiRDispPVal = 0.05;
    
    private double                      diffExpressionPVal = 0.05;
    
    private ArrayList<String>           cleanupFiles = new ArrayList<>();
    private int                         cleanupNoOfThreads = 4;
    
    private PipelineData                pipelineData = new PipelineData();
    private ArrayList<SampleDataEntry>  SampleData = new ArrayList<>();
    private ArrayList<NGSStep>          NGSSteps = new ArrayList<>();

    
    
    /**
     * here we load the information from the three specified files
     * (run, pipeline and dataset) 
     * 
     * @throws IOException 
     */
    public void prepare_pipeline() throws IOException {      
        
        
        
        if(new File(this.getConfigurationFile()).exists()== false)
        {
            throw new IOException("configuration file <" + this.getConfigurationFile() + "> does not exist");
        }
        
        if(new File(this.getPipelineFile()).exists()== false)
        {
            throw new IOException("pipeline file <" + this.getPipelineFile() + "> does not exist");
        }
        
        if(new File(this.getDataFile()).exists()== false)
        {
            throw new IOException("data file <" + this.getDataFile() + "> does not exist");
        }
        
        this.readConfigurationFile();
        this.readPipelineFile();
        this.readDataFile();
        
        
    }
    
    
    /**
     * 
     * build each step and execute
     * @throws IOException, Exception
     */
    public void runPipeline() throws IOException, Exception{
        
        for (NGSRunStepData stepData: this.getPipelineData().getStepsData()){
            
            switch (stepData.getStepType()){
                
                case StepUnzipInputFiles.STEP_ID_STRING:
                    this.executeStepUnzipInputFiles(stepData);                    
                    break;
                    
                case StepSingleTrimAdapters.STEP_ID_STRING:
                    this.executeStepSingleTrimAdapters(stepData);
                    break;
                    
                case "CollapseReads":                
                    this.executeStepCollapseReads(stepData);
                    break;
                    
                case "BowtieMapSingleReads":
                    this.executeStepBowtieMapSingleReads(stepData);
                    break;
                    
                case "BowtieMapPairedReads":
                    this.executeStepBowtieMapPairedReads(stepData);
                    break;
                    
                case "ParseSAMForMiRNAs":
                    this.executeStepParseSAMForMiRNAs(stepData);
                    break;
                    
                case "AnalyzeStartPositions":
                    this.executeStepAnalyzeStartPositions(stepData);
                    break;
                    
                case "AnalyzeIsomiRDispersion":
                    this.executeStepAnalyzeIsomiRDispersion(stepData);
                    break;

                case "DifferentialExpression":
                    this.executeStepDifferentialExpression(stepData);
                    break;
                    
                case "Cleanup":
                    this.executeStepCleanup(stepData);
                    break;
                    
                case "exit":
                    return;
                    
            }
            
        }
        
    }

    
    
    /**
     * build and execute step to unzip input (Fastq) files
     * 
     * @param stepData 
     */
    private void executeStepUnzipInputFiles(NGSRunStepData stepData) throws IOException, Exception{
        
        StepInputData sidUnzip = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepUnzipInputFiles ngsUnzipFastqStep = new StepUnzipInputFiles(sidUnzip);        
        ngsUnzipFastqStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepUnzipInputFiles.STEP_ID_STRING));        
        ngsUnzipFastqStep.execute();
        
    }

    
    
    /**
     * build and execute step to trim adapter sequences from single end FASTQ files
     * @param stepData 
     */
    private void executeStepSingleTrimAdapters(NGSRunStepData stepData) throws IOException, Exception{

        HashMap trimAdapterParams = new HashMap();
        trimAdapterParams.put("trimAdapterSoftware",    this.getSoftwareRootFolder() + FileSeparator + this.getAdapterTrimmingSoftware());
        trimAdapterParams.put("trimAdapterFile",        this.getTrimAdapterFile());
        trimAdapterParams.put("trimNoOfMismatches",     this.getTrimNoOfMismatches());
        trimAdapterParams.put("trimMinAlignScore",      this.getTrimMinAlignScore());
        trimAdapterParams.put("trimNoOfThreads",        this.getTrimNoOfThreads());
        trimAdapterParams.put("trimMinAvgReadQuality",  this.getTrimMinAvgReadQuality());

        StepInputData sidTrim = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepSingleTrimAdapters ngsSingleTrimStep = new StepSingleTrimAdapters(sidTrim);
        ngsSingleTrimStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepSingleTrimAdapters.STEP_ID_STRING));        

        ngsSingleTrimStep.execute();
        
    }
    
    
    
    /**
     * build and execute step to collapse supplied FASTQ files into FASTA files
     * and then identify identical reads
     * 
     * @param stepData 
     */
    private void executeStepCollapseReads(NGSRunStepData stepData) throws IOException{

        HashMap collapseReadsParams = new HashMap();
        collapseReadsParams.put("fastqTofasta", this.getFastq2fastaSoftware());
        collapseReadsParams.put("collapseFasta", this.getCollapseFastaSoftware());

        StepInputData sidCollapse = new StepInputData(collapseReadsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepCollapseReads ngsCollapseStep = new StepCollapseReads(sidCollapse);

        ngsCollapseStep.execute();
        
    }
    
    
    /**
     * build and execute step to map single reads using Bowtie
     * 
     * @param stepData 
     */
    private void executeStepBowtieMapSingleReads(NGSRunStepData stepData) throws IOException{

        HashMap bowtieMapSingleReadsParams = new HashMap();
        bowtieMapSingleReadsParams.put("bowtieMappingCommand", this.getBowtieMappingCommand());
        bowtieMapSingleReadsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
        bowtieMapSingleReadsParams.put("bowtieReferenceGenome", this.getBowtieSMappingReferenceGenome());
        bowtieMapSingleReadsParams.put("bowtieMapAlignMode", this.getBowtieSMappingAlignMode());
        bowtieMapSingleReadsParams.put("bowtieMapNoOfMismatches", this.getBowtieSMappingNoOfMismatches()); //getBowtieMappingNoOfMismatches
        bowtieMapSingleReadsParams.put("bowtieNoOfThreads", this.getBowtieSMappingNoOfThreads());

        StepInputData sidMapSR = new StepInputData(bowtieMapSingleReadsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapSingleReads ngsBowtieMapReads = new StepBowtieMapSingleReads(sidMapSR);

        ngsBowtieMapReads.execute();
        
    }
    
    
    /**
     * build and execute step to map paired end reads using bowtie
     * 
     * @param stepData 
     */
    private void executeStepBowtieMapPairedReads(NGSRunStepData stepData) throws IOException{

        HashMap bowtieMapPairedReadsParams = new HashMap();
        bowtieMapPairedReadsParams.put("bowtieMappingCommand", this.getBowtieMappingCommand());
        bowtieMapPairedReadsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
        bowtieMapPairedReadsParams.put("bowtieReferenceGenome", this.getBowtiePMappingReferenceGenome());
        bowtieMapPairedReadsParams.put("bowtieMapAlignMode", this.getBowtiePMappingAlignMode());
        bowtieMapPairedReadsParams.put("bowtieMapNoOfMismatches", this.getBowtiePMappingNoOfMismatches()); //getBowtieMappingNoOfMismatches
        bowtieMapPairedReadsParams.put("bowtieNoOfThreads", this.getBowtiePMappingNoOfThreads());

        StepInputData sidMapPR = new StepInputData(bowtieMapPairedReadsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapPairedReads ngsBowtieMapPairedReads = new StepBowtieMapPairedReads(sidMapPR);

        ngsBowtieMapPairedReads.execute();
        
    }
    
    
    /**
     * build and execute step to parse a SAM alignment file to identify miRNAs
     * 
     * @param stepData 
     */
    private void executeStepParseSAMForMiRNAs(NGSRunStepData stepData) throws IOException{
        
       HashMap parseSAMmiRNAsParams = new HashMap();
       parseSAMmiRNAsParams.put("inputFolder", stepData.getInputFileList());
       parseSAMmiRNAsParams.put("outputFolder", stepData.getOutputFileList());
       parseSAMmiRNAsParams.put("bleed", this.getSamParseForMiRNAsBleed());
       parseSAMmiRNAsParams.put("baseline_percent", this.getSamParseForMiRNAsBaselinePercent());
       parseSAMmiRNAsParams.put("host", this.getBowtieSMappingReferenceGenome());
       parseSAMmiRNAsParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());
       parseSAMmiRNAsParams.put("analyze_isomirs", this.getSamParseForMiRNAsAnalyzeIsomirs());

       StepInputData sidSAM = new StepInputData(parseSAMmiRNAsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
       StepParseSAMForMiRNAs ngsParseSAMForMiRNAs = new StepParseSAMForMiRNAs(sidSAM);

       ngsParseSAMForMiRNAs.execute();
       
    }
    
    
    /**
     * build and execute step to analyze start reads in SAM file to identify
     * reads that may be associated with novel smallRNAs
     * 
     * @param stepData 
     */
    private void executeStepAnalyzeStartPositions(NGSRunStepData stepData) throws IOException{
        
        HashMap analyzeSAMStartPositionsParams = new HashMap();
        analyzeSAMStartPositionsParams.put("bleed", this.getSamParseStartPosBleed());
        analyzeSAMStartPositionsParams.put("separation", this.getSamParseFeatureSeparation());
        analyzeSAMStartPositionsParams.put("shortest_feature", this.getSamParseShortestFeature());    
        analyzeSAMStartPositionsParams.put("longest_feature", this.getSamParseLongestFeature());    
        analyzeSAMStartPositionsParams.put("min_counts", this.getSamParseMinCounts());
        analyzeSAMStartPositionsParams.put("feature_types", this.getSamParseFeatureTypes());
        analyzeSAMStartPositionsParams.put("host", this.getBowtieSMappingReferenceGenome());
        analyzeSAMStartPositionsParams.put("genomeReferenceGFFFile", this.getGenomeAnnotationGFF());
        analyzeSAMStartPositionsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
        analyzeSAMStartPositionsParams.put("bowtieReferenceGenome", this.getBowtieSMappingReferenceGenome());

        StepInputData sidStart = new StepInputData(analyzeSAMStartPositionsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeMappedReads ngsAnalyzeSAMStartPos = new StepAnalyzeMappedReads(sidStart);

        ngsAnalyzeSAMStartPos.execute();
        
    }
    
    
    /**
     * build and execute step to analyze SAM file to identify isomiR populations
     * 
     * @param stepData 
     */
    private void executeStepAnalyzeIsomiRDispersion(NGSRunStepData stepData) throws IOException{
        
        HashMap isomiRDispersionAnalysisParams = new HashMap();
        isomiRDispersionAnalysisParams.put("pvalue", this.getAnalyzeIsomiRDispPVal());
        isomiRDispersionAnalysisParams.put("host", this.getBowtieSMappingReferenceGenome());
        isomiRDispersionAnalysisParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());

        StepInputData sidIsoDisp = new StepInputData(isomiRDispersionAnalysisParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeIsomiRDispersion analyzeIsomiRDispersions = new StepAnalyzeIsomiRDispersion(sidIsoDisp);

        analyzeIsomiRDispersions.execute();
        
    }
    
    
    /**
     * build and execute step to perform differential expression based on 
     * supplied count files and grouping information
     * 
     * @param stepData 
     */
    private void executeStepDifferentialExpression(NGSRunStepData stepData) throws IOException{
        /*
        Need to add:
            Number of Tags to report
            minimum counts for cut off
            plot size: width and height, resolution

        */
        HashMap diffExpressionAnalysisParams = new HashMap();
        diffExpressionAnalysisParams.put("rScriptCommand", this.getrScriptCommand());
        diffExpressionAnalysisParams.put("pvalue", this.getDiffExpressionPVal());
        diffExpressionAnalysisParams.put("host", this.getBowtieSMappingReferenceGenome());
        diffExpressionAnalysisParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());

        StepInputData siodDiffExpr = new StepInputData(diffExpressionAnalysisParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepDEwithEdgeR edgeRDE = new StepDEwithEdgeR(siodDiffExpr);

        edgeRDE.execute();
        
    }
    
    
    /**
     * build and execute step to perform clean up after an analysis is complete
     * 
     * @param stepData 
     */
    private void executeStepCleanup(NGSRunStepData stepData) throws IOException{

        HashMap cleanUpParams = new HashMap();
        cleanUpParams.put("zipSoftware",          this.getZipSoftware());
        cleanUpParams.put("trimNoOfThreads",        this.getTrimNoOfThreads());
        cleanUpParams.put("fileTypes",              this.getCleanupFiles());

        StepInputData sidCleanUp = new StepInputData(cleanUpParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepCleanUp cleanUp = new StepCleanUp(sidCleanUp);

        cleanUp.execute();
        
    }
    
    
    
    
    /**
     * read file containing information about the data to be analyzed
     * 
     * @throws IOException 
     * 
     */
    public void readDataFile() throws IOException{
        
        
        logger.info("\nread data file");
        String line = "";
        BufferedReader bwData = new BufferedReader(new FileReader(new File(this.getDataFile())));
            bwData.readLine(); // skip header line
            while((line=bwData.readLine())!= null){
                if (line.startsWith("#")) 
                    continue;
                
                String tokens[] = line.split("\t");
                String file1 = tokens[0].split(",")[0].trim();
                String file2 = tokens[0].split(",")[1].trim();
                String source = tokens[1];
                String condition = tokens[2];
                String time = tokens[3];
                String note = tokens[4];
                
                getSampleData().add(new SampleDataEntry(file1, file2, source, condition, time, note));
                SampleDataEntry x = new SampleDataEntry(file1, file2, source, condition, time, note);
                
                logger.info(line);
                
            }
        bwData.close();
        logger.info("\nread " + getSampleData().size() + " entries");
        
        
    }
    
    
    
    
    /**
     * read configuration settings for the pipeline
     * there should be settings for each step. If the data is missing for a step
     * the method should return a warning.
     * 
     * @throws IOException 
     */
    public void readConfigurationFile() throws IOException{
        
        logger.info("read pipeline configuration file");
        pipelineConfigurationDataHash = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));
        
        HashMap dataOptions = (HashMap) pipelineConfigurationDataHash.get("data");
        this.setGenomeRootFolder((String) dataOptions.get("genome_root_folder"));
        this.setMirbaseRootFolder((String) dataOptions.get("mirbase_folder"));

        HashMap softwareOptions = (HashMap) pipelineConfigurationDataHash.get("software");
        this.setSoftwareRootFolder((String) softwareOptions.get("root_folder"));
        this.setZipSoftware((String) softwareOptions.get("unzip"));
        this.setAdapterTrimmingSoftware((String) softwareOptions.get("adapter_trimming"));
        this.setFastq2fastaSoftware((String) softwareOptions.get("fastq_to_fasta"));
        this.setCollapseFastaSoftware((String) softwareOptions.get("fastx_collapser"));
        this.setBowtieMappingCommand((String) softwareOptions.get("bowtie_mapping_command"));
        this.setrScriptCommand((String) softwareOptions.get("rscript_command"));
        
        HashMap unzipOptions = (HashMap) pipelineConfigurationDataHash.get("unzip");
        this.setUnzipNoOfThreads((int) unzipOptions.get("no_of_threads"));
        
        HashMap trimAdapterOptions = (HashMap) pipelineConfigurationDataHash.get("single_adapter_trimming");
        this.setTrimAdapterFile((String) trimAdapterOptions.get("adapter_file"));
        this.setTrimMinAlignScore((int) trimAdapterOptions.get("no_of_mismatches"));
        this.setTrimNoOfMismatches((int) trimAdapterOptions.get("min_align_score"));
        this.setTrimNoOfThreads((int) trimAdapterOptions.get("no_of_threads"));
        this.setTrimMinAvgReadQuality((int) trimAdapterOptions.get("min_avg_read_qual"));
        
        HashMap mapSingleReadsOptions = (HashMap) pipelineConfigurationDataHash.get("bowtie_single_mapping");
        this.setBowtieSMappingAlignMode((String) mapSingleReadsOptions.get("alignment_mode"));
        this.setBowtieSMappingNoOfMismatches((int) mapSingleReadsOptions.get("no_of_mismatches"));
        this.setBowtieSMappingNoOfThreads((int) mapSingleReadsOptions.get("no_of_threads"));
        this.setBowtieSMappingReferenceGenome((String) mapSingleReadsOptions.get("host"));
        
        HashMap mapPairedReadsOptions = (HashMap) pipelineConfigurationDataHash.get("bowtie_paired_mapping");
        this.setBowtiePMappingAlignMode((String) mapPairedReadsOptions.get("alignment_mode"));
        this.setBowtiePMappingNoOfMismatches((int) mapPairedReadsOptions.get("no_of_mismatches"));
        this.setBowtiePMappingNoOfThreads((int) mapPairedReadsOptions.get("no_of_threads"));
        this.setBowtiePMappingReferenceGenome((String) mapPairedReadsOptions.get("host"));
        
        HashMap processSAMForMiRNAsOptions  = (HashMap) pipelineConfigurationDataHash.get("sam_mirna_processing");
        this.setSamParseForMiRNAsBleed((int) processSAMForMiRNAsOptions.get("bleed"));
        this.setSamParseForMiRNAsMiRBaseVersion( String.valueOf(processSAMForMiRNAsOptions.get("mirbase_release")));
        this.setSamParseForMiRNAsBaselinePercent((int) processSAMForMiRNAsOptions.get("baseline_percent"));
        this.setSamParseForMiRNAsAnalyzeIsomirs((Boolean) processSAMForMiRNAsOptions.get("analyze_isomirs"));
        
        HashMap processSAMStartPosOptions = (HashMap) pipelineConfigurationDataHash.get("sam_startpos_processing");
        this.setSamParseStartPosBleed((int) processSAMStartPosOptions.get("bleed"));
        this.setSamParseFeatureSeparation((int) processSAMStartPosOptions.get("separation"));
        this.setSamParseShortestFeature((int) processSAMStartPosOptions.get("shortest_feature"));
        this.setSamParseLongestFeature((int) processSAMStartPosOptions.get("longest_feature"));
        this.setSamParseMinCounts((int) processSAMStartPosOptions.get("min_counts"));
        this.setSamParseFeatureTypes((ArrayList<String>)processSAMStartPosOptions.get("feature_types"));
        
        HashMap analyzeIsomiRDispersionOptions = (HashMap) pipelineConfigurationDataHash.get("analyze_isomir_dispersion");
        this.setAnalyzeIsomiRDispPVal((double) analyzeIsomiRDispersionOptions.get("pvalue"));
        
        HashMap diffExpressionOptions = (HashMap) pipelineConfigurationDataHash.get("differential_expression");
        this.setDiffExpressionPVal((double) diffExpressionOptions.get("pvalue"));
        
        HashMap cleanupOptions = (HashMap) pipelineConfigurationDataHash.get("cleanup");
        this.setCleanupFiles((ArrayList<String>)cleanupOptions.get("file_types"));
        this.setCleanupNoOfThreads((int) cleanupOptions.get("no_of_threads"));
        
        
        logger.info("done\n");
        
    }
    
    
    
    public void writeSampleConfigurationFile() throws IOException{
        logger.info("write example pipeline configuration file");
        HashMap pipelineConfiguration = new HashMap();
        
    
        FileWriter sampleFileWriter = new FileWriter(new File("pipelineConfiguration.sample.yaml"));
        
        
        
        
        yaml.dump(pipelineConfiguration, sampleFileWriter);
        
    }
    /**
     * 
     */
    public void readPipelineFile(){
        
	ObjectMapper mapper = new ObjectMapper();
        mapper.setVisibility(JsonMethod.FIELD, JsonAutoDetect.Visibility.ANY);
        
	try { 
            
            
                logger.info("read pipeline file <" + this.getPipelineFile() + ">");
                pipelineData = mapper.readValue(new File(this.getPipelineFile()), PipelineData.class); 
                logger.info("done");
                logger.info(getPipelineData().toString());
                
 
	} catch (JsonGenerationException e) {
 
		e.printStackTrace();
 
	} catch (JsonMappingException e) {
 
		e.printStackTrace();
 
	} catch (IOException e) {
 
		e.printStackTrace();
 
	}
        
        
    }
    
    
    
    
    /**
     * return path to the file containing location information for the known miRNAs
     * for the specified reference genome
     * 
     * @return 
     */
    public String getMiRBaseHostGFF(){
        return this.mirbaseRootFolder + FileSeparator + this.samParseForMiRNAsMiRBaseVersion + FileSeparator + this.getBowtieSMappingReferenceGenome() + ".gff3";
    }
    
    
    
    
    /**
     * return path to the GFF/GTF file containing annotation information 
     * for the specified reference genome
     * 
     * @return 
     */
    public String getGenomeAnnotationGFF(){
        return this.getGenomeRootFolder() + FileSeparator + this.getBowtieSMappingReferenceGenome() + FileSeparator + "Annotation/Genes/genes.gtf";
    }
    
    
    
    
    /**
     * return path to root folder for MirBase version
     * 
     * @return 
     */
    public String getMirBaseRootFolder(){
        return this.mirbaseRootFolder + FileSeparator + this.samParseForMiRNAsMiRBaseVersion;
    }
    
    
    /**
     * @return the configurationFile
     */
    public String getConfigurationFile() {
        return configurationFile;
    }

    /**
     * @param configurationFile the configurationFile to set
     */
    public void setConfigurationFile(String configurationFile) {
        this.configurationFile = configurationFile;
    }

    /**
     * @return the pipelineFile
     */
    public String getPipelineFile() {
        return pipelineFile;
    }

    /**
     * @param pipelineFile the pipelineFile to set
     */
    public void setPipelineFile(String pipelineFile) {
        this.pipelineFile = pipelineFile;
    }

    /**
     * @return the dataFile
     */
    public String getDataFile() {
        return dataFile;
    }

    /**
     * @param dataFile the dataFile to set
     */
    public void setDataFile(String dataFile) {
        this.dataFile = dataFile;
    }

    /**
     * @return the pipelineData
     */
    public PipelineData getPipelineData() {
        return pipelineData;
    }

    /**
     * @return the softwareRootFolder
     */
    public String getSoftwareRootFolder() {
        return softwareRootFolder;
    }

    /**
     * @param softwareRootFolder the softwareRootFolder to set
     */
    public void setSoftwareRootFolder(String softwareRootFolder) {
        this.softwareRootFolder = softwareRootFolder;
    }

    /**
     * @return the adapterTrimmingSoftware
     */
    public String getAdapterTrimmingSoftware() {
        return adapterTrimmingSoftware;
    }

    /**
     * @param adapterTrimmingSoftware the adapterTrimmingSoftware to set
     */
    public void setAdapterTrimmingSoftware(String adapterTrimmingSoftware) {
        this.adapterTrimmingSoftware = adapterTrimmingSoftware;
    }

    /**
     * @return the NGSSteps
     */
    public ArrayList<NGSStep> getNGSSteps() {
        return NGSSteps;
    }

    /**
     * @return the trimAdapterFile
     */
    public String getTrimAdapterFile() {
        return trimAdapterFile;
    }

    /**
     * @param trimAdapterFile the trimAdapterFile to set
     */
    public void setTrimAdapterFile(String trimAdapterFile) {
        this.trimAdapterFile = trimAdapterFile;
    }

    /**
     * @return the trimNoOfMismatches
     */
    public int getTrimNoOfMismatches() {
        return trimNoOfMismatches;
    }

    /**
     * @param trimNoOfMismatches the trimNoOfMismatches to set
     */
    public void setTrimNoOfMismatches(int trimNoOfMismatches) {
        this.trimNoOfMismatches = trimNoOfMismatches;
    }

    /**
     * @return the trimMinAlignScore
     */
    public int getTrimMinAlignScore() {
        return trimMinAlignScore;
    }

    /**
     * @param trimMinAlignScore the trimMinAlignScore to set
     */
    public void setTrimMinAlignScore(int trimMinAlignScore) {
        this.trimMinAlignScore = trimMinAlignScore;
    }

    /**
     * @return the trimNoOfThreads
     */
    public int getTrimNoOfThreads() {
        return trimNoOfThreads;
    }

    /**
     * @param trimNoOfThreads the trimNoOfThreads to set
     */
    public void setTrimNoOfThreads(int trimNoOfThreads) {
        this.trimNoOfThreads = trimNoOfThreads;
    }

    /**
     * @return the SampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return SampleData;
    }

    /**
     * @return the trimMinReadQuality
     */
    public int getTrimMinAvgReadQuality() {
        return trimMinAvgReadQuality;
    }

    /**
     * @param trimMinReadQuality the trimMinReadQuality to set
     */
    public void setTrimMinAvgReadQuality(int trimMinReadQuality) {
        this.trimMinAvgReadQuality = trimMinReadQuality;
    }

    /**
     * @return the fastq2fastaSoftware
     */
    public String getFastq2fastaSoftware() {
        return fastq2fastaSoftware;
    }

    /**
     * @param fastq2fastaSoftware the fastq2fastaSoftware to set
     */
    public void setFastq2fastaSoftware(String fastq2fastaSoftware) {
        this.fastq2fastaSoftware = fastq2fastaSoftware;
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

    /**
     * @return the bowtieMappingAlignMode
     */
    public String getBowtieSMappingAlignMode() {
        return bowtieSMappingAlignMode;
    }

    /**
     * @param bowtieMappingAlignMode the bowtieMappingAlignMode to set
     */
    public void setBowtieSMappingAlignMode(String bowtieMappingAlignMode) {
        this.bowtieSMappingAlignMode = bowtieMappingAlignMode;
    }

    /**
     * @return the bowtieMappingNoOfMismatches
     */
    public int getBowtieSMappingNoOfMismatches() {
        return bowtieSMappingNoOfMismatches;
    }

    /**
     * @param bowtieMappingNoOfMismatches the bowtieMappingNoOfMismatches to set
     */
    public void setBowtieSMappingNoOfMismatches(int bowtieMappingNoOfMismatches) {
        this.bowtieSMappingNoOfMismatches = bowtieMappingNoOfMismatches;
    }

    /**
     * @return the bowtieMappingNoOfThreads
     */
    public int getBowtieSMappingNoOfThreads() {
        return bowtieSMappingNoOfThreads;
    }

    /**
     * @param bowtieMappingNoOfThreads the bowtieMappingNoOfThreads to set
     */
    public void setBowtieSMappingNoOfThreads(int bowtieMappingNoOfThreads) {
        this.bowtieSMappingNoOfThreads = bowtieMappingNoOfThreads;
    }

    /**
     * @return the bowtieMappingReferenceGenome
     */
    public String getBowtieSMappingReferenceGenome() {
        return bowtieSMappingReferenceGenome;
    }

    /**
     * @param bowtieMappingReferenceGenome the bowtieMappingReferenceGenome to set
     */
    public void setBowtieSMappingReferenceGenome(String bowtieMappingReferenceGenome) {
        this.bowtieSMappingReferenceGenome = bowtieMappingReferenceGenome;
    }

    /**
     * @return the genomeRootFolder
     */
    public String getGenomeRootFolder() {
        return genomeRootFolder;
    }

    /**
     * @param genomeRootFolder the genomeRootFolder to set
     */
    public void setGenomeRootFolder(String genomeRootFolder) {
        this.genomeRootFolder = genomeRootFolder;
    }

    /**
     * @return the samParseForMiRNAsBleed
     */
    public int getSamParseForMiRNAsBleed() {
        return samParseForMiRNAsBleed;
    }

    /**
     * @param samParseForMiRNAsBleed the samParseForMiRNAsBleed to set
     */
    public void setSamParseForMiRNAsBleed(int samParseForMiRNAsBleed) {
        this.samParseForMiRNAsBleed = samParseForMiRNAsBleed;
    }

    /**
     * @return the samParseForMiRNAsMiRBaseVersion
     */
    public String getSamParseForMiRNAsMiRBaseVersion() {
        return samParseForMiRNAsMiRBaseVersion;
    }

    /**
     * @param samParseForMiRNAsMiRBaseVersion the samParseForMiRNAsMiRBaseVersion to set
     */
    public void setSamParseForMiRNAsMiRBaseVersion(String samParseForMiRNAsMiRBaseVersion) {
        this.samParseForMiRNAsMiRBaseVersion = samParseForMiRNAsMiRBaseVersion;
    }

    /**
     * @return the mirbaseRootFolder
     */
    public String getMirBaseVersionRoot() {
        return mirbaseRootFolder;
    }

    /**
     * @param mirbaseRootFolder the mirbaseRootFolder to set
     */
    public void setMirbaseRootFolder(String mirbaseRootFolder) {
        this.mirbaseRootFolder = mirbaseRootFolder;
    }

    /**
     * @return the samParseForMiRNAsBaselinePercent
     */
    public int getSamParseForMiRNAsBaselinePercent() {
        return samParseForMiRNAsBaselinePercent;
    }

    /**
     * @param samParseForMiRNAsBaselinePercent the samParseForMiRNAsBaselinePercent to set
     */
    public void setSamParseForMiRNAsBaselinePercent(int samParseForMiRNAsBaselinePercent) {
        this.samParseForMiRNAsBaselinePercent = samParseForMiRNAsBaselinePercent;
    }

    /**
     * @return the samParseStartPosBleed
     */
    public int getSamParseStartPosBleed() {
        return samParseStartPosBleed;
    }

    /**
     * @param samParseStartPosBleed the samParseStartPosBleed to set
     */
    public void setSamParseStartPosBleed(int samParseStartPosBleed) {
        this.samParseStartPosBleed = samParseStartPosBleed;
    }

    /**
     * @return the samParseFeatureTypes
     */
    public ArrayList<String> getSamParseFeatureTypes() {
        return samParseFeatureTypes;
    }

    /**
     * @param samParseFeatureTypes the samParseFeatureTypes to set
     */
    public void setSamParseFeatureTypes(ArrayList<String> samParseFeatureTypes) {
        this.samParseFeatureTypes = samParseFeatureTypes;
    }

    /**
     * @return the analyzeIsomiRDispPVal
     */
    public double getAnalyzeIsomiRDispPVal() {
        return analyzeIsomiRDispPVal;
    }

    /**
     * @param analyzeIsomiRDispPVal the analyzeIsomiRDispPVal to set
     */
    public void setAnalyzeIsomiRDispPVal(double analyzeIsomiRDispPVal) {
        this.analyzeIsomiRDispPVal = analyzeIsomiRDispPVal;
    }

    /**
     * @return the diffExpressionPVal
     */
    public double getDiffExpressionPVal() {
        return diffExpressionPVal;
    }

    /**
     * @param diffExpressionPVal the diffExpressionPVal to set
     */
    public void setDiffExpressionPVal(double diffExpressionPVal) {
        this.diffExpressionPVal = diffExpressionPVal;
    }

    /**
     * @return the mappingCommand
     */
    public String getBowtieMappingCommand() {
        return bowtieMappingCommand;
    }

    /**
     * @param mappingCommand the mappingCommand to set
     */
    public void setBowtieMappingCommand(String mappingCommand) {
        this.bowtieMappingCommand = mappingCommand;
    }

    /**
     * @return the samParseFeatureSeparation
     */
    public int getSamParseFeatureSeparation() {
        return samParseFeatureSeparation;
    }

    /**
     * @param samParseFeatureSeparation the samParseFeatureSeparation to set
     */
    public void setSamParseFeatureSeparation(int samParseFeatureSeparation) {
        this.samParseFeatureSeparation = samParseFeatureSeparation;
    }

    /**
     * @return the samParseLongestFeature
     */
    public int getSamParseLongestFeature() {
        return samParseLongestFeature;
    }

    /**
     * @param samParseLongestFeature the samParseLongestFeature to set
     */
    public void setSamParseLongestFeature(int samParseLongestFeature) {
        this.samParseLongestFeature = samParseLongestFeature;
    }

    /**
     * @return the samParseMinCounts
     */
    public int getSamParseMinCounts() {
        return samParseMinCounts;
    }

    /**
     * @param samParseMinCounts the samParseMinCounts to set
     */
    public void setSamParseMinCounts(int samParseMinCounts) {
        this.samParseMinCounts = samParseMinCounts;
    }

    /**
     * @return the samParseShortestFeature
     */
    public int getSamParseShortestFeature() {
        return samParseShortestFeature;
    }

    /**
     * @param samParseShortestFeature the samParseShortestFeature to set
     */
    public void setSamParseShortestFeature(int samParseShortestFeature) {
        this.samParseShortestFeature = samParseShortestFeature;
    }

    /**
     * @return the samParseForMiRNAsAnalyzeIsomirs
     */
    public Boolean getSamParseForMiRNAsAnalyzeIsomirs() {
        return samParseForMiRNAsAnalyzeIsomirs;
    }

    /**
     * @param samParseForMiRNAsAnalyzeIsomirs the samParseForMiRNAsAnalyzeIsomirs to set
     */
    public void setSamParseForMiRNAsAnalyzeIsomirs(Boolean samParseForMiRNAsAnalyzeIsomirs) {
        this.samParseForMiRNAsAnalyzeIsomirs = samParseForMiRNAsAnalyzeIsomirs;
    }

    /**
     * @return the unzipNoOfThreads
     */
    public int getUnzipNoOfThreads() {
        return unzipNoOfThreads;
    }

    /**
     * @param unzipNoOfThreads the unzipNoOfThreads to set
     */
    public void setUnzipNoOfThreads(int unzipNoOfThreads) {
        this.unzipNoOfThreads = unzipNoOfThreads;
    }

    /**
     * @return the cleanupFiles
     */
    public ArrayList<String> getCleanupFiles() {
        return cleanupFiles;
    }

    /**
     * @param cleanupFiles the cleanupFiles to set
     */
    public void setCleanupFiles(ArrayList<String> cleanupFiles) {
        this.cleanupFiles = cleanupFiles;
    }

    /**
     * @return the zipSoftware
     */
    public String getZipSoftware() {
        return zipSoftware;
    }

    /**
     * @param zipSoftware the zipSoftware to set
     */
    public void setZipSoftware(String zipSoftware) {
        this.zipSoftware = zipSoftware;
    }

    /**
     * @return the cleanupNoOfThreads
     */
    public int getCleanupNoOfThreads() {
        return cleanupNoOfThreads;
    }

    /**
     * @param cleanupNoOfThreads the cleanupNoOfThreads to set
     */
    public void setCleanupNoOfThreads(int cleanupNoOfThreads) {
        this.cleanupNoOfThreads = cleanupNoOfThreads;
    }

    /**
     * @return the bowtiePMappingAlignMode
     */
    public String getBowtiePMappingAlignMode() {
        return bowtiePMappingAlignMode;
    }

    /**
     * @param bowtiePMappingAlignMode the bowtiePMappingAlignMode to set
     */
    public void setBowtiePMappingAlignMode(String bowtiePMappingAlignMode) {
        this.bowtiePMappingAlignMode = bowtiePMappingAlignMode;
    }

    /**
     * @return the bowtiePMappingNoOfMismatches
     */
    public int getBowtiePMappingNoOfMismatches() {
        return bowtiePMappingNoOfMismatches;
    }

    /**
     * @param bowtiePMappingNoOfMismatches the bowtiePMappingNoOfMismatches to set
     */
    public void setBowtiePMappingNoOfMismatches(int bowtiePMappingNoOfMismatches) {
        this.bowtiePMappingNoOfMismatches = bowtiePMappingNoOfMismatches;
    }

    /**
     * @return the bowtiePMappingNoOfThreads
     */
    public int getBowtiePMappingNoOfThreads() {
        return bowtiePMappingNoOfThreads;
    }

    /**
     * @param bowtiePMappingNoOfThreads the bowtiePMappingNoOfThreads to set
     */
    public void setBowtiePMappingNoOfThreads(int bowtiePMappingNoOfThreads) {
        this.bowtiePMappingNoOfThreads = bowtiePMappingNoOfThreads;
    }

    /**
     * @return the bowtiePMappingReferenceGenome
     */
    public String getBowtiePMappingReferenceGenome() {
        return bowtiePMappingReferenceGenome;
    }

    /**
     * @param bowtiePMappingReferenceGenome the bowtiePMappingReferenceGenome to set
     */
    public void setBowtiePMappingReferenceGenome(String bowtiePMappingReferenceGenome) {
        this.bowtiePMappingReferenceGenome = bowtiePMappingReferenceGenome;
    }

    /**
     * @return the rScriptCommand
     */
    public String getrScriptCommand() {
        return rScriptCommand;
    }

    /**
     * @param rScriptCommand the rScriptCommand to set
     */
    public void setrScriptCommand(String rScriptCommand) {
        this.rScriptCommand = rScriptCommand;
    }

    /**
     * @return the generateSampleFiles
     */
    public Boolean getGenerateSampleFiles() {
        return generateSampleFiles;
    }

    /**
     * @param generateSampleFiles the generateSampleFiles to set
     */
    public void setGenerateSampleFiles(Boolean generateSampleFiles) {
        this.generateSampleFiles = generateSampleFiles;
    }
    
    
    
    
}
