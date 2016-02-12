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
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import static no.uio.medisin.bag.ngssmallrna.pipeline.SmallNGSCmd.logger;
import no.uio.medisin.bag.ngssmallrna.steps.NGSBase;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeIsomiRDispersion;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapSingleReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCollapseReads;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepData;
import no.uio.medisin.bag.ngssmallrna.steps.StepParseSAMForMiRNAs;
import no.uio.medisin.bag.ngssmallrna.steps.StepDEwithEdgeR;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeSAMforStartPositions;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapPairedReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCleanUp;
import no.uio.medisin.bag.ngssmallrna.steps.StepInputData;
import no.uio.medisin.bag.ngssmallrna.steps.StepUnzipInputFiles;
import no.uio.medisin.bag.ngssmallrna.steps.StepSingleTrimAdapters;
import org.yaml.snakeyaml.Yaml;
import org.codehaus.jackson.annotate.JsonAutoDetect;
import org.codehaus.jackson.annotate.JsonMethod;
import org.codehaus.jackson.map.ObjectMapper;
import org.yaml.snakeyaml.DumperOptions;

/**
 *
 * @author sr
 */
public class SmallNGSPipeline {
    
    static  final String                FILE_SEPARATOR = System.getProperty("file.separator");
    
    private String                      configurationFile = "";
    private String                      pipelineFile = "";
    private String                      dataFile = "";
    static  Yaml                        yaml;
    
    private Boolean                     generateExampleConfigurationFile;
    
    private HashMap                     pipelineConfigurationDataHash;    
    private ReferenceDataLocations      refDataLocations;
    
    
    
    private ArrayList<String>           cleanupFiles    = new ArrayList<>();
    private PipelineData                pipelineData    = new PipelineData();
    private ArrayList<SampleDataEntry>  SampleData      = new ArrayList<>();

    private ArrayList<NGSBase>          ngsSteps        = new ArrayList<>();
    
    /**
     * here we load the information from the three specified files
     * (run, pipeline and dataset) 
     * 
     * @throws IOException 
     */
    public void prepare_pipeline() throws IOException, Exception {      
        
        
        
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
     * This preprocesses the steps to ensure as far as possible that the 
     * pipeline is runnable. 
     * @throws IOException, Exception
     */
    public void prepareSteps() throws IOException, Exception{
        
        for (NGSRunStepData stepData: this.getPipelineData().getStepsData()){
            
            switch (stepData.getStepType()){
                
                case StepUnzipInputFiles.STEP_ID_STRING:
                    this.addStepUnzipInputFiles(stepData);                    
                    break;
                    
                case StepSingleTrimAdapters.STEP_ID_STRING:
                    this.addStepSingleTrimAdapters(stepData);
                    break;
                    
                case StepCollapseReads.STEP_ID_STRING:                
                    this.addStepCollapseReads(stepData);
                    break;
                    
                case StepBowtieMapSingleReads.STEP_ID_STRING:
                    this.addStepBowtieMapSingleReads(stepData);
                    break;
                    
                case "BowtieMapPairedReads":
                    break;
                    
                case StepParseSAMForMiRNAs.STEP_ID_STRING:
                    this.addStepParseSAMForMiRNAs(stepData);
                    break;
                    
                case StepAnalyzeSAMforStartPositions.STEP_ID_STRING:
                    this.addStepAnalyzeStartPositions(stepData);
                    break;
                    
                case StepAnalyzeIsomiRDispersion.STEP_ID_STRING:
                    this.addStepAnalyzeIsomiRDispersion(stepData);
                    break;

                case StepDEwithEdgeR.STEP_ID_STRING:
                    this.addStepDifferentialExpression(stepData);
                    break;
                    
                case StepCleanUp.STEP_ID_STRING:
                    this.addStepCleanup(stepData);
                    break;
                    
                case "exit":
                    return;
                    
            }
            
        }
        
    }

    
    
    public void executePipeline() throws IOException, Exception{
        for (NGSBase ngsStep: ngsSteps){
            ngsStep.execute();
        }
        
    }
    
    
    /**
     * verify and add step to unzip input files for pipeline
     * 
     * @param stepData 
     */
    private void addStepUnzipInputFiles(NGSRunStepData stepData) throws IOException, Exception{
        
        StepInputData sidUnzip = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepUnzipInputFiles ngsUnzipFastqStep = new StepUnzipInputFiles(sidUnzip);        
        ngsUnzipFastqStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepUnzipInputFiles.STEP_ID_STRING));   
        ngsUnzipFastqStep.verifyInputData();
        ngsSteps.add(ngsUnzipFastqStep);
        
    }

    
    
    /**
     * build and execute step to trim adapter sequences from single end FASTQ files
     * @param stepData 
     */
    private void addStepSingleTrimAdapters(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidTrim = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepSingleTrimAdapters ngsSingleTrimStep = new StepSingleTrimAdapters(sidTrim);
        ngsSingleTrimStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepSingleTrimAdapters.STEP_ID_STRING));        
        ngsSingleTrimStep.verifyInputData();

        
    }
    
    
    
    /**
     * build and execute step to collapse supplied FASTQ files into FASTA files
     * and then identify identical reads
     * 
     * @param stepData 
     */
    private void addStepCollapseReads(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidCollapse = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepCollapseReads ngsCollapseStep = new StepCollapseReads(sidCollapse);
        ngsCollapseStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepCollapseReads.STEP_ID_STRING));
        ngsCollapseStep.verifyInputData();

        
    }
    
    
    /**
     * build and execute step to map single reads using Bowtie
     * 
     * @param stepData 
     */
    private void addStepBowtieMapSingleReads(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidMapSR = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapSingleReads ngsBowtieMapReads = new StepBowtieMapSingleReads(sidMapSR);
        ngsBowtieMapReads.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepBowtieMapSingleReads.STEP_ID_STRING));
        ngsBowtieMapReads.verifyInputData();
        
    }
    
    
    /**
     * build and execute step to map paired end reads using bowtie
     * 
     * @param stepData 
     */
    private void addStepBowtieMapPairedReads(NGSRunStepData stepData) throws IOException, Exception{


        StepInputData sidMapPR = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapPairedReads ngsBowtieMapPairedReads = new StepBowtieMapPairedReads(sidMapPR);
        ngsBowtieMapPairedReads.verifyInputData();
        
    }
    
    
    /**
     * build and execute step to parse a SAM alignment file to identify miRNAs
     * 
     * @param stepData 
     */
    private void addStepParseSAMForMiRNAs(NGSRunStepData stepData) throws IOException, Exception{
        

       StepInputData sidSAM = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
       StepParseSAMForMiRNAs ngsParseSAMForMiRNAs = new StepParseSAMForMiRNAs(sidSAM);
       ngsParseSAMForMiRNAs.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepParseSAMForMiRNAs.STEP_ID_STRING));
       ngsParseSAMForMiRNAs.verifyInputData();
       
    }
    
    
    /**
     * build and execute step to analyze start reads in SAM file to identify
     * reads that may be associated with novel smallRNAs
     * 
     * @param stepData 
     */
    private void addStepAnalyzeStartPositions(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidStart = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeSAMforStartPositions ngsAnalyzeSAMStartPos = new StepAnalyzeSAMforStartPositions(sidStart);
        ngsAnalyzeSAMStartPos.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepAnalyzeSAMforStartPositions.STEP_ID_STRING));
        ngsAnalyzeSAMStartPos.verifyInputData();
        
    }
    
    
    /**
     * build and execute step to analyze SAM file to identify isomiR populations
     * 
     * @param stepData 
     */
    private void addStepAnalyzeIsomiRDispersion(NGSRunStepData stepData) throws IOException, Exception{
        

        StepInputData sidIsoDisp = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeIsomiRDispersion analyzeIsomiRDispersions = new StepAnalyzeIsomiRDispersion(sidIsoDisp);
        analyzeIsomiRDispersions.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepAnalyzeIsomiRDispersion.STEP_ID_STRING));
        analyzeIsomiRDispersions.verifyInputData();
        
    }
    
    
    /**
     * build and execute step to perform differential expression based on 
     * supplied count files and grouping information
     * 
     * @param stepData 
     */
    private void addStepDifferentialExpression(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData siodDiffExpr = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepDEwithEdgeR edgeRDE = new StepDEwithEdgeR(siodDiffExpr);

        edgeRDE.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepDEwithEdgeR.STEP_ID_STRING));
        edgeRDE.verifyInputData();
        
    }
    
    
    /**
     * build and execute step to perform clean up after an analysis is complete
     * 
     * @param stepData 
     */
    private void addStepCleanup(NGSRunStepData stepData) throws IOException, Exception{


        StepInputData sidCleanUp = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepCleanUp cleanUp = new StepCleanUp(sidCleanUp);
        cleanUp.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepCleanUp.STEP_ID_STRING));
        cleanUp.verifyInputData();
        
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
        refDataLocations = new ReferenceDataLocations();
        refDataLocations.setGenomeRootFolder((String) dataOptions.get(ReferenceDataLocations.ID_GENOME_FOLDER));
        refDataLocations.setMirbaseFolder((String)dataOptions.get(ReferenceDataLocations.ID_MIRBASE_FOLDER));
        

        HashMap softwareOptions = (HashMap) pipelineConfigurationDataHash.get("software");
        
                
        
        
        logger.info("done\n");
        
    }
    
    
    /**
     * generate a sample configuration file. this is primarily to inform
     * users what can be specified for each step
     * @throws IOException 
     */
    public void generateSampleConfigurationFile() throws IOException{
        
        logger.info("write example pipeline configuration file");
        DumperOptions options = new DumperOptions();
        options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);
        yaml = new Yaml(options);   

        
        StepInputData emptySID = new StepInputData("", "", new ReferenceDataLocations(), "", "", new ArrayList<SampleDataEntry>());
        Map <String, Object> pipelineExampleConfiguration = new HashMap();
        StepUnzipInputFiles stepUnzip = new StepUnzipInputFiles(emptySID);
        pipelineExampleConfiguration.put(StepUnzipInputFiles.STEP_ID_STRING, stepUnzip.generateExampleConfigurationData());
        
        StepSingleTrimAdapters stepSingleAdapterTrim = new StepSingleTrimAdapters(emptySID);
        pipelineExampleConfiguration.put(StepSingleTrimAdapters.STEP_ID_STRING, stepSingleAdapterTrim.generateExampleConfigurationData());
        
        StepCollapseReads stepcollapseReads = new StepCollapseReads(emptySID);
        pipelineExampleConfiguration.put(StepCollapseReads.STEP_ID_STRING, stepcollapseReads.generateExampleConfigurationData());
        
        StepBowtieMapSingleReads stepBowtieSingleMap = new StepBowtieMapSingleReads(emptySID);
        pipelineExampleConfiguration.put(StepBowtieMapSingleReads.STEP_ID_STRING, stepBowtieSingleMap.generateExampleConfigurationData());
        
        StepAnalyzeSAMforStartPositions stepAnalyzeStartPos = new StepAnalyzeSAMforStartPositions(emptySID);
        pipelineExampleConfiguration.put(StepAnalyzeSAMforStartPositions.STEP_ID_STRING, stepAnalyzeStartPos.generateExampleConfigurationData());
        
        StepParseSAMForMiRNAs stepParseSAMforMirs = new StepParseSAMForMiRNAs(emptySID);
        pipelineExampleConfiguration.put(StepParseSAMForMiRNAs.STEP_ID_STRING, stepParseSAMforMirs.generateExampleConfigurationData());
        
        StepDEwithEdgeR stepDEwithEdgeR = new StepDEwithEdgeR(emptySID);
        pipelineExampleConfiguration.put(StepDEwithEdgeR.STEP_ID_STRING, stepDEwithEdgeR.generateExampleConfigurationData());
        
        StepCleanUp stepCleanUp = new StepCleanUp(emptySID);
        pipelineExampleConfiguration.put(StepCleanUp.STEP_ID_STRING, stepCleanUp.generateExampleConfigurationData());
        
        

        String dumpFilename = System.getProperty("user.dir") + FILE_SEPARATOR + "pipelineConfiguration.sample.yaml";
        StringWriter writer = new StringWriter();
        yaml.dump(pipelineExampleConfiguration, writer);
        logger.info(writer);
        try (FileWriter sampleFileWriter = new FileWriter(new File(dumpFilename))) {
            yaml.dump(pipelineExampleConfiguration, sampleFileWriter);    
        }
    }
    
    
    
    /**
     * read pipeline file
     * 
     * @throws Exception 
     */
    public void readPipelineFile() throws Exception{
        
	ObjectMapper mapper = new ObjectMapper();
        mapper.setVisibility(JsonMethod.FIELD, JsonAutoDetect.Visibility.ANY);

        logger.info("read pipeline file <" + this.getPipelineFile() + ">");
        pipelineData = mapper.readValue(new File(this.getPipelineFile()), PipelineData.class); 
        logger.info("done");
        logger.info(getPipelineData().toString());
                        
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
     * @return the SampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return SampleData;
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
     * @return the dataLocations
     */
    public ReferenceDataLocations getDataLocations() {
        return refDataLocations;
    }

    /**
     * @param dataLocations the dataLocations to set
     */
    public void setDataLocations(ReferenceDataLocations dataLocations) {
        this.refDataLocations = dataLocations;
    }

    /**
     * @return the generateExampleConfigurationFile
     */
    public Boolean getGenerateExampleConfigurationFile() {
        return generateExampleConfigurationFile;
    }

    /**
     * @param generateExampleConfigurationFile the generateExampleConfigurationFile to set
     */
    public void setGenerateExampleConfigurationFile(Boolean generateExampleConfigurationFile) {
        this.generateExampleConfigurationFile = generateExampleConfigurationFile;
    }

    
    
    
    
}
