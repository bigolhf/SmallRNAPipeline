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
    private ReferenceDataLocations      refDataLocations;
    
    
    
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
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepUnzipInputFiles ngsUnzipFastqStep = new StepUnzipInputFiles(sidUnzip);        
        ngsUnzipFastqStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepUnzipInputFiles.STEP_ID_STRING));        
        ngsUnzipFastqStep.execute();
        
    }

    
    
    /**
     * build and execute step to trim adapter sequences from single end FASTQ files
     * @param stepData 
     */
    private void executeStepSingleTrimAdapters(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidTrim = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
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
    private void executeStepCollapseReads(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidCollapse = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepCollapseReads ngsCollapseStep = new StepCollapseReads(sidCollapse);
        ngsCollapseStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(StepCollapseReads.STEP_ID_STRING));

        ngsCollapseStep.execute();
        
    }
    
    
    /**
     * build and execute step to map single reads using Bowtie
     * 
     * @param stepData 
     */
    private void executeStepBowtieMapSingleReads(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidMapSR = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapSingleReads ngsBowtieMapReads = new StepBowtieMapSingleReads(sidMapSR);
        ngsBowtieMapReads.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepBowtieMapSingleReads.STEP_ID_STRING));
        ngsBowtieMapReads.execute();
        
    }
    
    
    /**
     * build and execute step to map paired end reads using bowtie
     * 
     * @param stepData 
     */
    private void executeStepBowtieMapPairedReads(NGSRunStepData stepData) throws IOException, Exception{


        StepInputData sidMapPR = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepBowtieMapPairedReads ngsBowtieMapPairedReads = new StepBowtieMapPairedReads(sidMapPR);

        ngsBowtieMapPairedReads.execute();
        
    }
    
    
    /**
     * build and execute step to parse a SAM alignment file to identify miRNAs
     * 
     * @param stepData 
     */
    private void executeStepParseSAMForMiRNAs(NGSRunStepData stepData) throws IOException, Exception{
        

       StepInputData sidSAM = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
       StepParseSAMForMiRNAs ngsParseSAMForMiRNAs = new StepParseSAMForMiRNAs(sidSAM);
       ngsParseSAMForMiRNAs.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepParseSAMForMiRNAs.STEP_ID_STRING));

       ngsParseSAMForMiRNAs.execute();
       
    }
    
    
    /**
     * build and execute step to analyze start reads in SAM file to identify
     * reads that may be associated with novel smallRNAs
     * 
     * @param stepData 
     */
    private void executeStepAnalyzeStartPositions(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData sidStart = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeMappedReads ngsAnalyzeSAMStartPos = new StepAnalyzeMappedReads(sidStart);
        ngsAnalyzeSAMStartPos.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepAnalyzeMappedReads.STEP_ID_STRING));
        ngsAnalyzeSAMStartPos.execute();
        
    }
    
    
    /**
     * build and execute step to analyze SAM file to identify isomiR populations
     * 
     * @param stepData 
     */
    private void executeStepAnalyzeIsomiRDispersion(NGSRunStepData stepData) throws IOException, Exception{
        

        StepInputData sidIsoDisp = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepAnalyzeIsomiRDispersion analyzeIsomiRDispersions = new StepAnalyzeIsomiRDispersion(sidIsoDisp);
        analyzeIsomiRDispersions.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepAnalyzeIsomiRDispersion.STEP_ID_STRING));
        analyzeIsomiRDispersions.execute();
        
    }
    
    
    /**
     * build and execute step to perform differential expression based on 
     * supplied count files and grouping information
     * 
     * @param stepData 
     */
    private void executeStepDifferentialExpression(NGSRunStepData stepData) throws IOException, Exception{

        StepInputData siodDiffExpr = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
        StepDEwithEdgeR edgeRDE = new StepDEwithEdgeR(siodDiffExpr);

        edgeRDE.parseConfigurationData((HashMap)pipelineConfigurationDataHash.get(StepDEwithEdgeR.STEP_ID_STRING));
        edgeRDE.execute();
        
    }
    
    
    /**
     * build and execute step to perform clean up after an analysis is complete
     * 
     * @param stepData 
     */
    private void executeStepCleanup(NGSRunStepData stepData) throws IOException, Exception{


        StepInputData sidCleanUp = new StepInputData(this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), 
                 refDataLocations, stepData.getInputFileList(), stepData.getOutputFileList(), this.getSampleData());
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
        refDataLocations = new ReferenceDataLocations();
        refDataLocations.setGenomeRootFolder((String) dataOptions.get(ReferenceDataLocations.ID_GENOME_FOLDER));
        refDataLocations.setMirbaseFolder((String)dataOptions.get(ReferenceDataLocations.ID_MIRBASE_FOLDER));
        

        HashMap softwareOptions = (HashMap) pipelineConfigurationDataHash.get("software");
        
                
        
        
        logger.info("done\n");
        
    }
    
    
/*
    Map<String, Object> data = new HashMap<String, Object>();
    data.put("name", "Silenthand Olleander");
    data.put("race", "Human");
    data.put("traits", new String[] { "ONE_HAND", "ONE_EYE" });
    Yaml yaml = new Yaml();
    String output = yaml.dump(data);
    System.out.println(output);        


        */        
    
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
     * @return the NGSSteps
     */
    public ArrayList<NGSStep> getNGSSteps() {
        return NGSSteps;
    }


    /**
     * @return the SampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return SampleData;
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

    
    
    
    
}
