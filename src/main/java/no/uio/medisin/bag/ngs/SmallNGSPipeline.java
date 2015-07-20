/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import static no.uio.medisin.bag.ngs.SmallNGSCmd.logger;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStepData;
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
    
    static  Yaml                    yaml = new Yaml();
    static  String                  FileSeparator = System.getProperty("file.separator");
    
    private String                  configurationFile = "";
    private String                  pipelineFile = "";
    private String                  dataFile = "";
    
    private String                  softwareRootFolder = "";
    private String                  adapterTrimmingSoftware = "";
    
    private PipelineData            pipelineData = new PipelineData();
    private ArrayList<NGSDataEntry> NGSData = new ArrayList<>();
    private ArrayList<NGSStep>      NGSSteps = new ArrayList<>();

    
    
    /**
     * here we will read the steps into a Hash and then step through
     * specific checks are carried out within each step
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
    
    
    
    public void buildPipeline(){
        
        for (NGSStepData stepData: this.getPipelineData().getSteps()){
            
            getNGSSteps().add(null);
        }
        
    }

    
    public void runPipeline(){
        
        
        for (NGSStep ngsStep: this.getNGSSteps()){
            
        }
        
    }
    
    
    /**
     * read file containing information about the data to be analyzed
     * 
     * @throws IOException 
     * 
     */
    public void readDataFile() throws IOException{
        
        
        logger.info("read data file");
        String line = "";
        BufferedReader bwData = new BufferedReader(new FileReader(new File(this.getDataFile())));
            while((line=bwData.readLine())!= null){
                
                String tokens[] = line.split("\t");
                String file = tokens[0];
                String source = tokens[1];
                String condition = tokens[2];
                String time = tokens[3];
                
                NGSData.add(new NGSDataEntry(file, source, condition, time));
                
            }
        bwData.close();
        logger.info("read " + NGSData.size() + " entries");
        
        
    }
    
    
    
    
    /**
     * read configuration settings for the pipeline
     * 
     * @throws IOException 
     */
    public void readConfigurationFile() throws IOException{
        
        logger.info("read pipeline configuration file");
        HashMap pipelineConfiguration = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));

        HashMap softwareOptions = (HashMap) pipelineConfiguration.get("software");
        this.setSoftwareRootFolder((String) softwareOptions.get("root_folder"));
        this.setAdapterTrimmingSoftware((String) softwareOptions.get("adapter_trimming"));
                        
        logger.info("done\n");
        
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
    
    
    
    
}
