/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import java.util.List;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepData;

/**
 *
 * @author sr
 */
public class PipelineData {
 
        private String pipelineName;
        private String projectID;
        private String projectRoot;
        private String dataRoot;
        
        private final List <NGSRunStepData> stepsData = new ArrayList<>();
        
 

        
        
        
 
	@Override
	public String toString() {
            String str = "User [" + pipelineName + "]\n";
                for(NGSRunStepData step: getStepsData()){
                    str += step.toString() + "\n";
                }
		return str;
	}

    /**
     * @return the steps
     */
    public List <NGSRunStepData> getStepsData() {
        return stepsData;
    }

    /**
     * @return the projectID
     */
    public String getProjectID() {
        return projectID;
    }

    /**
     * @return the projectRoot
     */
    public String getProjectRoot() {
        return projectRoot;
    }

    /**
     * @return the dataRoot
     */
    public String getDataRoot() {
        return dataRoot;
    }
}
