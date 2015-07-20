/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngs;

import java.util.ArrayList;
import java.util.List;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStepData;

/**
 *
 * @author sr
 */
public class PipelineData {
 
        private String pipelineName;
        private final List <NGSStepData> steps = new ArrayList<>();
        
 

        
        
        
 
	@Override
	public String toString() {
            String str = "User [" + pipelineName + "]\n";
                for(NGSStepData step: getSteps()){
                    str += step.toString() + "\n";
                }
		return str;
	}

    /**
     * @return the steps
     */
    public List <NGSStepData> getSteps() {
        return steps;
    }
}
