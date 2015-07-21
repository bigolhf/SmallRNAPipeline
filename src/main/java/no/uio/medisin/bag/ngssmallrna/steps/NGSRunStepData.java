/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

/**
 *
 * @author sr
 */
public class NGSRunStepData {
    
    private String stepType;
    private String inputFileList;
    private String outputFileList;

    
    @Override
    	public String toString() {
            return "stepType=" + stepType + "\t in:" + inputFileList + "\t out:" + outputFileList;
        }
    /**
     * @return the name
     */
    public String getStepType() {
        return stepType;
    }

    /**
     * @param name the name to set
     */
    public void setStepType(String name) {
        this.stepType = name;
    }

    /**
     * @return the inputFileList
     */
    public String getInputFileList() {
        return inputFileList;
    }

    /**
     * @param inputFileList the inputFileList to set
     */
    public void setInputFileList(String inputFileList) {
        this.inputFileList = inputFileList;
    }

    /**
     * @return the outputFileList
     */
    public String getOutputFileList() {
        return outputFileList;
    }

    /**
     * @param outputFileList the outputFileList to set
     */
    public void setOutputFileList(String outputFileList) {
        this.outputFileList = outputFileList;
    }
    
}
