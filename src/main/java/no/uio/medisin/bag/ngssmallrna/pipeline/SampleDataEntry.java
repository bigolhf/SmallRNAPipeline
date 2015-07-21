/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 *
 * @author sr
 */
public class SampleDataEntry {
    
    private String dataFile;
    private String dataSource;
    private String condition;
    private String time;

    
    
    public SampleDataEntry(String fname, String source, String cond, String t){
        
        dataFile = fname;
        dataSource = source;
        condition = cond;
        time = t;
        
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
     * @return the dataSource
     */
    public String getDataSource() {
        return dataSource;
    }

    /**
     * @param dataSource the dataSource to set
     */
    public void setDataSource(String dataSource) {
        this.dataSource = dataSource;
    }

    /**
     * @return the Condition
     */
    public String getCondition() {
        return condition;
    }

    /**
     * @param Condition the Condition to set
     */
    public void setCondition(String Condition) {
        this.condition = Condition;
    }

    /**
     * @return the time
     */
    public String getTime() {
        return time;
    }

    /**
     * @param time the time to set
     */
    public void setTime(String time) {
        this.time = time;
    }
}
