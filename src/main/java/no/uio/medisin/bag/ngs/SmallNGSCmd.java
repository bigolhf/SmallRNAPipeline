/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngs;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author sr
 */
public class SmallNGSCmd {
    
    static Logger logger = LogManager.getRootLogger();
    static Options options = new Options();
    SmallNGSPipeline smallPipeline;
    
    public void main(String args[]){
        
        smallPipeline = new SmallNGSPipeline();
        parseArguments(args);
        
        
    }
    
    
    
    /**
     * 
     * get pipeline settings
     * User must provide the following
     *   run config file    : all the information needed to run the pipeline (e.g. data folders)
     *   data file          : contains information about samples and grouping
     *   steps file         : contains information about what type of analyses should be performed
     * 
     * @param args 
     * 
     */
    public void parseArguments(String args[]){
        
        
        logger.info("parse arguments");
        options.addOption("h", "help", false, "view help");
        options.addOption("c", "config", true, "pipeline configuration file");
        options.addOption("d", "data", true, "data file list with grouping");
        options.addOption("s", "steps", true, "step file in JSON format");
        
        CommandLineParser clParser = new BasicParser();
        CommandLine cmd = null;
        
        
        try{
            cmd = clParser.parse(options, args);
                
            if(cmd.hasOption("h")){
                this.printHelp();
            }
            
            if (cmd.hasOption("c")) {
                logger.info("run configuration file set to " + cmd.getOptionValue("c"));
                smallPipeline.setConfigurationFile(cmd.getOptionValue("c"));
            }                    
            else
                throw new ParseException("no configuration file was specified") ;
            
            if (cmd.hasOption("p")) {
                logger.info("pipeline file set to " + cmd.getOptionValue("p"));
                smallPipeline.setPipelineFile(cmd.getOptionValue("p"));
            }                    
            else
                throw new ParseException("no pipeline file was specified") ;
            
            if (cmd.hasOption("d")) {
                logger.info("NGS data filelist set to " + cmd.getOptionValue("d"));
                smallPipeline.setDataFile(cmd.getOptionValue("d"));
            }                    
            else
                throw new ParseException("no NGS data filelist was specified") ;
            
        }
        catch(ParseException exPa){
            
        }
    }
    

    
    
    /**
     * print command line usage
     * 
     */
    public void printHelp(){
        printBanner();
        HelpFormatter formatter = new HelpFormatter();

        formatter.printHelp("command line options", options);
        
        System.exit(0);
        
    }
    
    
    
    /**
     * print program info
     * 
     */
    public void printBanner(){
        logger.info("\n\n\n"
                + "    =======================================================================\n"
                + "    |  NGS smallRNA pipeline:                                             |"
                + "    |    performs a range of analyses on smallRNA NGS data                |\n"
                + "    =======================================================================\n\n\n");
        
        logger.info("*** report bugs to simon.rayner@medisin.uio.no\n");
    }
    
    
}
