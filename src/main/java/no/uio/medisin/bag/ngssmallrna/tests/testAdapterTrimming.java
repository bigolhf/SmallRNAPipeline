/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.tests;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import no.uio.medisin.bag.ngssmallrna.steps.StepInputData;
import no.uio.medisin.bag.ngssmallrna.steps.TrimAdapters;


/**
 * test Adapter Trimming step
 * @author sr
 */
public class testAdapterTrimming {
    
    static Logger logger = LogManager.getRootLogger();
    public static void main(String args[]){
        StepInputData sid = new StepInputData();
        sid.setProjectRoot("/home/sr/research/");
        sid.setProjectID("sweden");
        TrimAdapters trimAdapters = new TrimAdapters(sid);
        trimAdapters.performTrimming();
        
    }
}
