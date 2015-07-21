/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.util.ArrayList;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;

/**
 *
 * @author sr
 */
abstract public class NGSStep {
    
    abstract void verifyInputData();
    abstract void outputResultData();
    abstract void execute();
    
}
