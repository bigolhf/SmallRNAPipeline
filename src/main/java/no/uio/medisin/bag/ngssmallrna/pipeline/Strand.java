/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 *
 * @author simonray
 */
public enum Strand {
    PLUS ("+"), 
    MINUS ("-"), 
    UNKNOWN("?");
    
    private final String strand;
    
    Strand(String strand){
        this.strand = strand;
    }

    @Override
    public String toString(){
        return strand;
    }
}


