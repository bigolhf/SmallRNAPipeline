/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author simonray
 */
public class TargetScanMirFamilyList {
    
    static Logger                       logger                          = LogManager.getLogger();
    private ArrayList<TargetScanMirFamily> targetScanMirFamilies;
    
    /**
     * 
     * @param miRFamilyDataFilename
     * @throws IOException 
     */
    public void loadConservedFamilyList(String miRFamilyDataFilename) throws IOException{
        targetScanMirFamilies = new ArrayList<>();
        
        String familyLine = "";
        try(BufferedReader brCF = new BufferedReader(new FileReader(new File(miRFamilyDataFilename)))){
            while((familyLine=brCF.readLine())!=null){
                targetScanMirFamilies.add(new TargetScanMirFamily(familyLine));
            }
            brCF.close();    
        }
        catch(IOException exIO){
            logger.error("error reading <miR_Family_Info> file < " + miRFamilyDataFilename);
            logger.error("error occurred on this line");
            logger.error(familyLine);
            
            throw new IOException("error reading " + "<miR_Family_Info> file < " + miRFamilyDataFilename + ">\n"
                    + "error occurred on this line\n" + familyLine);
        }
            
        
    }
    
    
            
    
}
