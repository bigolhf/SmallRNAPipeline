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
     * @return list size
     * @throws IOException 
     */
    public int loadConservedFamilyList(String miRFamilyDataFilename) throws IOException{
        targetScanMirFamilies = new ArrayList<>();
        
        String familyLine = "";
        try(BufferedReader brCF = new BufferedReader(new FileReader(new File(miRFamilyDataFilename)))){
            brCF.readLine();
            while((familyLine=brCF.readLine())!=null){
                logger.debug(familyLine);
                targetScanMirFamilies.add(new TargetScanMirFamily(familyLine));
            }
            brCF.close();   
            return targetScanMirFamilies.size();
        }
        catch(IOException exIO){
            logger.error("error reading <miR_Family_Info> file < " + miRFamilyDataFilename);
            logger.error("error occurred on this line");
            logger.error(familyLine);
            
            throw new IOException("error reading " + "<miR_Family_Info> file < " + miRFamilyDataFilename + ">\n"
                    + "error occurred on this line\n" + familyLine);
        }
            
    }
    
    
    
    /**
     * get a list of all the family entries matching the query seed
     * 
     * @param qSeedm8
     * @return 
     */
    public ArrayList<String> findSeedHits(String qSeedm8){
        ArrayList<String> familyHitList = new ArrayList<>();
        for(TargetScanMirFamily familyEntry: this.targetScanMirFamilies){
            if(familyEntry.getSeedm8().equals(qSeedm8)){
                if(!familyHitList.contains(familyEntry.getMiRFamily())){
                    familyHitList.add(familyEntry.getMiRFamily());
                }
            }
        }
        
        return familyHitList;
    }
    

        
        
        
     /**
     * return the number of entries in the dataset
     * 
     * @return 
     */
    public int getNumberOfEntries(){
        return targetScanMirFamilies.size();
    }
    
            
    
}
