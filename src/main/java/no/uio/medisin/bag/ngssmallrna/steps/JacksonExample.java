/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import no.uio.medisin.bag.ngssmallrna.pipeline.PipelineData;
import java.io.File;
import java.io.IOException;
import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.annotate.JsonAutoDetect.Visibility;
import org.codehaus.jackson.annotate.JsonMethod;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

/**
 *
 * @author sr
 */
public class JacksonExample {
    public static void main(String[] args) {
 
	User user = new User();
        PipelineData user2 = new PipelineData();
	ObjectMapper mapper = new ObjectMapper();
        ObjectMapper mapper2 = new ObjectMapper();
        
        mapper.setVisibility(JsonMethod.FIELD, Visibility.ANY);
        mapper2.setVisibility(JsonMethod.FIELD, Visibility.ANY);
        
	try {
 
		// convert user object to json string, and save to a file
//		mapper.writeValue(new File("/home/sr/NetBeansProjects/NGSsmallRNA/user.json"), user);
                user = mapper.readValue(new File("/home/sr/NetBeansProjects/NGSsmallRNA/user.json"), User.class);
		// display to console
		System.out.println(mapper.writeValueAsString(user));
                
                user2 = mapper2.readValue(new File("/home/sr/NetBeansProjects/NGSsmallRNA/test/testpipe.json"), PipelineData.class); 
               
               System.out.println(user2);
 
	} catch (JsonGenerationException e) {
 
		e.printStackTrace();
 
	} catch (JsonMappingException e) {
 
		e.printStackTrace();
 
	} catch (IOException e) {
 
		e.printStackTrace();
 
	}
 
  }
 
}