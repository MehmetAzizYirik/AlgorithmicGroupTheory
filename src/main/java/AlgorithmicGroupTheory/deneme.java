package AlgorithmicGroupTheory;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import com.google.common.io.Files;

public class deneme {
	public static IChemObjectBuilder builder =  SilentChemObjectBuilder.getInstance();
	public static IAtomContainer build(int size) {
		IAtomContainer ac = builder.newInstance(IAtomContainer.class); 
	    for(int i=0;i<size;i++) {
	    	ac.addAtom(builder.newInstance(IAtom.class, "C"));
	    }
	    return ac;
	}
	
	/**
	 * Molecule depiction 
	 */
	
	public static void depict(IAtomContainer molecule, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		//depiction.withSize(200, 200).withZoom(4).depict(molecule).writeTo(path);
		//depiction.withCarbonSymbols().withAtomValues().withAtomMapNumbers().withAtomColors().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
		depiction.withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
		//depiction.withAtomColors().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
	}
	public static IAtomContainer addSingleBond(IAtomContainer ac, String[] indices) {
		for(int i=0;i<indices.length;i++) {
			if(i%2!=0) {
				ac.addBond(((int)Integer.parseInt(indices[i-1])), (int)Integer.parseInt(indices[i]), Order.SINGLE);
			}
		}
		return ac;
	}
	
	public static void main(String[] args)throws Exception { 
	// We need to provide file path as the parameter: 
	// double backquote is to avoid compiler interpret words 
	// like \test as \t (ie. as a escape sequence) 
		
	FileReader fr=new FileReader("C:\\Users\\mehme\\Desktop\\output2.txt");   //reads the file  
	BufferedReader br=new BufferedReader(fr);  //creates a buffering character input stream  
	StringBuffer sb=new StringBuffer();    //constructs a string buffer with no characters  
	String line; 
	String[] say= new String[4];
	say[0]="1";
	say[1]="2";
	say[2]="1";
	say[3]="3";
	
	int n=0;
	for(int i=0;i<say.length;i++) {
		if(i%2!=0) {
			n=(int)Integer.parseInt(say[i]);			
		}
	}
	
	int i=0;
	while((line=br.readLine())!=null)  
	{  
		if(line.length()>21) {
			
			//System.out.println(line.length()+" "+line); 
			String[] splited = line.split("\\s+");
			IAtomContainer ac= build(10);
			ac=addSingleBond(ac,splited);
			i++;
			depict(ac,"C:\\Users\\mehme\\Desktop\\deneme\\molecule "+i+".png");
			//System.out.println(Arrays.toString(splited));
		}
	}  
	fr.close();    //closes the stream and release the resources  
	}
}
