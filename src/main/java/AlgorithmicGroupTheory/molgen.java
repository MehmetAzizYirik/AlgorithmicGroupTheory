/**
 * MIT License
 *
 * Copyright (c) 2019 Mehmet Aziz Yirik
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


/**
 * This class is the re-implementation of MOLGEN project. For more details;
 * please check Mathematical Chemistry and Cheminformatics book[1]. Chapter 1
 * and Chapter 5 for molecular structure generation.
 * 
 * [1] Kerber, A., Laue, R., Meringer, M., Rücker, C. and Schymanski, E., 2013.
 * Mathematical chemistry and chemoinformatics: structure generation, elucidation
 * and quantitative structure-property relationships. Walter de Gruyter.
 * 
 * @author Mehmet Aziz Yirik
 */

package AlgorithmicGroupTheory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.silent.Atom;


public class molgen {
	public static boolean verbose = false;
	public static String formula;
	static String filedir = null;
	static BufferedWriter fileWriter;
	public static Map<String, Integer> valences; 
	static {
		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
			
		valences.put("C", 4);
		valences.put("N", 5);
		valences.put("O", 2);
		valences.put("S", 6);
		valences.put("P", 5);
		valences.put("F", 1);
		valences.put("I", 1);
		valences.put("Cl", 1);
		valences.put("Br", 1);
		valences.put("H", 1);
	}
	
	public static IAtomContainer build(String mol) throws IOException, CloneNotSupportedException, CDKException {
		IAtomContainer atomcontainer = new AtomContainer();
		List<String> symbols = new ArrayList<String>();
	    List<Integer> occur = new ArrayList<Integer>();
	    String[] atoms = mol.split("(?=[A-Z])");
	    for (String atom : atoms) {
	    	String[] info = atom.split("(?=[0-9])", 2);   
	        symbols.add(info[0]);
	        occur.add(info.length > 1 ? Integer.parseInt(info[1]):1);
	    }
	    for(int i=0;i<symbols.size();i++) {
	    	for(int j=0;j<occur.get(i);j++) {
	    		atomcontainer.addAtom(new Atom(symbols.get(i)));
	        }
	    }
	    return atomcontainer;
	}
	
	public static ArrayList<Integer> degreeSeq(String formula){
		ArrayList<Integer> seq= new ArrayList<Integer>();
		List<String> symbols = new ArrayList<String>();
		List<Integer> occur  = new ArrayList<Integer>();
	    String[] atoms = formula.split("(?=[A-Z])");
	    for (String atom : atoms) {
	    	String[] info = atom.split("(?=[0-9])", 2);   
	        symbols.add(info[0]);
	        occur.add(info.length > 1 ? Integer.parseInt(info[1]):1);
	    }
	    for(int i=0;i<symbols.size();i++) {
	    	String symbol= symbols.get(i);
	        for(int j=0;j<occur.get(i);j++) {
	        	seq.add(valences.get(symbol));
	        }
	    }
		return seq;
	}
	
	 /** 
	  * Grund Thesis - Algorithm 3.2.3 (DONE)
	  * Simple multigraph generation
	  * @param degrees atom valences  
	  * @throws IOException
	  */
		 
	public static void canonicalMatrix(ArrayList<Integer> degrees) throws IOException{
		int size    = degrees.size();
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] L2  = upperTriangularL2(degrees,A);
		int[][] C   = upperTriangularC(degrees);
		int[][] C2  = upperTriangularC2(degrees,A);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forward(degrees,A,max,L,C,indices); //It was originally without partition entry.
	}
		
	//Set diagonal to zeros.
	public static void zeroDiagonal(int[][] mat) {
		for(int i=0; i<mat.length;i++) {
			for(int j=0;j<mat.length;j++) {
				if(i==j) {
					mat[i][j]=0;
				}
			}
		}
	}
		 
	/**
	* Forward step in Algorithm 3.2.3 (DONE)
	* @param degrees valence
	* @param A matrix in the generation process
	* @param max maximal matrix
	* @param L L Matrix	
	* @param C	C Matrix
	* @param indices initial indices (0,1)
	* @throws IOException
	*/
			
	public static void forward(ArrayList<Integer> degrees,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
		int minimal= Math.min(max[i][j],Math.min(l2,c2));
		for(int h=minimal;h>=0;h--) {
			if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
				A[i][j]=A[j][i]=h;
			 	if(i==(max.length-2) && j==(max.length-1)) {
			 		backward(degrees,A, max,L, C, indices);
			 	}else {
			 		ArrayList<Integer> modified=successor(indices,max.length);
			 		forward(degrees,A, max, L, C, modified);
			 	}
			 }else {
			 	backward(degrees,A, max, L, C, indices);
			 }
		 }
	}
		
	/**
	* Backward step in Algorithm 3.2.3 (DONE)
	* @param degrees valence
	* @param A matrix in the generation process
	* @param max maximal matrix
	* @param L L Matrix	
	* @param C	C Matrix
	* @param indices initial indices (0,1)
	* @throws IOException
	**/
	public static ArrayList<int[][]> output = new ArrayList<int[][]>();
	public static void backward(ArrayList<Integer> degrees,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
		if(i==max.length-2 && j==max.length-1) {
			int[][] mat2= new int[A.length][A.length];
			for(int k=0;k<A.length;k++) {
				for(int l=0;l<A.length;l++) {
					mat2[k][l]=A[k][l];
				}
			}
			output.add(mat2);
			writeMatrix(mat2);
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					ArrayList<Integer> modified2=successor(modified,max.length);
					forward(degrees,A, max, L, C, modified2);
				}else {
					backward(degrees, A, max, L, C, modified);
				}
			}
		}
	}
	
	/**
	* Grund Thesis 3.2.1 - Page 35 (DONE)
	* L and C Inverse formulas given at the bottom of the page.
	* @param degrees atom valences
	* @param i row index
	* @param j column index
	* @param A matrix in the generation process
	* @return
	*/
			
	public static int LInverse(ArrayList<Integer> degrees, int i, int j, int[][]A) {
		int sum=0;
		for(int s=0;s<j;s++) {
			sum=sum+A[i][s];
		}
		return degrees.get(i)-sum;
	}
				 
	public static int CInverse(ArrayList<Integer> degrees, int i, int j, int[][]A) {
		int sum=0;
		for(int s=0;s<i;s++) {
			sum=sum+A[s][j];
		}
		return degrees.get(j)-sum;
	}
			
	/**
	* Successor and predecessor functions for forward and backward functions. (DONE)
	* @param indices indices of the matrix
	* @param size matrix dimension
	* @return
	*/
			
	public static ArrayList<Integer> successor(ArrayList<Integer> indices, int size) {
		int i0= indices.get(0);
		int i1= indices.get(1);
		ArrayList<Integer> modified= new ArrayList<Integer>();
		//if(i0!=size-2 && i1!=size-1) {
		if(i1==size-1) {
			modified.add(i0+1);
			modified.add(i0+2);
		}else if(i1<size-1) {
			modified.add(i0);
			modified.add(i1+1);
		}
		//}
		return modified;
	}
			 
	public static ArrayList<Integer> predecessor(ArrayList<Integer> indices, int size) {
		int i0= indices.get(0);
		int i1= indices.get(1);
		ArrayList<Integer> modified= new ArrayList<Integer>();
		if(i0==i1-1) {
			modified.add(i0-1);
			modified.add(size-1);
		}else {
			modified.add(i0);
			modified.add(i1-1);
		}
		return modified;
	}
		
	//Page 35, the bottom upper triangular matrices.
	public static int Csum2(int[][] A, int i, int j) {
		int sum=0;
		for(int k=0;k<i;k++) {
			sum=sum+A[k][j];
		}
		return sum;
	}
		  
	/**
	* C; upper triangular matrix like given in 3.2.1.
	* @param degrees
	* @return upper triangular matrix
	*/
	public static int[][] upperTriangularC2(ArrayList<Integer> degrees, int[][]A){
		int size= degrees.size();
		int[][] C= new int[size][size]; 	
		for(int i=0;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				C[i][j]= (degrees.get(j-1) - Csum2(A,i,j));
			}
		}
		return C;
	}
		 
	//TODO: Thesis or my code is wrong. I modified by minus 1.
	/**
	* C; upper triangular matrix like given in 3.2.1.
	* @param degrees
	* @return upper triangular matrix
	*/
	
	public static int[][] upperTriangularC(ArrayList<Integer> degrees){
		int size= degrees.size();
		int[][] C= new int[size][size]; 	
		int[][] max= maximalMatrix(degrees);
		for(int i=0;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				C[i][j]= Math.min(degrees.get(j-1), Csum(max,i+1,j,size)); //TODO: WHY -1 I dont know.
			}
		}
		return C;
	}
		 	 
	public static int Csum(int[][] max, int i, int j, int size) {
		int sum=0;
		for(int k=i;k<size;k++) {
			sum=sum+max[k][j];
		}
		return sum;
	}
	
	/**
	* Theorem 3.2.1
	* @param degrees degree sequence
	*/
	
	public static int[][] maximalMatrix(ArrayList<Integer> degrees) {
		int[][] maxMatrix= new int[degrees.size()][degrees.size()];
		for(int i=0;i<degrees.size();i++) {
			for(int j=0; j<degrees.size();j++) {
				int di= degrees.get(i);
				int dj= degrees.get(j);
				if(i==j) {
					maxMatrix[i][j]=0;
				}else {
					if(di!=dj) {
						maxMatrix[i][j]=Math.min(di, dj);
					}else {
						maxMatrix[i][j]=di-1;
					}
				}
			}
		}
		return maxMatrix;
	}
		 
	/**
	* L; upper triangular matrix like given in 3.2.1.
	* @param degrees
	* @return upper triangular matrix
	*/
	public static int[][] upperTriangularL(ArrayList<Integer> degrees){
		int size= degrees.size();
		int[][] L= new int[size][size]; //TODO: Maybe zeros matrix ?	
		int[][] max= maximalMatrix(degrees);
		for(int i=0;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				L[i][j]= Math.min(degrees.get(i), Lsum(max,i,j+1,size));
			}
		}
		return L;
	}
		 
	public static int Lsum(int[][] max, int i, int j, int size) {
		int sum=0;
		for(int k=j;k<size;k++) {
			sum=sum+max[i][k];
		}
		return sum;
	}
		 
	/**
	* L; upper triangular matrix like given at the bottom of page 35.
	* @param degrees
	* @return upper triangular matrix
	*/
		 
	public static int[][] upperTriangularL2(ArrayList<Integer> degrees, int[][] A){	 
		int size= degrees.size();
		int[][] L= new int[size][size]; //TODO: Maybe zeros matrix ?	
		for(int i=0;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				L[i][j]= (degrees.get(i) - Lsum2(A,i,j));
			}
		}
		return L;
	}
		 
	//Page 35, the bottom upper triangular matrices.
	public static int Lsum2(int[][] A, int i, int j) {
		int sum=0;
		for(int k=0;k<j;k++) {
			sum=sum+A[i][k];
		}
		return sum;
	}
		 
	public static String inchiGeneration(IAtomContainer molecule) throws CDKException {
		String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(molecule).getInchi();	
		return inchi;
	}
		 
	public static void depict(IAtomContainer molecule, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		//depiction.withSize(200, 200).withZoom(4).depict(molecule).writeTo(path);
		depiction.withAtomColors().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
	}
		
	public static void run(String formula) throws IOException, CDKException, CloneNotSupportedException {
		if(verbose) System.out.println("For the formula, "+formula+",start generating all possible connectivity matrices...");
		canonicalMatrix(degreeSeq(formula));
		if(verbose) System.out.println("The number of connectivity matrices is: "+output.size());
		
	}
	
	public static void setFileWriter() throws IOException {
		 molgen.fileWriter = new BufferedWriter(new FileWriter(filedir+"output.txt"));
	}
	
	public static void writeMatrix(int[][] mat) throws IOException {
		 fileWriter.write(String.format("Connecticity matrix - %d", output.size()));
		 fileWriter.newLine();
		 for (int i = 0; i < mat.length; i++) {
			 for (int j = 0; j < mat[i].length; j++) {
				 fileWriter.write(mat[i][j] + ((j == mat[i].length-1) ? "" : ","));
			 }
			 fileWriter.newLine();
		 }
		 fileWriter.write("----------");
		 fileWriter.newLine();
		 fileWriter.flush();
	 }
	private void parseArgs(String[] args) throws ParseException, IOException
	{
		Options options = setupOptions(args);	
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			molgen.formula = cmd.getOptionValue("formula");
			molgen.filedir = cmd.getOptionValue("filedir");
			setFileWriter();
			if (cmd.hasOption("verbose")) molgen.verbose = true;		
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.setOptionComparator(null);
			String header = "\nGenerates connectivity matrices for a given molecular formula."
					+ " The input is a molecular formula string."
					+ "For example 'C2H2O'."
					+ "Besides this formula, the directory is needed to be specified for the output"
					+ "file. \n\n";
			String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory";
			formatter.printHelp( "java -jar AlgorithmicGroupTheory.jar", header, options, footer, true );
			throw new ParseException("Problem parsing command line");
		}
	}
		
	private Options setupOptions(String[] args)
	{
		Options options = new Options();
		Option formula = Option.builder("f")
				.required(true)
				.hasArg()
				.longOpt("formula")
				.desc("formula (required)")
				.build();
		options.addOption(formula);
		Option verbose = Option.builder("v")
				.required(false)
				.longOpt("verbose")
				.desc("print message")
				.build();
		options.addOption(verbose);	
		Option filedir = Option.builder("d")
			    .required(true)
			    .hasArg()
			    .longOpt("filedir")
			    .desc("Creates and store the output txt file in the directory (required)")
			    .build();
		options.addOption(filedir);
		return options;
	}
	
	public static void main(String[] args) throws CloneNotSupportedException, CDKException, IOException {
		molgen gen = null;
		//String[] args1= {"-f","C2OH2","-v","-d","C:\\Users\\mehme\\Desktop\\"};
		try {
			gen = new molgen();
			gen.parseArgs(args);
			molgen.run(molgen.formula);
		} catch (Exception e) {
			if (molgen.verbose) e.getCause(); 
		}
	}
}
