package AlgorithmicGroupTheory;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.interfaces.IBond.Order;

import AlgorithmicGroupTheory.HydrogenDistributor;

public class buNe extends groupFunctions {
	public static int size=0;
	public static int hIndex=0;
	public static int count=0;
	public static boolean verbose = false;
	public static SDFWriter outFile;
	public static String formula;
	public static String filedir;
	public static ArrayList<int[][]> output= new ArrayList<int[][]>();
	public static ArrayList<String> inchis= new ArrayList<String>();
	public static ArrayList<Integer> initialDegrees= new ArrayList<Integer>();
	public static ArrayList<Integer> inputPartition= new ArrayList<Integer>();	
	public static IChemObjectBuilder builder=DefaultChemObjectBuilder.getInstance();
	public static IAtomContainer atomContainer= builder.newInstance(IAtomContainer.class);
	public static ArrayList<ArrayList<Integer>> partitionList= new ArrayList<ArrayList<Integer>>();
	public static ArrayList<ArrayList<Permutation>> representatives = new ArrayList<ArrayList<Permutation>>();
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
	
	/**
	 * InChI Generator from CDK.
	 * @param molecule IAtomContainer
	 * @return
	 * @throws CDKException
	 */
	
	public static String inchiGeneration(IAtomContainer molecule) throws CDKException {
		String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(molecule).getInchi();	
		return inchi;
	}
	
	/**
	 * The main function for the block based canonical matrix generation.
	 * @param degrees ArrayList<Integer> atom valences 
	 * @param partition ArrayList<Integer> atom partitioning
	 * @param filedir String The directory to store the output sdf file
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	
	public static void canonicalBlockbasedGenerator(ArrayList<Integer> degrees, ArrayList<Integer> partition, String filedir) throws IOException, CloneNotSupportedException, CDKException{
		long startTime = System.nanoTime(); 
		ArrayList<Integer> den= new ArrayList<Integer>();
		den.add(4);
		den.add(4);
		den.add(4);
		den.add(4);
		den.add(4);
		den.add(4);
		
		den.add(1);
		den.add(1);
		den.add(1);
		den.add(1);
		den.add(1);
		den.add(1);
		inputPartition=partition;
		initialDegrees=den;
		size=12;
		//List<ArrayList<Integer>> newDegrees= distributeHydrogens(partition, degrees);
		partitionList.add(0,inputPartition);
		gen(degrees);
		//for(ArrayList<Integer> degree: newDegrees) {
			//gen(degree);
		//}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		System.out.println("time"+" "+d.format(seconds)+" "+count);
	}
	
	/**
	 * Distribute the hydrogens in all possible ways to other atoms. 
	 * @param partition ArrayList<Integer> atom partitioning 
	 * @param degrees ArrayList<Integer> atom valences 
	 * @return
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	
	public static List<ArrayList<Integer>> distributeHydrogens(ArrayList<Integer> partition, ArrayList<Integer> degrees) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException{
		List<ArrayList<Integer>> degreeList= new ArrayList<ArrayList<Integer>>();
		List<int[]> distributions= HydrogenDistributor.run(partition,degrees);
		inputPartition=redefineThePartition(partition);
		hIndex= distributions.get(0).length;
		for(int[] dist: distributions) {
			ArrayList<Integer> newDegree= new ArrayList<Integer>();
			for(int i=0;i<dist.length;i++) {
				newDegree.add(i,(degrees.get(i)-dist[i]));
			}
			degreeList.add(newDegree);
		}
		return degreeList;
	}
	
	/**
	 * Once the hydrogens are distributed, we need to redefine the input partition.
	 * @param partition ArrayList<Integer> atom partitioning
	 * @return
	 */
	
	public static ArrayList<Integer> redefineThePartition(ArrayList<Integer> partition){
		ArrayList<Integer> newPart= new ArrayList<Integer>();
		for(int i=0;i<partition.size()-1;i++) {
			newPart.add(partition.get(i));
		}
		return newPart;
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
	
	public static int Csum(int[][] max, int i, int j, int size) {
		int sum=0;
		for(int k=i;k<size;k++) {
			sum=sum+max[k][j];
		}
		return sum;
	}
	 
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
	 
	/**
	 * Make diagonal entries zeros
	 * @param mat int[][]
	 */
	
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
	 * Generates structures for the degrees, set after hydrogen distribution.
	 * @param degrees ArrayList<Integer> degree
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	
	public static void gen(ArrayList<Integer> degrees) throws IOException, CloneNotSupportedException, CDKException {
		int r=0;
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] C   = upperTriangularC(degrees);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		forwardBuild(r,degrees,inputPartition,A,max,L,C,indices);
	}
	
	/**
	 * Functions to detect x,y,z row indices of a block as explained in Grund' thesis.
	 * @param r int block index
	 * @return
	 */
	
	public static int findX(int r) {
		return (sum(inputPartition,(r-1))-1);
	}
	
	public static int findY(int r) {
		return (sum(inputPartition,(r-1)));
	}
	
	public static int findZ(int r) {
		return (sum(inputPartition,r)-1);
	}
	
	/**
	 * Grund Thesis 3.2.1 - Page 35 
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
	 * Successor and predecessor functions for forward and backward functions.
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
	 
	 /**
	  * Forward function of the matrix generator.
	  * @param r int block index
	  * @param degrees ArrayList<Integer> valences
	  * @param partition ArrayList<Integer> atom partition
	  * @param A int[][] adjacency matrix
	  * @param max int[][] maximal matrix
	  * @param L int[][] L matrix
	  * @param C int[][] C matrix
	  * @param indices ArrayList<Integer> entry indices
	  * @throws IOException
	  * @throws CloneNotSupportedException
	  * @throws CDKException
	  */
	 
	 public static void forwardBuild(int r, ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException, CloneNotSupportedException, CDKException {
		 int y=findY(r);
		 int z=findZ(r);
		 int i=indices.get(0);
		 int j=indices.get(1);
		 int l2= LInverse(degrees,i,j,A);
		 int c2= CInverse(degrees,i,j,A);
		 int minimal= Math.min(max[i][j],Math.min(l2,c2));
		 for(int h=minimal;h>=0;h--) {
			 if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
				 A[i][j]=A[j][i]=h;	
				 if(i==(max.length-2) && j==(max.length-1)) {
					 if(blockTest(i, r, A)) {
						 backwardBuild(r, degrees,partitionList.get(y),A, max, L, C, indices);
			 		 }
			     }else {
			    	 ArrayList<Integer> modified=successor(indices,max.length);
			 		 if(modified.get(0)>z) { // If its higher than z value. 
			 			 if(blockTest(i, r, A)) {
			 				 r++; //z+1
			 				 forwardBuild(r, degrees, partitionList.get(z+1), A, max, L, C, modified);
			 			 }
			 		 }else {
			 			 forwardBuild(r, degrees, partitionList.get(y), A, max, L, C, modified);
			 		 }
			 	 }
			 	 }else {
			 		 if(i==max.length-2 && j==max.length-1) {
			 			 ArrayList<Integer> modified=predecessor(indices, max.length);
			 			 indices=modified;
			 		 }
			 		 backwardBuild(r, degrees, partitionList.get(y),A, max, L, C, indices);
			 	 }
		}	 
	}
	 
	/**public static void forwardBuild(int r, ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices, PrintWriter pWriter) throws IOException, CloneNotSupportedException, CDKException {
		int y=findY(r);
		int z=findZ(r);
		pWriter.print("forward build r y z"+" "+r+" "+y+" "+z+"\n");
		int i=indices.get(0);
		int j=indices.get(1);
		pWriter.print("forward indices"+" "+i+" "+j+"\n");
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
	 	int minimal= Math.min(max[i][j],Math.min(l2,c2));
	 	for(int h=minimal;h>=0;h--) {
	 		if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
	 			A[i][j]=A[j][i]=h;
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				pWriter.print("last entry"+"\n");
	 				backwardBuild(r, degrees,partitionList.get(y),A, max, L, C, indices, pWriter);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				pWriter.print("successor"+" "+modified.get(0)+" "+modified.get(1)+"\n");
	 				if(modified.get(0)>z) { // If its higher than z value. 
	 					pWriter.print("forward block is done i r y z"+" "+i+" "+r+" "+y+" "+z+" "+"matrix"+" "+Arrays.deepToString(A)+"\n");
	 					if(blockTest(i, r, A,pWriter)) {
	 						r++;
	 						forwardBuild(r, degrees, partitionList.get(z+1), A, max, L, C, modified, pWriter);
	 					}
	 				}else {
	 					forwardBuild(r, degrees, partitionList.get(y), A, max, L, C, modified, pWriter);
	 				}
	 			}
	 		}else {
	 			if(i==max.length-2 && j==max.length-1) {
	 				ArrayList<Integer> modified=predecessor(indices, max.length);
		 			indices=modified;
	 			}
	 			backwardBuild(r, degrees, partitionList.get(y),A, max, L, C, indices, pWriter);
	 		}
	 	}
	}**/
	
	/**
	 * In backward function, updating the block index.
	 * @param r int block index
	 * @param indices ArrayList<Integer> atom valences
	 * @return
	 */
	 
	public static int updateR(int r, ArrayList<Integer> indices) {
		int y=findY(r);
		if(indices.get(0)<y) {
			return (r-1);
		}else {
			return r;
		}
	}
	
	/**
	 * After generating matrices, adding the hydrogen with respect to
	 * the pre-hydrogen distribution. 
	 * @param A int[][] adjacency matrix
	 * @param index int beginning index for the hydrogen setting
	 * @return
	 */
	
	public static int[][] addHydrogens(int[][] A, int index){
		int hIndex=index;
		int hydrogen=0;
		for(int i=0;i<index;i++) {
			hydrogen=initialDegrees.get(i)-sum(A[i]);
			for(int j=hIndex;j<(hIndex+hydrogen);j++) {
				A[i][j]=1;
				A[j][i]=1;
			}
			if(hydrogen==0) {
				//hIndex++;
			}else {
				hIndex=hIndex+hydrogen;
			}
		}
		return A;
	}
	
	/**
	 * Summing entries of an array.
	 * @param array int[]
	 * @return int sum
	 */
	
	public static int sum(int[] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	/**
	 * Building an atom container for an adjacency matrix
	 * @param mat int[][] adjacency matrix
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static IAtomContainer buildC(int[][] mat) throws CloneNotSupportedException {
		IAtomContainer ac2= atomContainer.clone();
		for(int i=0;i<mat.length;i++) {
			for(int j=i+1;j<mat.length;j++) {
				if(mat[i][j]==1) {
					ac2.addBond(i, j, Order.SINGLE);;
				}else if(mat[i][j]==2) {
					ac2.addBond(i, j, Order.DOUBLE);
				}else if(mat[i][j]==3) {
					ac2.addBond(i, j, Order.TRIPLE);
				}
			}
		}
		
		ac2=AtomContainerManipulator.removeHydrogens(ac2);
		return ac2;
	}
	 
	/**
	 * Building an atom container from a string of atom-implicit hydrogen information.
	 * If provided, fragments are also added.
	 * @param mol molecular information
	 * @return atom container new atom container
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	
	public static void build(String mol) throws IOException, CloneNotSupportedException, CDKException {
		List<String> symbols = new ArrayList<String>();
	    List<Integer> number = new ArrayList<Integer>();
	    String[] atoms = mol.split("(?=[A-Z])");
	    for (String atom : atoms) {
	    	String[] info = atom.split("(?=[0-9])", 2);   
	        symbols.add(info[0]);
	        number.add(info.length > 1 ? Integer.parseInt(info[1]):1);
	    }
	    
	    for(int i=0;i<symbols.size();i++) {
	    	for(int j=0;j<number.get(i);j++) {
	    		atomContainer.addAtom(new Atom(symbols.get(i)));
	    	}
	    }
	    
	    for(IAtom atom: atomContainer.atoms()) {
	    	atom.setImplicitHydrogenCount(0);
	    }
	}
	
	/**
	 * Cloning an atom container
	 * @param mol IAtomContainer 
	 * @return
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	
	public static IAtomContainer clone(IAtomContainer mol) throws IOException, CloneNotSupportedException, CDKException {
		IAtomContainer mol2 = builder.newInstance(IAtomContainer.class);
		for(IAtom a: mol.atoms()) {
			mol2.addAtom(a);
		}
		for(IBond b: mol.bonds()) {
			mol2.addBond(b);
		}
		return mol2;
	}
	
	/**
	 * Molecule depiction 
	 */
		
	public static void depict(IAtomContainer molecule, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		depiction.withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
	}
		
	/**
	 * Backward function in the matrix generator.
	 * @param r int block index
	 * @param degrees ArrayList<Integer> valences
	 * @param partition ArrayList<Integer> atom partition
	 * @param A int[][] adjacency matrix
	 * @param max int[][] maximal matrix
	 * @param L int[][] L matrix
	 * @param C int[][] C matrix
	 * @param indices ArrayList<Integer> entry indices
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 * @throws CDKException	
	 */
	
	public static void backwardBuild(int r, ArrayList<Integer> degrees,ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException, CloneNotSupportedException, CDKException {
		r=updateR(r,indices);
		int y=findY(r);
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
		if(i==max.length-2 && j==max.length-1) {
			output.add(A);
			int[][] mat2= new int[A.length][A.length]; 
			for(int k=0;k<A.length;k++) {
				for(int l=0;l<A.length;l++) {
					mat2[k][l]=A[k][l];
				}
			}
			A=addHydrogens(A,hIndex);
			if(connectivityTest(hIndex,A)){
				IAtomContainer mol= buildC(A);
				String inchi= inchiGeneration(mol);
				if(!inchis.contains(inchi)) {
					outFile.write(clone(mol));
					count++;
					inchis.add(inchi);
				}
			}

			for(int h=0; h<representatives.size(); h++) {
				representatives.get(h).removeAll(representatives.get(h));
			}
			for(int k=(1); k<partitionList.size(); k++) {
				partitionList.get(k).removeAll(partitionList.get(k));
			}
			
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					ArrayList<Integer> modified2=successor(modified,max.length);
					forwardBuild(r, degrees, partitionList.get(y), A, max, L, C, modified2);
				}else {
					backwardBuild(r, degrees,partitionList.get(y), A, max, L, C, modified);
				}
			}
		}
	}
	
	/**
	 * Testing every row in the block with the canonical test. 
	 * @param index int row index
	 * @param r int the block number
	 * @param A int[][] adjacency matrix
	 * @param partition ArrayList<Integer> atom partition
	 * @param newPartition ArrayList<Integer> canonical partition
	 * @return boolean 
	 */
	 
	public static boolean blockTest(int index, int r, int[][] A) {
		boolean check=true;
		int y=findY(r);
		int z=findZ(r);
		 
		for(int i=y;i<=z;i++) {
			if(!block(i,r, A, partitionList.get(i),canonicalPartition(i,partitionList.get(i)))){
				check=false;
				break;
			}
		}
		clear(check, y);
		return check;
	}
	 
	/**
	 * If the blockTest is false, clear the reps and representatives 
	 * of the tested block. 
	 * @param check boolean blockTest
	 * @param y int beginning index of the block.
	 */
	 
	 public static void clear(boolean check, int y) {
		 if(check==false) {
			 for(int i=y;i<representatives.size();i++) {
				 representatives.get(i).removeAll(representatives.get(i));
			 }
			 for(int i=y+1;i<partitionList.size();i++) {
				 partitionList.get(i).removeAll(partitionList.get(i));
			 }
		 }
	 }
	 
	 /**
	  * Grund Thesis 3.3.3.
	  * To calculate the number of conjugacy classes, used in cycle 
	  * transposition calculation.
	  * @param partEx ArrayList<Integer> former atom partition
	  * @param degree ArrayList<Integer> atom valences
	  * @return
	  */
	 
	 public static int LValue(ArrayList<Integer> partEx, int degree) {
		 return (sum(partEx,(degree))-(degree)); //In Grund, the numeration starts with 1. Since we start with 0, it should be -(degree+1)+1, so just -degree
	 }
		
	 /**
	  * Grund Thesis 3.3.3.
	  * To generate cycle transpositions for two Youn subgroups.
	  * @param index int row index
	  * @param partition ArrayList<Integer> atom partition
	  * @return
	  */
	 
	 public static ArrayList<Permutation> cycleTranspositions(int index, ArrayList<Integer> partition) {
		 int total=sum(partition);
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 int lValue = LValue(partition,index);
		 for(int i=0;i<lValue;i++) {
	    	 int[] values= idValues(total);
	    	 int former =  values[index];
			 values[index] =  values[index+i];
			 values[index+i] = former;
			 Permutation p = new Permutation(values);
			 perms.add(p);
		 }
		 return perms;
	 }
	 
	 /**
	  * Within block test, testing whether the row is in maximal form
	  * or not.
	  * @param index int row index
	  * @param r int block index
	  * @param A int[][] adjacency matrices 
	  * @param partition ArrayList<Integer> atom partition
	  * @param newPartition ArrayList<Integer> new atom partition after canonical repartitioning
	  * @return
	  */
	 public static boolean block(int index, int r, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition) {
		 int y = findY(r);
		 boolean check= true;
		 int total = sum(partition);
		 ArrayList<Permutation> cycleTrans= cycleTranspositions(index,partition);
		 for(Permutation cycle: cycleTrans) {
			 boolean formerCheck = formerPermutationsCheck(index,y, total, A, partition, newPartition, cycle,false);
			 if(!formerCheck) {
				 check=false;
				 break;
			 }else {
				 addPartition(index, newPartition, A);
			 }
		 } 
		 return check;
	 }
	 
	 /**
	  * Updating canonical partition list. 
	  * @param index row index
	  * @param newPartition atom partition
	  * @param A int[][] adjacency matrix
	  */
	 
	 public static void addPartition(int index, ArrayList<Integer> newPartition, int[][] A) {
		 ArrayList<Integer> refinedPartition= refinedPartitioning(newPartition,A[index]);
		 if(partitionList.size()==(index+1)) {
			 partitionList.add(refinedPartition);
		 }else {
			 partitionList.set(index+1, refinedPartition);
		 }
	 }
	 
	 /**
	  * Grund 3.3.8
	  * In case if the row is canonical, with respect to former partition, 
	  * we build refined partitioning.
	  * @param partition ArrayList<Integer> atom partition
	  * @param row int[] row
	  * @return
	  */
	 
	 public static ArrayList<Integer> refinedPartitioning(ArrayList<Integer> partition, int[] row){
		 ArrayList<Integer> refined= new ArrayList<Integer>();
		 int index=0;
		 int count=1;
		 for(Integer p:partition) {
			 if(p!=1) {
				 for(int i=index;i<p+index-1;i++) {
					 if(i+1<p+index-1) { 
						 if(row[i]==row[i+1]) {
							 count++;
						 }else{
							 refined.add(count);
							 count=1;
						 } 
					 }else {
						 if(row[i]==row[i+1]) {
							 count++;
							 refined.add(count);
							 count=1;
						 }else{
							 refined.add(count);
							 refined.add(1);
							 count=1;
						 }
					 }
				 }
				 index=index+p;
			 }else {
				 index++;
				 refined.add(1);
				 count=1;
			 }
		 }
		 return refined;
	 }
	 
	 /**
	  * To calculate former permutations within the block or the permutations of former block.
	  * @param formerBlocks boolean
	  * @param index int row index
	  * @param y int y index of the block
	  * @param total int number of atoms
	  * @return
	  */
	 
	 public static ArrayList<Permutation> formerPermutations(boolean formerBlocks,int index, int y, int total){
		 if(formerBlocks) {
			 return formerPermutationsFromFormerBlocks(y);
		 }else {
			 return formerPermutationsInBlock(index,y,total);
		 }
	 }
	 
	 /**
	  * The multiplication of the former permutations from the former block.
	  * @param y int y index of the block.
	  * @return
	  */
	 
	 public static ArrayList<Permutation> formerPermutationsFromFormerBlocks(int y) {
		 ArrayList<Permutation> list= representatives.get(0);
		 ArrayList<Permutation> nList= new ArrayList<Permutation>();
		 for(int i=1;i<y;i++) {
			 for(int l=0;l<list.size();l++) {
				 for(Permutation perm: representatives.get(i)) {
					 Permutation mult= list.get(l).multiply(perm);
					 if(!nList.contains(mult)) {
						 nList.add(mult); 
					 }
				 }
			 }
			 list.clear();
			 list.addAll(nList); // list=nList does not work.
			 nList.clear();
		 }
		 return list; 
	 }
	 
	 /**
	  * Like in Grund Example 3.3.19; We need to check with the permutations 
	  * of the former representatives. This function helps us to calculate
	  * the multiplication of the former representatives. So, before
	  * the beginning of the current strip (or block).
	  * 
	  * These are M(i) members.
	  * 
	  * @param index int current row index
	  * @param y int beginning of the block
	  * @return
	  */
	 
	 public static ArrayList<Permutation> formerPermutationsInBlock(int index, int y, int total) {
		 if(index==y) {
			 ArrayList<Permutation> form= new ArrayList<Permutation>();
			 form.add(idPermutation(total));
			 return form;
		 }else {
			 ArrayList<Permutation> list= representatives.get(y);
			 ArrayList<Permutation> nList= new ArrayList<Permutation>();
			 for(int i=(y+1);i<index;i++) {
				 for(int l=0;l<list.size();l++) {
					 for(Permutation perm: representatives.get(i)) {
						 Permutation mult= list.get(l).multiply(perm);
						 //if(!mult.isIdentity()) {
						 if(!nList.contains(mult)) {
							 nList.add(mult); 
						 }
						 //}
					 }
				 }
				 list.clear();
				 list.addAll(nList); // list=nList does not work.
				 nList.clear();
			 }
			 return list; 
		 }
	 }
	 
	 /**
	  * To add a permutation to representatives.
	  * @param index representatives table index
	  * @param perm  Permutation
	  */
	 
	 public static void addRepresentatives(int index, Permutation perm) {
		 if(representatives.size()==index) {
			 ArrayList<Permutation> perms = new ArrayList<Permutation>();
			 perms.add(perm);
			 representatives.add(perms);
		 }else {
			 if(!representatives.get(index).contains(perm)) {
				 representatives.get(index).add(perm);
			 }
		 }
	 }
	 
	 /**
	  * To calculated the automorphism of a row, so it is still
	  * in maximal form and the permutation maps row to itself.
	  * @param index int row index
	  * @param y int y index of the block
	  * @param total int number of atoms
	  * @param A int[][] adjacency matrices
	  * @param partition ArrayList<Integer> atom partition
	  * @param newPartition ArrayList<Integer> canonical atom partition
	  * @param former Permutation former permutation
	  * @param cycle Permutation cycle coming from Grund 3.3.3.
	  * @return
	  */
	 
	 public static Permutation getCanonicalPermutatiom(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation former, Permutation cycle) {
		 Permutation cycleM= former.multiply(cycle);
		 Permutation canonical = idPermutation(total);
		 //if(descBlockwiseCheck2(index,y,A,newPartition)) {
		 if(descBlockwiseCheck2(index,y,A,newPartition)) {
			 if(!equalBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total))) {
				 canonical=getEqualPermutation(cycleM, index, y, total,A, partition, newPartition);
				 canonical=cycle.multiply(canonical);
				 /**if(canonical.isIdentity()) {
					 if(!descBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total),2,pWriter)) {
						 canonical = idPermutation(total);
					 }else {
						 canonical = cycle.multiply(canonical);
					 }
				 }else {
					 canonical=cycle.multiply(canonical);
				 }**/
			 }else {
				 canonical=cycle.multiply(canonical);
			 }
		 }
		 return canonical;
	 }
	 
	 /**
	  * This function performs the canonical test with the former representatives.
	  * That can be the reps from the former representatives in the block. 
	  * If it outputs the canonical permutation as id, then not
	  * passed the canonical test and remove already added new entries from the
	  * list of representatives.
	  * 
	  * @param index int row index
	  * @param y int the first index of the block
	  * @param total int number of atoms
	  * @param A int[][] adjacency matrix
	  * @param partition ArrayList<Integer> atom partition
	  * @param newPartition ArrayList<Integer> canonical partition
	  * @param perm Permutation perm
	  * @param formerBlocks boolean true if permutations are from former blocks
	  * @return Permutation
	  */
	 
	 public static boolean formerPermutationsCheck(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation perm, boolean formerBlocks) {
		 boolean check=true;
		 ArrayList<Permutation> formerPerms= formerPermutations( formerBlocks,index, y, total);
		 Permutation canonical = idPermutation(total);
		 for(Permutation former: formerPerms) {
			 Permutation test=former.multiply(perm);
			 //if(getCanonicalTest(index, y, total, A, partition, newPartition, test, pWriter)) {
			 canonical=getCanonicalPermutatiom(index, y, total, A, partition, newPartition, former,perm);
			 if(canonical.isIdentity()) {
				 //if(descBlockwiseCheck2(index,y,A,newPartition) && descBlockCheck(test,newPartition,index,y,A,idPermutation(total))) { //TODO: Not enought for the example. Not des but id so if the array is desc but the comparison is not. Then we cant say continue.
				 if(descBlockwiseCheck2(index,y,A,newPartition)) {	 
					 addRepresentatives(index,idPermutation(total));
					 /**boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
						check=false;
						break;
					 }else {
						addRepresentatives(index,idPermutation(total), pWriter); 
					 }**/
				 }else {
					 check=false;
					 break;
				 }
			 }else {
				 if(test.equals(former.multiply(canonical))) {
					 addRepresentatives(index, idPermutation(total)); 
					 /**boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
					 	check=false;
						break;
				 	 }else {
						addRepresentatives(index, idPermutation(total), pWriter); 
				  	 }**/
				 }else {
					 addRepresentatives(index, canonical);
					 /**boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
						check=false;
						break;
					 }else {
						pWriter.print("add rep"+" "+index+" "+representatives.size()+" "+canonical.toCycleString()+" "+"\n");
						addRepresentatives(index, canonical,pWriter);
					 }**/
				 }
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Former block representatives: After the test within the block,
	  * we test the row with the former blocks representatives. 
	  * @param index int row index
	  * @param y int first row index in block
	  * @param A int[][] adjacency matrix
	  * @param newPartition ArrayList<Integer> atom partition
	  * @param perm Permutation
	  * @return boolean
	  */
	 
	 public static boolean formerBlocksRepresentatives(int index, int y, int[][] A, ArrayList<Integer> newPartition, Permutation perm) {
		 boolean check=true;
		 ArrayList<ArrayList<Permutation>> removeList = describeRemovalList(y);
		 /**
		  * If we break within  the nested loop, need to put a label to 
		  * the outer loop also to break.
		  */
		 outer: for(int i=(y-1);i>=0;i--) {	
			 
			 /**
			  * No need to multiply them all since we test line by line
			  * in the list of representatives.
			  */
			 
			 for(Permutation rep: representatives.get(i)) {
				 Permutation test= rep.multiply(perm);
				 if(!descBlockCheck(test,newPartition,index,y,A,test)) {
					 check=false;
					 break outer;
				 }else {
					 if(!equalBlockCheck(test,newPartition,index,y,A,test)) {
						 if(!rep.isIdentity()) {
							 removeList.get(i).add(rep); 
						 }
					 }
				 }
			 }
		 }
		 if(check) {
			 removeFormerPermutations(y,removeList);
		 }
		 return check; 
	 }
	 
	 /**
	  * Grund Thesis 3.3.19.
	  * If the permutations of the former block keep the row in descending order
	  * or equal, we can avoid that permutations and remove from the list.
	  * @param y int y index of the block
	  * @param removeList ArrayList<ArrayList<Permutation>> Permutations to remove from the list
	  */
	 
	 public static void removeFormerPermutations(int y,ArrayList<ArrayList<Permutation>> removeList) {
		 ArrayList<Permutation> check= new ArrayList<Permutation>();
		 for(int i=(y-1);i>=0;i--) {
			 check=removeList.get(i);
			 if(check.size()!=0) {
				 representatives.get(i).removeAll(check);
			 }
		 }
	 }
	 
	 /**
	  * Check whether the modified and original rows are equal or not.
	  * @param partition atom partition
	  * @param index row index
	  * @param y beginning index of a block
	  * @param A adjacency matrix
	  * @param perm permutation
	  * @return
	  */
	 
	 public static boolean equalBlockCheck(Permutation cycleM,ArrayList<Integer> partition, int index, int y, int[][] A, Permutation perm){
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] original= A[index];
		 Permutation mult=cycleM.multiply(perm);
		 original=actArray(A[mult.get(index)],mult);
		 if(!Arrays.equals(canonical,original)) {
			 check=false;
		 }
		 return check;
	 }
	 
	 /**
	  * Describing ArrayList<ArrayList<Permutation>>
	  * for formerBlocksRepresentatives function.
	  * 
	  * @param y int index of the first row in block
	  * @return ArrayList<ArrayList<Permutation>> permutations to remove.
	  */
	 
	 public static ArrayList<ArrayList<Permutation>> describeRemovalList(int y){
		 ArrayList<ArrayList<Permutation>> removeList = new ArrayList<ArrayList<Permutation>>();
		 for(int i=0;i<y;i++) {
			 ArrayList<Permutation> list = new ArrayList<Permutation>();
			 removeList.add(i,list);
		 }
		 return removeList;
	 }	
		
	 public static Integer [] getBlocks(int[] row, int begin, int end) {
		 Integer[] sub=subArray(row,begin,end);
		 return sub;
	 }
		
	 public static Integer[] subArray(int[] array, int begin,int end) {
		 return IntStream.range(begin, end).mapToObj(i->array[i]).toArray(Integer[]::new);
	 }
		
	 /**
	  * Descending order check for a integer array
	  * @param array int[] 
	  * @return boolean
	  */
		
	 public static boolean desOrderCheck(Integer[] array) {
		 boolean check=true;
		 int length=array.length;
		 for(int i=0;i<length-1;i++) {
			 if(array[i]<array[i+1]) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }	
		
	 /**
	  * int array to int representation. This is needed for blockwise descending order check. 
	  * @param array
	  * @return
	  */
			
	 public static int toInt(Integer[] array) {
		 int result=0;
		 for(int i=0;i<array.length;i++) {
			 result=result*10;
			 result=result+array[i];
		 }
		 return result;
	 }
	 
	 /**
	  * To check whether the block is in descending order after the permutation actions.
	  * @param cycle Permutation cycle
	  * @param partition ArrayList<Integer> atom partition
	  * @param index int row index
	  * @param y int y index of the block
	  * @param A int[][] adjacency matrix
	  * @param perm Permutation
	  * @return
	  */
	 public static boolean descBlockCheck(Permutation cycle,ArrayList<Integer> partition, int index, int y, int[][] A, Permutation perm){
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] original= A[index];
		 original=actArray(A[cycle.get(index)],cycle);
		 int[] temp= new int[sum(partition)];
		 if(cycle.get(index)<index) {
			 temp=canonical;
			 canonical=original;
			 original=temp;
		 }

		 //int[] original=actArray(A[index],perm);
		 
		 int i=0;
		 for(Integer p:partition) {
			 Integer[] can= getBlocks(canonical,i,p+i);
			 Integer[] org= getBlocks(original,i,p+i);
			 //if(desOrderCheck(can) && desOrderCheck(org)) { not both the first should be desc
			 if(desOrderCheck(can)) {	 
			 	if(!Arrays.equals(can,org)) {
					 if(toInt(can)==toInt(org)) {
						 continue;
					 }else if(toInt(can)>toInt(org)) {
						 check=true;
						 break; //If bigger than already in lexico order.
					 }else if(toInt(can)<toInt(org)) {
						 check=false;
						 break;
					 }
				 }
				 /**else {
					 if(!descendingOrderCheck(can)) {
						 check=false;
						 break;
					 }
				 }**/
				 i=i+p;
			 }else {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Within blockwise checking the row is descending or not after the permutation action.
	  * @param index int row index
	  * @param y int y index of the block
	  * @param A int[][] adjacency matrix
	  * @param partition ArrayList<Integer> partition
	  * @return
	  */
	 
	 public static boolean descBlockwiseCheck2(int index, int y, int[][] A, ArrayList<Integer> partition){
		 int[] array= A[index];
		 boolean check=true;
		 int i=0;
		 for(Integer p:partition) {
			 Integer[] can= getBlocks(array,i,p+i);
			 if(desOrderCheck(can)) {
				 i=i+p;
				 continue;
			 }else {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**public static Permutation getCanonicalPermutation(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation cycleM, PrintWriter pWriter) {
		 pWriter.print("getCanonicalPermutation"+"\n");
		 pWriter.print("permutation"+" "+cycleM+"\n");
		 //PermutationGroup group= getYoungGroup(newPartition,total);
		 Permutation canonical = idPermutation(total);
		 pWriter.print("des"+" "+descBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total))+" "+"\n");
		 pWriter.print("equal"+" "+equalBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total), pWriter)+" "+"\n");
		 if(descBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total))){
			 if(!equalBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total), pWriter)) {
				 canonical=getEqualPermutation(cycleM, index, y, total,A, partition, newPartition, pWriter);
				 pWriter.print("equal perm"+" "+canonical.toCycleString()+"\n");
				 canonical=cycleM.multiply(canonical);
			 }else {
				 /**
				  * If id then not canonical but if cycle then 
				  * no need to add anything to reps and ignore
				  * that cycle for the further steps.
				  */
				 //canonical=cycleM.multiply(canonical);
			 //}
		 //}else {
			 
			 /**
			  * des already check equal or descending.
			  */
			 
			 /**canonical=getEqualPermutation(cycleM, index, y, total, A, partition, newPartition, pWriter);
			 if(!canonical.isIdentity()) {
				 canonical=cycleM.multiply(canonical);
			 }**/
		 //}
		 //return canonical;
	 //}
	 
	 /**
	  * Checks whether there is an automorphism permutation mapping
	  * the row itself and keeping it in maximal form.
	  * @param cycle Permutation cycle
	  * @param index int row index
	  * @param y int y index of the block
	  * @param total int number of atoms
	  * @param A int[][] adjacency matrix
	  * @param partition ArrayList<Integer> atom partition
	  * @param newPartition ArrayList<Integer> new partition after canonical re-partitioning
	  * @return
	  */
	 
	 public static Permutation getEqualPermutation(Permutation cycle, int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition) {
		 PermutationGroup group= getYoungGroup(newPartition,total);
		 Permutation canonical = idPermutation(total);
		 for(Permutation perm : group.all()) {
			 if(!perm.isIdentity()) {
				 if(equalBlockCheck(cycle,newPartition,index,y,A,perm)){
					 canonical = perm;
					 break;
				 }
			 }
		 }
		 return canonical;
	 }
		
	 public static ArrayList<ArrayList<Integer>> subPartitionsList(ArrayList<Integer> partition){
		 ArrayList<ArrayList<Integer>> parts = new ArrayList<ArrayList<Integer>>();
		 for(int i=0;i<partition.size();i++) {
			 parts.add(subPartition(partition,i)); 
		 }
		 return parts;
	 }
	
	 /**
	  * To get the canonical partition like in Grund Thesis 3.3.11
	  * @param i int row index
	  * @param partition ArrayList<Integer> partition
	  * @return
	  */
	 
	 public static ArrayList<Integer> canonicalPartition(int i, ArrayList<Integer> partition){
		 if(i==0) {
			 return partitionCriteria(partition,1);
		 }else {
			 return partitionCriteria(partition,i+1);
		 }
	 }
	 
	 /**
	  *	Grund Thesis 3.3.2 Partitioning criteria (DONE)  
	  * @param partEx the former partition
	  * @param degree degree of the partitioning.
	  * @return
	  */
		 
	 public static ArrayList<Integer> partitionCriteria(ArrayList<Integer> partEx, int degree){
		 ArrayList<Integer> partNew = new ArrayList<Integer>();
		 addOnes(partNew,degree);
		 if(partEx.get(degree-1)>1) {
			 partNew.add(partEx.get(degree-1)-1);
			 for(int k=degree;k<partEx.size();k++) {
				 partNew.add(partEx.get(k));
			 }
		 }else if(partEx.get(degree-1)==1){
			 for(int k=degree;k<partEx.size();k++) {
				 partNew.add(partEx.get(k));
			 }
		 }
		 return partNew;
	 }
		 
	 /**
	  * add number of 1s into a ArrayList
	  */
		 
	 public static void addOnes(ArrayList<Integer> list, int number) {
		 for(int i=0;i<number;i++) {
			 list.add(1);
		 }
	 }
	
	 
	 /**
	  * 3.6.2. Connectivity Test
	  */
		
	 /**
	  * Finding the N values as explained in 3.6.2.
	  * @param index int row index
	  * @param total int number of atoms.
	  * @param mat int[][] adjacency matrix
	  * @return Set<Integer>
	  */
		 
	 public static Set<Integer> nValues(int index, int total, int[][] mat) {
		 Set<Integer> nValues= new HashSet<Integer>();
		 nValues.add(index);
		 int[] theRow = mat[index]; 
		 for(int i=(index+1);i<total;i++) {
			 if(theRow[i]>0) {
				 nValues.add(i);
			 }
		 }
		 return nValues;
	 }
		 
	 /**
	  * Finding the W values as explained in 3.6.2.
	  * @param nValues ArrayList<Integer> N values 
	  * @param Kformer ArrayList<Integer> the K values of the former step
	  * @return Set<Integer>
	  */
		 
	 public static Set<Integer> wValues(Set<Integer> nValues, ArrayList<Integer> Kformer){
		 Set<Integer> wValues= new HashSet<Integer>();
		 for(Integer i:nValues) {
			 wValues.add(Kformer.get(i));
		 }
		 return wValues;
	 }
		 
	 /**
	  * Finding the K values as explained in 3.6.2.
	  * @param wValues Set<Integer> wValues  
	  * @param kFormer ArrayList<Integer> the K values of the former step
	  * @return ArrayList<Integer>
	  */
		 
	 public static ArrayList<Integer> kValues(int total, Set<Integer> wValues, ArrayList<Integer> kFormer){
		 ArrayList<Integer> kValues= new ArrayList<Integer>();
		 int min= Collections.min(wValues);
		 for(int i=0;i<total;i++) {
			 if(wValues.contains(kFormer.get(i))) {
				 kValues.add(i,min);
			 }else {
				 kValues.add(i,kFormer.get(i));
			 }
		 }
		 return kValues;
	 }
		 
	 /**
	  * Test whether an adjacency matrix is connected or disconnected.
	  * @param p int the index where the atoms with degree 1 starts.
	  * @param mat int[][] adjacency matrix
	  * @return boolean
	  */
		 
	 public static boolean connectivityTest(int p, int[][] mat) {
		 boolean check=false;
		 int total= mat.length;
		 ArrayList<Integer> kValues=initialKList(total);
		 Set<Integer> nValues= new HashSet<Integer>();
		 Set<Integer> wValues= new HashSet<Integer>();
		 int zValue= 0;
		 for(int i=0;i<p;i++) {
			 nValues= nValues(i, total, mat);
			 wValues=wValues(nValues,kValues);
			 zValue = Collections.min(wValues);
			 kValues= kValues(total, wValues, kValues);
		 }
		 if(zValue==0 && allIs0(kValues)) {
			 check=true;
		 }
		 return check;
	 }
		
		 
	 /**
	  * Checks whether all the entries are equal to 0 or not. (Grund 3.6.2)
	  * @param list ArrayList<Integer>
	  * @return boolean
	  */
		 
	 public static boolean allIs0(ArrayList<Integer> list) {
		 boolean check=true;
		 for(int i=0;i<list.size();i++) {
			 if(list.get(i)!=0) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
		 
	 /**
	  * Initializing the first K values list.
	  * @param total int number of atoms.
	  * @return ArrayList<Integer>
	  */
		 
	 public static ArrayList<Integer> initialKList(int total){
		 ArrayList<Integer> k= new ArrayList<Integer>();
		 for(int i=0;i<total;i++) {
			 k.add(i);
		 }
		 return k;
	 }
		
	 /**
	  * The initializer function, reading the formula to set the degrees,
	  * partition and file directory variables.
	  * @param formula String molecular formula
	  * @param filedir String file directory
	  * @throws IOException
	  * @throws CDKException
	  * @throws CloneNotSupportedException
	  */
	 
	 public static void run(String formula, String filedir) throws IOException, CDKException, CloneNotSupportedException {
		 if(verbose) System.out.println("For the formula, "+formula+",start generating all possible connectivity matrices...");
		 build(formula);
		 outFile = new SDFWriter(new FileWriter(filedir+"output.sdf"));
		 canonicalBlockbasedGenerator(degreeSeq(formula),distribution(formula),filedir);
		 if(verbose) System.out.println("The number of structures is: "+count);	
		 outFile.close();
	 }
		
	 /**
	  * Degree sequence is set from formula
	  * @param formula String
	  * @return
	  */
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
	  * Degree sequence is set from formula
	  * @param formula String
	  * @return
	  */
	 
	 public static ArrayList<Integer> distribution(String formula){
		 List<String> symbols = new ArrayList<String>();
		 ArrayList<Integer> occur  = new ArrayList<Integer>();
		 String[] atoms = formula.split("(?=[A-Z])");
		 for (String atom : atoms) {
			 String[] info = atom.split("(?=[0-9])", 2);   
		     symbols.add(info[0]);
		     occur.add(info.length > 1 ? Integer.parseInt(info[1]):1);
		 }
		 return occur;
	 } 
	
	 private void parseArgs(String[] args) throws ParseException, IOException{
		 Options options = setupOptions(args);	
		 CommandLineParser parser = new DefaultParser();
		 try {
			 CommandLine cmd = parser.parse(options, args);
			 buNe.formula = cmd.getOptionValue("formula");
			 buNe.filedir = cmd.getOptionValue("filedir");
			 //setFileWriter();
				 
			 if (cmd.hasOption("verbose")) buNe.verbose = true;		
		 } catch (ParseException e) {
			 HelpFormatter formatter = new HelpFormatter();
			 formatter.setOptionComparator(null);
			 String header = "\nGenerates 	molecular structures for a given molecular formula."
						 + " The input is a molecular formula string."
						 + "For example 'C2OH4'."
						 +"For this version, the hydrogens should be kept at the end of the string."
						 + "Besides this formula, the directory is needed to be specified for the output"
						 + "file. \n\n";
			 String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory";
			 formatter.printHelp( "java -jar AlgorithmicGroupTheory.jar", header, options, footer, true );
			 throw new ParseException("Problem parsing command line");
		 }
	 }
				
	 private Options setupOptions(String[] args){
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
		 ArrayList<Integer> degrees= new ArrayList<Integer>();
		 degrees.add(3);
		 degrees.add(3);
		 degrees.add(3);
		 degrees.add(3);
		 degrees.add(3);
		 degrees.add(3);
		
			
		 ArrayList<Integer> partition= new ArrayList<Integer>();
		 partition.add(6);
		 
		 buNe.formula="C6H6";
		 buNe.filedir="C:\\Users\\mehme\\Desktop\\";
		 build(formula);
		 outFile = new SDFWriter(new FileWriter(filedir+"output.sdf"));
		 canonicalBlockbasedGenerator(degrees,partition,filedir);
		 outFile.close();
		 
		 /**buNe gen = null;
		 //String[] args1= {"-f","C4OH10","-v","-d","C:\\Users\\mehme\\Desktop\\"};
		 try {
			 gen = new buNe();
			 gen.parseArgs(args);
			 buNe.run(buNe.formula,buNe.filedir);
		 } catch (Exception e) {
			 if (buNe.verbose) e.getCause(); 
		 }**/
	 }
}
