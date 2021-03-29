package AlgorithmicGroupTheory;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class MORGEN {
	public static int size=0;
	public static int[] lernenIndices = new int[2];
	public static int lernenC=0;
	public static int hIndex=0;
	public static int count=0;
	public static int matrixSize=0;
	public static boolean verbose = false;
	public static SDFWriter outFile;
	public static String formula;
	public static String filedir;
	public static boolean flag=true;
	public static boolean lernen=false;
	public static boolean biggest=true;
	public static List<List<Permutation>> formerPermutations= new ArrayList<List<Permutation>>();
	public static int[] degrees;
	public static int[] initialDegrees;
	public static ArrayList<Integer> initialPartition;
	public static IChemObjectBuilder builder=DefaultChemObjectBuilder.getInstance();
	public static IAtomContainer atomContainer= builder.newInstance(IAtomContainer.class);
	public static List<ArrayList<Integer>> partitionList= new ArrayList<ArrayList<Integer>>();
	public static List<String> symbols = new ArrayList<String>();
	public static ArrayList<Integer> occurrences  = new ArrayList<Integer>();
	public static Map<String, Integer> valences; 
	public static PrintWriter pWriter;
	public static int[][] max;
	public static int[][] L;
	public static int[][] C;
	public static int r=0; 
	public static int y=0;
	public static int z=0;
	static {
		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
			
		valences.put("C", 4);
		valences.put("N", 3);
		valences.put("O", 2);
		valences.put("S", 2);
		valences.put("P", 3);
		valences.put("F", 1);
		valences.put("I", 1);
		valences.put("Cl", 1);
		valences.put("Br", 1);
		valences.put("H", 1);
	}
	
	/**
	 * Basic functions
	 */
	
	/**
	 * Molecule depiction 
	 */
		
	public static void depict(IAtomContainer molecule, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		depiction.withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
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
	 
	 public static Integer[] permuteArray(Integer[] array, int i, int j) {
    	 int temp=0;
    	 temp= array[i];
    	 array[i]=array[j];
    	 array[j]=temp;
    	 return array; 
     }
	 
	/**
	 * Summing entries of a list.
	 * @param list List<Integer>
	 * @return int sum
	 */
	
	public static int sum(List<Integer> list) {
		int sum=0;
		for(int i=0;i<list.size();i++) {
			sum=sum+list.get(i);
		}
		return sum;
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
	
	public static int sum(ArrayList<Integer> list, int index) {
		int sum=0;
		for(int i=0;i<=index;i++) {
			sum=sum+list.get(i);
		}
		return sum;
	}
	
	public static int atomOccurrunce(String[] info) {
		 int number=1;
		 if(info.length >1) {
			 number=Integer.parseInt(info[1]);
		 }
		 return number;
	}
	
	public static int[] actArray(int[] strip, Permutation p) {
		int permLength= p.size();
		int arrayLength=strip.length;
		int[] modified = new int[arrayLength];
		for(int i=0; i<permLength;i++) {
			modified[p.get(i)]=strip[i]; 
		}
		return modified;
	}
	
	/**
	  * Values for an id permutation for a given size
	  * @param size 
	  * @return
	  */
	 public static int[] idValues(int size) {
		 int[] id= new int[size];
		 for(int i=0;i<size;i++) {
			 id[i]=i;
		 }
		 return id;
	 }
	 
	 /**
	  * Build id permutation
	  */
	 
	 public static Permutation idPermutation(int size) {
		 Permutation perm= new Permutation(size);
		 return perm;
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
	 
	
	 /**
	  * Setting the symbols and occurrences global variables.
	  * @param formula  String molecular formula
	  */
	 
	 public static List<String> firstSymbols= new ArrayList<String>();
	 public static ArrayList<Integer> firstOccurrences = new ArrayList<Integer>();
	 public static void getSymbolsOccurrences(String formula) {
		 String[] atoms = formula.split("(?=[A-Z])");
		 for(String atom : atoms) {
			 String[] info = atom.split("(?=[0-9])", 2); 
			 matrixSize=matrixSize+atomOccurrunce(info);
			 firstSymbols.add(info[0]);
			 firstOccurrences.add(atomOccurrunce(info));
			 if(!info[0].equals("H")) {
				 symbols.add(info[0]);
				 occurrences.add(atomOccurrunce(info));
				 hIndex=hIndex+atomOccurrunce(info);
			 }
		 }
	 }
	
	 /**
	  * Degree sequence is set from formula
	  * @param formula String molecular formula
	  * @return int[] valences 
	  */
	 public static int[] firstDegrees;
	 public static void initialDegrees(){
		 int size= sum(occurrences);
		 initialDegrees= new int[size];
		 firstDegrees= new int[sum(firstOccurrences)];
		 int index=0;
		 int firstIndex=0;
		 for(int i=0;i<symbols.size();i++) {
			 String symbol= symbols.get(i);
		     for(int j=0;j<occurrences.get(i);j++) {
		    	 if(!symbol.equals("H")) {
		    		 initialDegrees[index]=valences.get(symbol);
			    	 index++; 
		    	 }
		     }
		 }
		 for(int i=0; i<firstSymbols.size();i++) {
			 String symbol= firstSymbols.get(i);
			 for(int j=0;j<firstOccurrences.get(i);j++) {
		    	 firstDegrees[firstIndex]=valences.get(symbol);
		    	 firstIndex++; 
		     }
		 }
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
		 for(int i=0;i<firstSymbols.size();i++) {
			 for(int j=0;j<firstOccurrences.get(i);j++) {
				 atomContainer.addAtom(new Atom(firstSymbols.get(i)));
			 }
		 }
		    
		 for(IAtom atom: atomContainer.atoms()) {
			 atom.setImplicitHydrogenCount(0);
		 }
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
	  * ********************************************************************
	  */
	 
	 /**
	  * Sorting and Comparison Functions
	  */
	
	 public static boolean equalBlockCheck(ArrayList<Integer> partition, int index, int y, int[][] A, Permutation cycleTransposition, Permutation former, Permutation perm){
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] array=row2compare(index, A, cycleTransposition, former);
		 if(!Arrays.equals(canonical,array)) {
			 check=false;
		 }
		 return check;
	 }
	 
	 public static boolean equalSetCheck(int[] original, int[] permuted, ArrayList<Integer> partition) {
		 int[] temp= cloneArray(permuted);
		 temp=descendingSortWithPartition(temp, partition);
		 return equalSetCheck(partition, original, temp);
	 }
	 
	 public static Integer [] getBlocks(int[] row, int begin, int end) {
		 return IntStream.range(begin, end).mapToObj(i->row[i]).toArray(Integer[]::new);
	 }
	 
	 public static boolean equalSetCheck(ArrayList<Integer> partition, int[] original, int[] permuted){
		 boolean check=true;		 
		 int i=0;
		 for(Integer p:partition) {
			 Integer[] org= getBlocks(original,i,p+i);
			 Integer[] perm= getBlocks(permuted,i,p+i);
			 if(Arrays.equals(org, perm)) {	 
				 i=i+p;
			 }else {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static boolean equalBlockCheck(ArrayList<Integer> partition, int index, int y, int[][] A, Permutation cycleTransposition,Permutation perm){
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] original= A[index];
		 original=actArray(actArray(A[cycleTransposition.get(index)],cycleTransposition),perm); 
		 if(!Arrays.equals(canonical,original)) {
			 check=false;
		 }
		 return check;
	 }
	 
	 public static boolean equalCheck(int[] former, int[] current, ArrayList<Integer> partition) {
		 boolean check=true;
		 int i=0;
		 for(int k=0;k<partition.size();k++) {
			 Integer[] can= getBlocks(former,i,partition.get(k)+i);
			 Integer[] org= getBlocks(current,i,partition.get(k)+i); 
			 if(!Arrays.equals(can,org)) {
				 check=false;
				 break;
			 }else {
				 i=i+partition.get(k);
				 continue;
			 }
		 }
		 return check;
	 }
	 
	 public static int[] descendingSort(int[] array, int index0, int index1) {
		 int temp=0;
		 for (int i = index0; i < index1; i++) {     
			 for (int j = i+1; j < index1; j++) {     
				 if(array[i] < array[j]) {    
					 temp = array[i];    
	                 array[i] = array[j];    
	                 array[j] = temp;    
	             }     
	         }     
		 }
		 return array;
	 }
	 
	 public static int[] descendingSortWithPartition(int[] array, ArrayList<Integer> partition) {
		 int i=0;
		 for(Integer p:partition) {
			 array=descendingSort(array,i,i+p);
			 i=i+p;
		 }
		 return array;
	 }
	 
	 public static boolean biggerCheck(int index, int[] original, int[] permuted, ArrayList<Integer> partition) {
		 int[] sorted= cloneArray(permuted);
		 sorted=descendingSortWithPartition(sorted, partition);
		 return descendingOrderUpperMatrix(index,partition, original, sorted);
	 }
	 
	 public static void biggerCheck(int index, int[][] A, Permutation permutation, ArrayList<Integer> partition) {
		 biggest=true;
		 int[] check = row2compare(index, A, permutation);
		 if(!biggerCheck(index, A[index],check,partition)) {
			 biggest=false;
		 }
	 }
	 
	 public static boolean firstCheck(int index, int[][]A, ArrayList<Integer> partition) {
		 boolean check = true;
		 if(!descendingOrdercomparison(index, partition,A[index])) {
			 check=false;
			 int[] array = cloneArray(A[index]);
			 array = descendingSortWithPartition(array, partition);
			 Permutation canonicalPermutation = getCanonicalPermutation(array,A[index],partition);
			 lernen=true;  
			 lernenIndices = upperIndex(index, index, A, canonicalPermutation); 
		 }
		 return check;
	 }
	 
	 /**
	  * Learning from the canonical test - Grund  Chapter 3.4.1
	  */
		 
	 /**
	  * By a given permutation, checking which entry is mapped to the index.
	  * @param perm Permutation
	  * @param index int entry index in the row
	  * @return int
	  */
	 
	 public static int getPermutedIndex(Permutation perm, int index) {
		 int out=0;
		 for(int i=0; i<perm.size();i++) {
			 if(perm.get(i)==index) {
				 out+=i;
				 break;
			 }
		 }
		 return out;
	 }
	 /**
	  * Looking for the upper limit where the original entry
	  * is smaller. 
	  * @param index int row index
	  * @param size int row size
	  * @param original int[] original row to compare
	  * @param perm Permutation permutation from canonical test
	  * @return int upper limit.
	  */
		
	 
	 public static int[] limit(int index, int nextRowIndex, int[][] A, Permutation permutation) {
		 int[] original= A[index];
		 int[] permuted= A[nextRowIndex];
		 int[] limit= new int[2];
		 limit[0]=index;
		 int newIndex, value, newValue;
		 for(int i=index+1; i<size;i++) {
			 newIndex= getPermutedIndex(permutation, i);
			 value= original[i];
			 newValue= permuted[newIndex];
			 if(value!= newValue) {
				 if(value< newValue) {
					 limit[1]=i;
				 }
				 break;
			 }
		 }	
		 return limit;
		}
		
		/**
		 * Looking for the maximum index where the entry is not zero.
		 * @param index int row index
		 * @param size int row size
		 * @param original int[] original row to compare
		 * @param perm Permutation permutation from canonical test
		 * @return int lowerIndex
		 */
		
		public static int[] lowerIndex(int index, int nextRowIndex, int[][] A, Permutation permutation) {
			int max=0;
			int upperLimit=limit(index, nextRowIndex, A, permutation)[1];
			int[] permuted= A[nextRowIndex];
			int newIndex, newValue;
			for(int i=index+1; i<upperLimit;i++) {
				newIndex = getPermutedIndex(permutation, i);
				newValue = permuted[newIndex];
				if(newValue>0) {
					if(max<newIndex) {
						max=newIndex;
					}
				}
			}
			int[] pair= new int[2];
			pair[0]=nextRowIndex;
			pair[1]=max;
			return pair;
		}
		
		/**
		 * As explained in Grund 3.4.1; we need to calculate upperIndex.
		 * First, we need our j index where the original row become 
		 * smaller, then calculating the lower index. Based on these two
		 * values and the value of j in the permutation, we calculate
		 * our upper index. 
		 * 
		 * This upper index is used for the 'learning from canonical test'
		 * method. 
		 * 
		 * @param index int row index
		 * @param size  int row size
		 * @param original int[] original row to compare
		 * @param perm Permutation permutation from canonical test
		 * @return int upper index
		 */
		
		public static int[] upperIndex(int index, int nextRowIndex, int[][] A, Permutation permutation) {
			int[] limit  = limit(index, nextRowIndex, A, permutation);
			int[] lowerLimit = lowerIndex(index, nextRowIndex, A, permutation);
			int[] upperLimit= new int[2];
			upperLimit[0]=nextRowIndex;
			upperLimit[1]=getPermutedIndex(permutation, limit[1]);
			int[] maximalIndices = getMaximumPair(upperLimit,getMaximumPair(limit,lowerLimit));			
			maximalIndices=maximalIndexWithNonZeroEntry(A, maximalIndices);		
			return getTranspose(maximalIndices, A);
		}
		
		public static boolean ascending(int[] indices) {
			boolean check=true;
			if(indices[0]>=indices[1]) {
				check=false;
			}
			return check;
		}
		
		public static int[] getOne(int[] indices, int[][] A) {
			int[] output= new int[2];
			output[0]=indices[0];
			output[1]=indices[0]+1;
			return output;
		}
		
		public static int[] maximalIndexWithNonZeroEntry(int[][] A, int[] maximalIndices) {
			int rowIndex= maximalIndices[0];
			int columnIndex= maximalIndices[1];
			if((columnIndex>rowIndex) && A[rowIndex][columnIndex]!=0) {
				return maximalIndices;
			}else {
				int[] output= new int[2];
				for(int i=columnIndex; i<size;i++) {
					if(A[rowIndex][i]>0) {
						output[0]=rowIndex;
						output[1]=i;
						break;
					}
				}
				return output;
			}
		}
		
		public static int[] getBetween(int[][] A, int[] maximalIndices) {
			int rowIndex= maximalIndices[0];
			int columnIndex= maximalIndices[1];
			int[] output= new int[2];
			for(int i=columnIndex+1; i<size;i++) {
				if(A[rowIndex][i]>0) {
					output[0]=rowIndex;
					output[1]=i;
					break;
				}
			}
			return output;
		}
		
		public static int[] getTranspose(int[] indices, int[][] A) {
			int[] out= new int[2];
			if(indices[0]>indices[1]) {
				out[0]=indices[1];
				out[1]=indices[0];
				return out;
			}else {
				return indices;
			}
		}
		
		public static int[] getMaximumPair(int[] a, int[] b) {
			if(a[0]>b[0]) {
				return a;
			}else if(b[0]>a[0]) {
				return b;
			}else {
				if(a[1]>b[1]) {
					return a;
				}else if(b[1]>a[1]) {
					return b;
				}else {
					return a;
				}
			}
		}
	
	 public static boolean descendingOrderUpperMatrix(int index, ArrayList<Integer> partition, int[] original, int[] permuted){
		 boolean check=true;		 
		 int i=index+1;
		 for(int k=index+1;k<partition.size();k++) {
			 int p = partition.get(k);
			 Integer[] org= getBlocks(original,i,p+i);
			 Integer[] perm= getBlocks(permuted,i,p+i);
			 if(desOrderCheck(org)) {	 
			 	if(!Arrays.equals(org,perm)) {
					 if(toInt(org)==toInt(perm)) {
						 continue;
					 }else if(toInt(org)>toInt(perm)) {
						 check=true;
						 break; //If bigger than already in lexico order.
					 }else if(toInt(org)<toInt(perm)) {
						 check=false;
						 break;
					 }
				 }
				 i=i+p;
			 }else {
				 check=false;
				 break;
			 }
		 }
		 return check;
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
	 	 
	 public static boolean descendingOrdercomparison(ArrayList<Integer> partition, int[] original, int[] permuted){
		 boolean check=true;		 
		 int i=0;
		 for(Integer p:partition) {
			 Integer[] org = getBlocks(original,i,p+i);
			 Integer[] perm = getBlocks(permuted,i,p+i);
			 if(desOrderCheck(org)) {	 
			 	if(!Arrays.equals(org,perm)) {
					 if(toInt(org)==toInt(perm)) {
						 continue;
					 }else if(toInt(org)>toInt(perm)) {
						 check=true;
						 break; //If bigger than already in lexico order.
					 }else if(toInt(org)<toInt(perm)) {
						 check=false;
						 break;
					 }
				 }
				 i=i+p;
			 }else {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static boolean descendingOrdercomparison(int index, ArrayList<Integer> partition, int[] row){
		 boolean check=true;		 
		 int i=0;
		 for(Integer p:partition) {
			 Integer[] org= getBlocks(row,i,p+i);
			 if(!desOrderCheck(org)) {
				 check=false;
				 break;
			 }else {
				 i=i+p;
			 }
		 }
		 return check;
	 }
	 
	/**
	 *********************************************************************
	 */
	
	/**
	 * Candidate Matrix Generation Functions
	 */
	
	 /**
	  * L; upper triangular matrix like given in 3.2.1.
	  * For (i,j), after the index, giving the maximum line capacity.
	  * @param degrees
	  * @return upper triangular matrix
	  */
		
	 public static void upperTriangularL(){
		 L= new int[hIndex][hIndex]; 
		 for(int i=0;i<hIndex;i++) {
			 for(int j=i+1;j<hIndex;j++) {
				 L[i][j]= Math.min(degrees[i], Lsum(i,j+1));
			 }
		 }
	 }
	 
	 /**
	  * C; upper triangular matrix like given in 3.2.1.
	  * For (i,j), after the index, giving the maximum column capacity.
	  * @param degrees
	  * @return upper triangular matrix
	  */
		
	 public static int[][] upperTriangularC(){
		 C= new int[hIndex][hIndex]; 	
		 for(int i=0;i<hIndex;i++) {
			 for(int j=i+1;j<hIndex;j++) {
				 C[i][j]= Math.min(degrees[j], Csum(i+1,j));
			 }
		 }
		 return C;
	 }
		
	 public static int Lsum(int i, int j) {
		 int sum=0;
		 for(int k=j;k<hIndex;k++) {
			 sum=sum+max[i][k];
		 }
		 return sum;
	 }
		
	 public static int Csum(int i, int j) {
		 int sum=0;
		 for(int k=i;k<hIndex;k++) {
			 sum=sum+max[k][j];
		 }
		 return sum;
	 }
	 		
	 /**
	  * Theorem 3.2.1
	  * @param degrees degree sequence
	  */
		
	 public static void maximalMatrix() {
		 max= new int[hIndex][hIndex];
		 for(int i=0;i<hIndex;i++) {
			 for(int j=0; j<hIndex;j++) {
				 int di= degrees[i];
				 int dj= degrees[j];	
				 if(i==j) {
					 max[i][j]=0;
				 }else {
					 if(di!=dj) {
						 max[i][j]=Math.min(di, dj);
					 }else if (di==dj && i!=j) {
						 max[i][j]=(di-1);
					 }
				 }
			 }
		 }
	 }
		
	 public static void genStrip(int[] degreeList) throws IOException, CloneNotSupportedException, CDKException {
		 int[][] A   = new int[matrixSize][matrixSize];
		 degrees=degreeList;
		 flag=true;
		 maximalMatrix();
		 upperTriangularL();
		 upperTriangularC();
		 int[] indices= new int[2];
		 indices[0]=0;
		 indices[1]=1;
		 boolean callForward=true;
		 r=0;
		 y=findY(r);
		 z=findZ(r);
		 while(flag) {
			 nextStep(A,indices,callForward);
	 		 if(!flag) {
	 			 break;
	 		 }
	 		 if(connectLernen) {
	 			indices=connect;
	 			findR(indices);
	 			clearFormers(false, findY(r));
	 			connectLernen=false;
	 			callForward=false;
	 		 }else {
	 			if(lernen) {
	 				indices=successor(lernenIndices,max.length);
	 				findR(indices);
	 				lernen=false;
	 				callForward=false;
	 			}else {
	 				indices=successor(indices,max.length);
	 				updateR(indices);
	 				callForward=true; 
	 			}
	 		 }
		 }
	}
	
    public static int[] successor(int[] indices, int size) {
    	int i0= indices[0];
    	int i1= indices[1];
    	if(i1<(size-1)) {
    		 indices[0]=i0;
    		 indices[1]=(i1+1);
    	}else if(i0<(size-2) && i1==(size-1)) {
    		 indices[0]=(i0+1);
    		 indices[1]=(i0+2);
    	}
    	return indices;
    }
    
    public static int[] predecessor(int[] indices, int size) {
    	int i0= indices[0];
   	 	int i1= indices[1];
   	 	if(i0==i1-1) {
   	 		indices[0]=(i0-1);
   	 		indices[1]=(size-1);
   	 	}else {
   	 		indices[0]=i0;
   	 		indices[1]=(i1-1);
   	 	}
   	 	return indices;
    }
    
	public static int[][] nextStep(int[][] A, int[] indices,boolean callForward) throws IOException, CloneNotSupportedException, CDKException{
		if(callForward) {
			return forwardRow(A, indices,callForward);
		}else {
			return backward(A,indices,callForward);
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
			hydrogen=initialDegrees[i]-sum(A[i]);
			for(int j=hIndex;j<(hIndex+hydrogen);j++) {
				A[i][j]=1;
				A[j][i]=1;
			}
			if(hydrogen!=0) {
				hIndex=hIndex+hydrogen;
			}
		}
		return A;
	}
	
	/**
	 * In backward function, updating the block index.
	 * @param r int block index
	 * @param indices ArrayList<Integer> atom valences
	 * @return
	 */
	 
	 public static void updateR(int[] indices) {
		 int y=findY(r);
		 int z=findZ(r);
		 if(indices[0]<y) {
			 r--;
		 }else if(indices[0]>z) {
			 r++;
		 }
	 }
	 
	 public static void findR(int[] indices) {
		 int block=0;
		 int index=0;
		 int rowIndex= indices[0];
		 for(Integer part: initialPartition) {
			 if(index<=rowIndex && rowIndex<(index+part)) {
				 break;
			 }
			 block++;
			 index=index+part;
		 }
		 r=block;
	 }
	 /**
	  * The third line of the backward method in Grund 3.2.3. The criteria 
	  * to decide which function is needed: forward or backward.
	  * 
	  * @param x 	the value in the adjacency matrix A[i][j]
	  * @param lInverse   lInverse value of indices {i,j}
	  * @param l
	  * @return
	  */
	 
	 public static boolean backwardCriteria(int x, int lInverse, int l) {
		boolean check=false;
		int newX= (x-1);
		if((lInverse-newX)<=l) {
			check=true;
		}
		return check;
	 }
	 
	/**
	  * The third step in Grund 3.2.3.
	  * 
	  * Backward step in the algorithm.
	  * 
	  * @param A 			int[][] adjacency matrix
	  * @param indices 		ArrayList<Integer> indices
	  * @throws IOException 
	  * @throws CloneNotSupportedException
	  * @throws CDKException
	  */
	 
	 public static int[][] backward(int[][] A, int[] indices, boolean callForward) throws IOException, CloneNotSupportedException, CDKException {
		int i=indices[0];
		int j=indices[1];			
		if(i==0 && j==1) {
			flag=false;
			return A;
		}else {
			indices=predecessor(indices, max.length);
			updateR(indices);
			i= indices[0];
			j= indices[1];
			int x= A[i][j];
			int l2= LInverse(degrees,i,j,A);
			int c2= CInverse(degrees,i,j,A);
			
			/**
			 * I changed in backcriteria from x to x-1 but then I had error c6h6 was 98 not 217
			 */
			if(x>0 && (backwardCriteria((x),l2,L[i][j]) && backwardCriteria((x),c2,C[i][j]))){
				A[i][j]=(x-1);
				A[j][i]=(x-1);
				indices = successor(indices,max.length);
				updateR(indices);
				callForward=true;
				return nextStep(A,indices,callForward);
			}else {
				callForward=false;
				return nextStep(A,indices,callForward);
			}
		}
	 }
	 
	public static int[][] forwardRow(int[][] A, int[] indices, boolean callForward) throws IOException, CloneNotSupportedException, CDKException {
		int i=indices[0];
		int j=indices[1];
		int lInverse= LInverse(degrees,i,j,A);
		int cInverse= CInverse(degrees,i,j,A);
		int minimal = Math.min(max[i][j],Math.min(lInverse,cInverse));
		 
		
		/**
		 * First step in the forward method.
		 */
		 
		int maximumValue = forwardMaximal(minimal, lInverse, L[i][j], cInverse, C[i][j]); 
		callForward=true;
		if(j==(max.length-1)) {
			int[][] newMat = forwardSubRow(lInverse, cInverse, maximumValue, i, j,A,indices,callForward);
			return newMat;
		}else {
			return forwardSubRow(lInverse, cInverse, maximumValue, i, j,A,indices,callForward); 
		}  
	} 
	
	
	public static int[][] forwardSubRow(int lInverse, int cInverse, int maximalX,int i, int j, int[][] A, int[] indices, boolean callForward) throws CloneNotSupportedException, CDKException, IOException {
		 if(((lInverse-maximalX)<=L[i][j]) && ((cInverse-maximalX)<=C[i][j])) {
			 A[i][j]=maximalX;
			 A[j][i]=maximalX;
			 if(i==(max.length-2) && j==(max.length-1)) {
				if(canonicalTest(A)) {
					int[][] mat2= new int[A.length][A.length]; 
					for(int k=0;k<A.length;k++) {
						for(int l=0;l<A.length;l++) {
							mat2[k][l]=A[k][l];
						}
					}
					if(connectivityTest(mat2)){
						//IAtomContainer mol= buildC(addHydrogens(mat2,hIndex));
						//outFile.write(mol);
						count++;
						callForward=false;
						return nextStep(A, indices, callForward);
					}else {
						connectLernen=true;
						return A;
					}
				}else {
					if(lernen) {
						return A;
					}else {
						callForward=false;
						return nextStep(A, indices, callForward);
					}
				}
			}else {	
				if(indices[0]==findZ(r) && indices[1]==(max.length-1)) {
					callForward=canonicalTest(A);
					if(callForward) {
						indices=successor(indices,max.length);
						updateR(indices);
						return nextStep(A, indices, callForward);
					}else {
						if(lernen) {
							return A;
						}else {
							return nextStep(A, indices, callForward);
						}
					}
			    }else {
			    	 if(i<findZ(r) && j==(max.length-1)) {
			    		 return A;
			    	 }else {
			    		 indices=successor(indices,max.length);
			    		 updateR(indices);
			    		 return nextStep(A, indices, true);
			    	 }
			    }
			}
		 }else {
			 callForward=false;
			 return nextStep(A, indices, callForward);
		 }
	 }
	
	 /**
	  * First step in the forward method.
	  * @param min minimum of L, C amd maximal matrices for {i,j} indices.
	  * @param lInverse Linverse value of {i,j}
	  * @param l        L value of {i,j}
	  * @param cInverse Cinverse value of {i,j}
	  * @param c		C value of {i,j}
	  * @return int max
	  */
	 
	 public static int forwardMaximal(int min, int lInverse, int l, int cInverse, int c) {
		 int max=0;
		 for(int v=min;v>=0;v--) {
			 if(((lInverse-v)<=l) && ((cInverse-v)<=c)) {
				 max=max+v;
				 break;
			 }
		 }
		 return max;
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
	
	public static int LInverse(int[] degrees, int i, int j, int[][]A) {
		int sum=0;
		for(int s=0;s<j;s++) {
			sum=sum+A[i][s];
		}
		return degrees[i]-sum;
	}
		 
	public static int CInverse(int[] degrees, int i, int j, int[][]A) {
		int sum=0;
		for(int s=0;s<i;s++) {
			sum=sum+A[s][j];
		}
		return degrees[j]-sum;
	}
	
	 /**
	  * Main functions
	  */
	 
	 public static ArrayList<Integer> getPartition(int[] degrees, ArrayList<Integer> partition){
    	 ArrayList<Integer> newPartition = new ArrayList<Integer>();
		 int i=0;
    	 for(Integer p:partition) {
    		 Integer[] subArray= getBlocks(degrees,i,p+i);
    		 newPartition.addAll(getSubPartition(subArray));
    		 i=i+p; 
    	 }
    	 return newPartition;
	 }
	 
	
	 public static ArrayList<Integer> getSubPartition(Integer[] degrees){
		 ArrayList<Integer> partition = new ArrayList<Integer>();
		 int i=0;
	     int size= degrees.length;
	     int count=0;
	     int next=0;
		 while(i<size) {
			 count=nextCount(i,size,degrees,partition);
			 next=(i+count);
			 if(next==size) {
				 break;
			 }else {
				 i=next;
			 }
	     }  
	     return partition;
	 }
	 
	 public static int nextCount(int i, int size, Integer[] degrees, ArrayList<Integer> partition) {
		 int count=1;
		 if(i==(size-1)) {
			 partition.add(1);
		 }else {
			 for(int j = i+1; j < size; j++){  
	         	 if(degrees[i] == degrees[j]){  
	        		 count++;  
	        		 if(j==(size-1)){
	        			 partition.add(count);
	        			 break;
	        		 }
	             }else{
	            	 partition.add(count);
	            	 break;
	             }
	         }   
		 }
		 return count;
	 }
	 public static void run() throws IOException, CDKException, CloneNotSupportedException {
		 long startTime = System.nanoTime(); 
		 if(verbose) System.out.println("MORGEN is generating isomers of "+formula+"...");
		 getSymbolsOccurrences(formula);
		 initialDegrees();
		 //build(formula);
		 outFile = new SDFWriter(new FileWriter(filedir+"output.sdf"));
		 canonicalBlockbasedGenerator();
		 if(verbose) System.out.println("The number of structures is: "+count);
		 outFile.close();
		 long endTime = System.nanoTime() - startTime;
		 double seconds = (double) endTime / 1000000000.0;
	     DecimalFormat d = new DecimalFormat(".###");
	     if(verbose) System.out.println("Time: "+d.format(seconds));
	}
	 
	 /**
	  * The main function for the block based canonical matrix generation.
	  * 
	  * @param molecularFormula String molecular formula
	  * @param degrees			ArrayList<Integer> initial degree set. 
	  * @param newDegrees		ArrayList<Integer> new degree set after hydrogen partition.
	  * @param numberOfAtoms	int total number of atoms.
	  	@param hydrogenIndex	int the beginning index of hydrogen atoms in A.
	  * @param filePath 	    String path to store the output files.
	  * @throws IOException
	  * @throws CloneNotSupportedException
	  * @throws CDKException
	  */
		
	 public static List<int[]> distributeHydrogens(ArrayList<Integer> partition, int[] degrees) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException{
		 List<int[]> degreeList= new ArrayList<int[]>();
		 List<int[]> distributions= HydrogenDistributor.run(partition,degrees);
		 for(int[] dist: distributions) {
			 int[] newDegree= new int[size];
			 for(int i=0;i<size;i++) {
				 newDegree[i]=(degrees[i]-dist[i]);
			 }
			 degreeList.add(newDegree);
		 }
		 return degreeList;
	 }
	 
	 public static int[] opc;
	 public static void canonicalBlockbasedGenerator() throws IOException, CloneNotSupportedException, CDKException{
		size=initialDegrees.length;
		List<int[]> newDegrees= distributeHydrogens(firstOccurrences, firstDegrees);
		lernenIndices= new int[2];
		lernen=false;
		connectLernen=false;
		for(int[] degree: newDegrees) {
			lernenIndices= new int[2];
			connect= new int[2];
			connectLernen=false;
			lernen=false;
			partitionList.clear();
			formerPermutations.clear();
			initialPartition=getPartition(degree,occurrences);
			partitionList.add(0,initialPartition);
			genStrip(degree);
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
		 
	 public static Set<Integer> wValues(Set<Integer> nValues, int[] Kformer){
		 Set<Integer> wValues= new HashSet<Integer>();
		 for(Integer i:nValues) {
			 wValues.add(Kformer[i]);
		 }
		 return wValues;
	 }
		 
	 /**
	  * Finding the K values as explained in 3.6.2.
	  * @param wValues Set<Integer> wValues  
	  * @param kFormer ArrayList<Integer> the K values of the former step
	  * @return ArrayList<Integer>
	  */
		 
	 public static int[] kValues(int total, Set<Integer> wValues, int[] kFormer){
		 int[] kValues= new int[total];
		 int min= Collections.min(wValues);
		 for(int i=0;i<total;i++) {
			 if(wValues.contains(kFormer[i])) {
				 kValues[i]=min;
			 }else {
				 kValues[i]=kFormer[i];
			 }
		 }
		 return kValues;
	 }
	 
	 /**
	  * Initializing the first K values list.
	  * @param total int number of atoms.
	  * @return ArrayList<Integer>
	  */
		 
	 public static int[] initialKList(int total){
		 int[] k= new int[total];
		 for(int i=0;i<total;i++) {
			 k[i]=i;
		 }
		 return k;
	 }
	 
	 /**
	  * Test whether an adjacency matrix is connected or disconnected.
	  * @param p int the index where the atoms with degree 1 starts.
	  * @param mat int[][] adjacency matrix
	  * @return boolean
	  */
		
	 public static int[] connect=new int[2];
	 public static boolean connectLernen=false; 
	 public static boolean connectivityTest(int[][] mat) {
		 connectLernen=false;
		 boolean check=false;
		 int[] kValues=initialKList(hIndex);
		 Set<Integer> nValues= new HashSet<Integer>();
		 Set<Integer> wValues= new HashSet<Integer>();
		 int zValue= 0;
		 for(int i=0;i<hIndex;i++) {
			 nValues= nValues(i, hIndex, mat);
			 wValues=wValues(nValues,kValues);
			 zValue = Collections.min(wValues);
			 kValues= kValues(hIndex, wValues, kValues);
		 }
		 if(zValue==0 && allIs0(kValues)) {
			 check=true;
		 }else {
			 setLernenConnectivity(mat, kValues);
		 }
		 return check;
	 }
	 
	 public static void setLernenConnectivity(int[][] mat, int[] kValues) {
		 int value=kValues[hIndex-1];
		 for(int i=hIndex-2;i>=0;i--) {
			 if(kValues[i]!=value) {
				 connectLernen=true;
				 connect[0]=i;
				 connect[1]=hIndex-1;
				 break;
			 }
		 }
	 }
	 
	 public static int findNonZeroEntry(int[] row) {
		 int index=0;
		 int lastEntry=row[hIndex-1];
		 for(int i=hIndex-2;i>0;i--) {
			 if(row[i]!=lastEntry) {
				 index=i;
				 break;
			 }
		 }
		 return index;
	 }
	 /**
	  * Checks whether all the entries are equal to 0 or not. (Grund 3.6.2)
	  * @param list ArrayList<Integer>
	  * @return boolean
	  */
		 
	 public static boolean allIs0(int[] list) {
		 boolean check=true;
		 for(int i=0;i<list.length;i++) {
			 if(list[i]!=0) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 	
	 /**
	  * Canonical Test Functions
	  */
		
	 public static int findY(int r) {
		 return (sum(initialPartition,(r-1)));
	 }
		
	 public static int findZ(int r) {
		 return (sum(initialPartition,r)-1);
	 }
	 
	 public static boolean canonicalTest(int[][] matrix) throws IOException, CloneNotSupportedException, CDKException {
		lernen=false;
		if(blockTest(r,matrix)) {
			return true;
		}else {
			return false;
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
	 
	public static boolean blockTest(int r, int[][] A) {
		boolean check=true;
		int y=findY(r);
		int z=findZ(r);
		for(int i=y;i<=z;i++) {
			if(!blockTestPerm(i, r, A, partitionList.get(i),canonicalPartition(i,partitionList.get(i)))){	
				check=false;
				break;
			}
		}
		clearFormers(check, y);
		return check;
	}
	
	public static void clearFormers(boolean check, int y) {
		 if(check==false) {
			 for(int i=y;i<formerPermutations.size();i++) {
				 formerPermutations.get(i).removeAll(formerPermutations.get(i));
			 }
			 for(int i=y+1;i<partitionList.size();i++) {
				 partitionList.get(i).removeAll(partitionList.get(i));
			 }
		 }
	}
	
	public static void candidatePermutations(int index, int y, int total, List<Permutation> cycles) {
		 List<Permutation> newList= new ArrayList<Permutation>();
		 newList.addAll(cycles);
		 formerPermutations.add(index,newList);
		 if(index!=0) {
			 List<Permutation> formers = formerPermutations.get(index-1); 
			 for(Permutation form: formers) {
				 if(!form.isIdentity()) {
					 formerPermutations.get(index).add(form);
				 }
			 }
			 List<Permutation> newForm = new ArrayList<Permutation>();
			 for(Permutation frm: formers) {
				 if(!frm.isIdentity()) {
					 newForm.add(frm);
				 }
			 }
			 List<Permutation> newCycles = new ArrayList<Permutation>();
			 for(Permutation cyc: cycles) {
				 if(!cyc.isIdentity()) {
					 newCycles.add(cyc);
				 }
			 }
			 for(Permutation perm: newForm) {
				 for(Permutation cycle: newCycles) {
					 Permutation newPermutation =cycle.multiply(perm);
					 if(!newPermutation.isIdentity()) {
						 formerPermutations.get(index).add(newPermutation); 
					 }
				 }
			 }
		 }
	 }
	
	public static boolean blockTestPerm(int index, int r,int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition) {
		 boolean check= true;
		 if(!firstCheck(index, A, newPartition)){
		 	check=false;
	     }else {
	    	 int y = findY(r);
	    	 int total = sum(partition);
	    	 List<Permutation> cycleTrans=cycleTranspositions(index, partition);
	    	 candidatePermutations(index,y,total,cycleTrans);		 
	    	 boolean cycleCheck=check(index, y, total, A, partition, newPartition); 
	    	 if(!cycleCheck) {
	    		 check=false;
	    	 }else {
	    		 addPartition(index, newPartition, A);  
	    	 }
	     }
		 return check;
	 }
	
	public static boolean allis1(ArrayList<Integer> partition) {
		 boolean check=true;
		 for(int i=0;i<partition.size();i++) {
			 if(partition.get(i)!=1) {
				 check=false;
				 break;
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
		 ArrayList<Integer> refinedPartition= new ArrayList<Integer>();
		 if(allis1(newPartition)) {
			 refinedPartition=newPartition;
		 }else {
			 refinedPartition= refinedPartitioning(newPartition,A[index]);
		 }
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
	 
	 public static int[] row2compare(int index, int[][] A, Permutation cycleTransposition) {
		 int[] array = cloneArray(A[findIndex(index, cycleTransposition)]);
		 array=actArray(array, cycleTransposition);
		 return array;
	 }
	 
	 /**
	  * Checking whether the modified and original rows are equal or not.
	  * @param cycleTransposition Permutation - cycle transposition
	  * @param partition atom partition
	  * @param index row index
	  * @param y beginning index of a block
	  * @param A adjacency matrix
	  * @param testPermutation Permutation
	  * @return
	  */
	 	 
	 public static int findIndex(int index, Permutation cycle, Permutation former) {
		 int size= cycle.size();
		 int output=0;
		 Permutation perm= cycle.multiply(former);		 
		 for(int i=0;i<size;i++) {
			 if(perm.get(i)==index) {
				 output=i;
				 break;
			 }
		 }
		 return output;
	 }
	 
	 public static int findIndex(int index, Permutation cycle) {
		 int size= cycle.size();
		 int output=0;
		 for(int i=0;i<size;i++) {
			 if(cycle.get(i)==index) {
				 output=i;
				 break;
			 }
		 }
		 return output;
	 }
	 
	 public static int[] row2compare(int index, int[][] A, Permutation cycleTransposition, Permutation formerPermutation) {
		 int[] array = cloneArray(A[findIndex(index, cycleTransposition,formerPermutation)]);
		 Permutation perm = formerPermutation.multiply(cycleTransposition);
		 array= actArray(array,perm);
		 return array;
	 }
	 
	 public static int[] cloneArray(int[] array) {
		 int length= array.length;
		 int[] cloned = new int[length];
		 for(int i=0;i<length;i++) {
			 cloned[i]=array[i];
		 }
		 return cloned;
	 }
	 
	 public static 	Permutation getCanonicalPermutation(int[] max, int[] check,ArrayList<Integer> partition) {
    	 //set size 
		 int[] cycles= getCyclePermutation(max, check, partition);
    	 int[] perm= new int[size];
		 for(int i=0;i<size;i++) {
			 for(int j=0;j<size;j++) {
				 if(i==cycles[j]) {
					 perm[i]=j;
				 } 
			 }
		 }
		 return new Permutation(perm);
     }
	 	 	 
	 public static int[] getCyclePermutation(int[] max, int[] check,ArrayList<Integer> partition) {
    	 int[] values= idValues(sum(partition));
    	 int i=0;
    	 if(!equalSetCheck(max,check,partition)) {
    		 return values;
    	 }else {
    		 for(Integer p:partition) {
    			 Integer[] can= getBlocks(max,i,p+i);
    			 Integer[] non= getBlocks(check,i,p+i);
    			 values = getCyclesList(can, non, i, values);
    			 i=i+p; 
    		 }
        	 return values;
    	 }
     }
	 
	 public static int findMatch(Integer[] max, Integer[] non, int value, int start) {
    	 int size=non.length;
    	 int index=start;
    	 for(int i=start;i<size;i++) {
    		 if(non[i]==value) {
    			 if(max[i]!=non[i]) {
    				 index=i;
        			 break; 
    			 }
    		 }
    	 }
    	 return index;
     }
	     
	 public static int[] getCyclesList(Integer[] max, Integer[] non, int index, int[] values) {
    	 int i=0;
    	 int permutationIndex=0;    	 
    	 while(i<max.length && max[i]!=0) {
    		 if(max[i]!=non[i]) {
        		 permutationIndex = findMatch(max,non, max[i],i);
        		 if(i!=permutationIndex) {
        			 non=permuteArray(non, i, permutationIndex);
        		 }
        		 int temp=values[i+index];
        		 values[i+index]=values[permutationIndex+index];
        		 values[permutationIndex+index]=temp;
    		 }
    		 i++;
    	 }
    	 return values;
     }
	 
	 public static Permutation getEqualPerm(Permutation cycleTransposition, int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition) {
		 int[] check=row2compare(index, A, cycleTransposition);
		 Permutation canonicalPermutation = getCanonicalPermutation(A[index],check,newPartition);
		 return canonicalPermutation;
	 }
	 
	 	 
	public static Permutation getCanonicalCycle(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition,Permutation cycleTransposition) {
		 biggest=true;
		 Permutation canonicalPermutation = idPermutation(total);
		 if(!equalBlockCheck(newPartition,index,y,A,cycleTransposition, canonicalPermutation)) {
			 canonicalPermutation = getEqualPerm(cycleTransposition, index, y, total,A, partition, newPartition);
			 int[] check= row2compare(index, A, cycleTransposition);
			 check = actArray(check,canonicalPermutation);
		 }			 
		 return canonicalPermutation;
	 }
	
	public static boolean check(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition) {
		 boolean check=true;
		 List<Permutation> formerList= new ArrayList<Permutation>();
		 for(Permutation permutation:formerPermutations.get(index)) { 
			 biggerCheck(index, A, permutation, newPartition);
			 if(biggest) {
				 Permutation canonicalPermutation = getCanonicalCycle(index, y, total, A, partition, newPartition, permutation);
				 int[] test=row2compare(index, A, permutation);
				 test= actArray(test,canonicalPermutation);	
				 //if(descendingOrderUpperMatrix(index, newPartition, A[index], test)) {
				 if(descendingOrdercomparison(newPartition,A[index], test)) {	
				 	 if(canonicalPermutation.isIdentity()) {
						 if(equalCheck(A[index],test,partition)) {
							 formerList.add(permutation); 
						 }
					 }else {
						 Permutation newPermutation=canonicalPermutation.multiply(permutation);
						 formerList.add(newPermutation);
					 }
				 }else {
					 formerList.clear();
					 check=false;
					 break;
				 }	  
			 }else {
				 formerList.clear();
				 check=false;
				 break;
			 }	  
		 }
		 if(check) {
			 formerPermutations.get(index).clear();
			 formerPermutations.set(index, formerList);
		 }
		 return check;
	 }
	public static List<Permutation> cycleTranspositions(int index, ArrayList<Integer> partition) {
		 int total=sum(partition);
		 List<Permutation> perms= new ArrayList<Permutation>();
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
	  * To get the canonical partition like in Grund Thesis 3.3.11
	  * @param i int row index
	  * @param partition ArrayList<Integer> partition
	  * @return
	  */
	 
	public static ArrayList<Integer> canonicalPartition(int i, ArrayList<Integer> partition){
		 return partitionCriteria(partition,i+1);
	}
	 
	 /**
	  * add number of 1s into a ArrayList
	  */
		 
	 public static ArrayList<Integer> addOnes(ArrayList<Integer> list, int number) {
		 for(int i=0;i<number;i++) {
			 list.add(1);
		 }
		 return list;
	 }
	 
	/**
	  *	Grund Thesis 3.3.2 Partitioning criteria (DONE)  
	  * @param partEx the former partition
	  * @param degree degree of the partitioning.
	  * @return
	  */
		 
	 public static ArrayList<Integer> partitionCriteria(ArrayList<Integer> partEx, int degree){
		 ArrayList<Integer> partNew = new ArrayList<Integer>();
		 partNew=addOnes(partNew,degree);
		 /**
		  * I had (degree-1)
		  */
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
	 
	 private void parseArgs(String[] args) throws ParseException, IOException, org.apache.commons.cli.ParseException{
		 Options options = setupOptions(args);	
		 CommandLineParser parser = new DefaultParser();
		 try {
			 CommandLine cmd = parser.parse(options, args);
			 MORGEN.formula = cmd.getOptionValue("formula");
			 MORGEN.filedir = cmd.getOptionValue("filedir");				 
			 if (cmd.hasOption("verbose")) MORGEN.verbose = true;		
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
			 formatter.printHelp( "java -jar MORGEN.jar", header, options, footer, true );
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
	
	public static void main(String[] args) throws IOException, CDKException, CloneNotSupportedException {	
		MORGEN gen = new MORGEN();
		//String[] arg= {"-f", "C10H16", "-d", "C:\\Users\\mehme\\Desktop\\", "-v"};
		try {
			gen.parseArgs(args);
			MORGEN.run();
		} catch (Exception e) {
			if (MORGEN.verbose) e.getCause(); 
		}
	}
}
