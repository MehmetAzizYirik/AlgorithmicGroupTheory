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
 * and Chapter 5 for molecular structure generation.And also Grund's thesis. [2]
 * 
 * [1] Kerber, A., Laue, R., Meringer, M., RÃ¼cker, C. and Schymanski, E., 2013.
 * Mathematical chemistry and chemoinformatics: structure generation, elucidation
 * and quantitative structure-property relationships. Walter de Gruyter.
 * 
 * [2] Grund, Roland, and Reinhard Müller. Konstruktion molekularer Graphen mit 
 * gegebenen Hybridisierungen und überlappungsfreien Fragmenten. Lehrstuhl II für 
 * Mathematik, 1995.
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
import java.util.stream.IntStream;

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
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
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
		 
	/**
	 * These functions for the duplicate free canonical matrix generation as given in Grund 
	 * thesis 3.3.13. There are still some duplicates as mentioned also at the thesis. These
	 * functions are still under development and will be tested.  
	 */
	
	public static void canonicalMatrixGenerator(ArrayList<Integer> degrees, ArrayList<Integer> partition) throws IOException{
		int size    = degrees.size();
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] L2  = upperTriangularL2(degrees,A); //Where and how these L2 and C2
		int[][] C   = upperTriangularC(degrees);
		int[][] C2  = upperTriangularC2(degrees,A);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forwardCanonical(degrees,partition,A,max,L,C,indices); //It was originally without partition entry.
	}
	
	public static void forwardCanonical(ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
	 	int minimal= Math.min(max[i][j],Math.min(l2,c2));
	 	for(int h=minimal;h>=0;h--) {
	 		if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
	 			A[i][j]=A[j][i]=h;
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backwardCanonical(degrees,partition,A, max, L, C, indices);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				if(modified.get(0)>i && j==A.length-1) { 
	 					partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					if(canonicalRowCheck(A[i],i,partition,A.length)) {
	 						partition=refinedPartitioning(partition,A[i]);
	 						forwardCanonical(degrees, partition, A, max, L, C, modified);
	 					}
	 				}else {
	 					forwardCanonical(degrees, partition, A, max, L, C, modified);
	 				}
	 			}
	 		}else {
	 			backwardCanonical(degrees, partition,A, max, L, C, indices);
	 		}
	 	}
	}
	
	public static boolean desCheck=true; 
	 
	 /**
	  * Refined Partition - Grund 3.3.11
	  * Create the refined partition from the new row. The occurences of
	  * the numbers in the blocks of the row. 
	  * @param partition Partition
	  * @param row int[] row of a matrix to check canonicality. 
	  * @return
	  */
	 
	 public static ArrayList<Integer> refinedPartitioning(ArrayList<Integer> partition, int[] row){
		 ArrayList<Integer> refined= new ArrayList<Integer>();
		 int index=0;
		 int count=1;
		 for(Integer p:partition) {
			 for(int i=index;i<p+index;i++) {
				 if(i+1<p+index) {
					 if(row[i]==row[i+1]) {
						 count++;
					 }else if(row[i]>row[i+1]){
						 refined.add(count);
						 count=1;
					 } else {
						 desCheck=false;
						 break;
					 }
				 }else {
					 refined.add(count);
					 count=1;
				 }
			 }
			 index=index+p;
		 }
		 return refined;
	 }
	 
	/**
	 * Canonical row check. To check the row is in descending order or not.
	 **/
	
	public static boolean canonicalRowCheck(int[] row,int i, ArrayList<Integer> partition, int total) {
		boolean check=true;
		PermutationGroup group=getYoungGroup(partition,total);
		for(Permutation perm: group.all()) {
			if(!desBlockwiseCheck(partition,row,actArray(row,perm))) {
				check=false;
				break;
			}
		}
		return check;
	}
	
	/**
	  * Blockwise descending order check of riginal row and its canonical form.
	  * The blocks are repsented as integers to check the descending order.
	  * @param partition
	  * @param canonical
	  * @param original
	  * @return boolean true if descending order.
	  */
	 
	 public static boolean desBlockwiseCheck(ArrayList<Integer> partition, int[] canonical, int[] original){
		 boolean check=true;
		 int index=0;
		 for(Integer p:partition) {
			 if(toInt(getBlocks(canonical,index,p+index))<toInt(getBlocks(original,index,p+index))) {
				 check=false;
				 break;
			 }
			 index=index+p;
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
	 
	 public static Integer [] getBlocks(int[] row, int begin, int end) {
		 Integer[] sub=subArray(row,begin,end);
		 return sub;
	 }
	 
	public static Integer[] subArray(int[] array, int begin,int end) {
		return IntStream.range(begin, end).mapToObj(i->array[i]).toArray(Integer[]::new);
	}
	public static int[] actArray(int[] strip, Permutation p) {
		int length= strip.length;
		int[] modified = new int[length];
		for(int i=0; i<length;i++) {
			modified[p.get(i)]=strip[i]; 
		}
		return modified;
	}
	/**
	 * Youngsubgroup generation based on the generator permutations. These permutations are built for
	 * each sub partitions
	 * @param partition atom partition
	 * @param size number of atoms
	 * @return
	 */
	public static PermutationGroup getYoungGroup(ArrayList<Integer> partition, int size) {
		return generateGroup(getYoungsubgroupGenerators(partition,size),size); 
	}
	
	/**
	 * Generating a group from a list permutations
	 * @param perm a list of permutations	
	 * @return permutation group generated by a list of permutations
	 */
	
	public static PermutationGroup generateGroup(List<Permutation> generators,int size) {
        return new PermutationGroup(size, generators); 
	}
	
	/**
	 * Getting the list of generator permutations for all the sub partitions.
	 * @param partition Atom partition
	 * @param total number of atoms
	 * @return
	 */
	
	public static List<Permutation> getYoungsubgroupGenerators(ArrayList<Integer> partition,int total){
		ArrayList<ArrayList<Integer>> parts=subPartitionsList(partition);
		List<Permutation> perms= new ArrayList<Permutation>();
		for(ArrayList<Integer> part:parts) {
			perms.addAll(getYoungsubgroupGenerator(part,total));
		}
		return perms;
	}
	
	/**
	 * Example: For a set {5,6,7} we need to get S3. Generators are (12) and (123).
	 * This 1 2 3 are the indices of the set elements. So we need to replace them 
	 * with the values. The generators also should start with 0 1 2 until 5. So
	 * If the size of the total partition is 10; the generators are
	 * 
	 * [0,1,2,3,4,6,5,7,8,9]
	 * [0,1,2,3,4,6,7,5,8,9]
	 *  
	 * @param subPart sub partition
	 * @param total partition size
	 * @return
	 */
	
	public static List<Permutation> getYoungsubgroupGenerator(ArrayList<Integer> subPart, int total){
		List<Permutation> perms= new ArrayList<Permutation>();
		int setSize= subPart.size();
		int[] id= idValues(total);
		if(setSize==1) {
			perms.add(new Permutation(id));
		}else if(setSize==2) {
			id[subPart.get(0)]   = subPart.get(0)+1;
			id[subPart.get(0)+1] = subPart.get(0);
			perms.add(new Permutation(id));
		}else {
			
			/**
			 * The first generator with length 2 which is (1,2)
			 */
			
			int firstIndex= subPart.get(0);
			id[firstIndex]  = firstIndex+1;
			id[firstIndex+1]= firstIndex;
			Permutation gen1= new Permutation(id);
			perms.add(gen1);
			
			/**
			 * The second generator with the full ( 1,2,..,n)
			 */
			
			Permutation gen2= new Permutation(fullValues(total,firstIndex,firstIndex+setSize));
			perms.add(gen2);
			
		}
		return perms;
	}
	
	/**
	 * Values for the second generator permutation - with the full ( 1,2,..,n)
	 * @param id the id values
	 * @param begin first value in the sub partition
	 * @param end	the last value in the sub partition
	 * @return
	 */
	
	public static int[] fullValues(int total, int begin, int end) {
		int[] id=idValues(total);
		for(int i=begin;i<end-1;i++) {
			id[i]=i+1;
		}
		id[end-1]=begin;
		return id;
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
	
	public static ArrayList<ArrayList<Integer>> subPartitionsList(ArrayList<Integer> partition){
		 ArrayList<ArrayList<Integer>> parts = new ArrayList<ArrayList<Integer>>();
		 for(int i=0;i<partition.size();i++) {
			parts.add(subPartition(partition,i)); 
		 }
		 return parts;
	 }
	
	/**
	  * Since permutations starts with 0, the subpartitions should also start with 0.
	  * @param partition
	  * @param index
	  * @return
	  */
	 
	 public static ArrayList<Integer> subPartition(ArrayList<Integer> partition, int index){
		 ArrayList<Integer> sub= new ArrayList<Integer>();
		 int former = sum(partition,index-1)+1; //TODO: This function is not coded in a proper way. Check later.
		 int latter = sum(partition,index);
		 if(former==latter) {
			 sub.add(former-1);
		 }else {
			 for(int i=former-1;i<=latter-1;i++) {
				 sub.add(i);
			 }
		 }
		 return sub;
	 }
	 
	 /**
	  * Sum values until an index
	  * @param list
	  * @param i
	  */
	public static int sum(ArrayList<Integer> list, int index) {
		int sum=0;
		for(int i=0;i<=index;i++) {
			sum=sum+list.get(i);
		}
		return sum;
	}
	
	public static ArrayList<Integer> canonicalPartition(int i, ArrayList<Integer> partition){
		if(i==0) {
			return partitionWDegree(partition,2);
		}else {
			return partitionWDegree(partition,1);
		}
	}
	
	/**
	  * Grund 3.3.2: Repartitioning a partition  (of degree zero) for a given degree.  (DONE)
	  * @param partEx partition with zero 
	  * @param degree repartitioning degree
	  * @return
	  */
	 
	 public static ArrayList<Integer> partitionWDegree(ArrayList<Integer> partEx, int degree){
		 for(int i=1;i<=degree;i++) {
			 partEx=partitionCriteria(partEx,i);
		 }
		 return partEx;
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
		 
	public static void backwardCanonical(ArrayList<Integer> degrees,ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
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
			System.out.println(Arrays.deepToString(A));
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					ArrayList<Integer> modified2=successor(modified,max.length);
	 				forwardCanonical(degrees, partition, A, max, L, C, modified2);
				}else {
					backwardCanonical(degrees,partition, A, max, L, C, modified);
				}
			}
		}
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
