package AlgorithmicGroupTheory;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;

public class canonicalBlockwiseTest extends groupFunctions {
	public static int size=0;
	public static int hIndex=0;
	public static int hoppa=0;
	public static ArrayList<int[][]> output= new ArrayList<int[][]>();
	public static ArrayList<Integer> initialDegrees= new ArrayList<Integer>();
	public static ArrayList<Integer> inputPartition= new ArrayList<Integer>();	
	public static ArrayList<ArrayList<Integer>> partitionList= new ArrayList<ArrayList<Integer>>();
	public static ArrayList<ArrayList<Permutation>> representatives = new ArrayList<ArrayList<Permutation>>();
	
	
	public static void canonicalBlockbasedGenerator(ArrayList<Integer> degrees, ArrayList<Integer> partition) throws IOException, CloneNotSupportedException, CDKException{
		long startTime = System.nanoTime(); //Recording the duration time.
		size = 6;
		initialDegrees=degrees;
		inputPartition=partition;
		FileWriter fWriter = new FileWriter("C:\\Users\\mehme\\Desktop\\output.txt");
		PrintWriter pWriter = new PrintWriter(fWriter);
		pWriter.print("Result"+"\n");
		partitionList.add(0,inputPartition);
		gen(degrees,pWriter);
		fWriter.close();
		pWriter.close();
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		System.out.println("time"+" "+d.format(seconds));
	}
	
	public static List<ArrayList<Integer>> distributeHydrogens(ArrayList<Integer> partition, ArrayList<Integer> degrees) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException{
		List<ArrayList<Integer>> degreeList= new ArrayList<ArrayList<Integer>>();
		List<int[]> distributions= run(partition,degrees);
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
		 
	public static void gen(ArrayList<Integer> degrees, PrintWriter pWriter) throws IOException, CloneNotSupportedException, CDKException {
		int r=0;
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] C   = upperTriangularC(degrees);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forwardBuild(r,degrees,inputPartition,A,max,L,C,indices,pWriter);
	}
	
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
	 
	 public static void forwardBuild(int r, ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices, PrintWriter pWriter) throws IOException, CloneNotSupportedException, CDKException {
		 int y=findY(r);
		 int z=findZ(r);
		 pWriter.print("forward build r y z"+" "+r+" "+y+" "+z+"\n");
		 pWriter.print("mat"+" "+Arrays.deepToString(A)+"\n");
		 int i=indices.get(0);
		 int j=indices.get(1);
		 pWriter.print("forward indices"+" "+i+" "+j+"\n");
		 int l2= LInverse(degrees,i,j,A);
		 int c2= CInverse(degrees,i,j,A);
		 int minimal= Math.min(max[i][j],Math.min(l2,c2));
		 pWriter.print("l and c and minimal"+" "+l2+" "+c2+" "+minimal+" "+max[i][j]+"\n");
		 for(int h=minimal;h>=0;h--) {
			 if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
				 pWriter.print("h"+" "+h+" "+A[i][j]+"\n");
		 		 A[i][j]=A[j][i]=h;
		 		 if(i==(max.length-2) && j==(max.length-1)) {
		 			 if(blockTest(i, r, A,pWriter)) {
		 				 backwardBuild(r, degrees,partitionList.get(y),A, max, L, C, indices, pWriter);
		 			 }
		 		 }else {
		 			 ArrayList<Integer> modified=successor(indices,max.length);
		 			 pWriter.print("successor"+" "+modified.get(0)+" "+modified.get(1)+"\n");
		 			 forwardBuild(r, degrees, partitionList.get(y), A, max, L, C, modified, pWriter);
		 		 }
		 	 }else {
		 		 if(i==max.length-2 && j==max.length-1) {
		 			 ArrayList<Integer> modified=predecessor(indices, max.length);
		 			 indices=modified;
		 		 }
		 		 backwardBuild(r, degrees, partitionList.get(y),A, max, L, C, indices, pWriter);
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
	
	public static int updateR(int r, ArrayList<Integer> indices) {
		int y=findY(r);
		if(indices.get(0)<y) {
			return (r-1);
		}else {
			return r;
		}
	}
	
	public static int[][] addHydrogens(int[][] A, int index){
		hIndex=index;
		int hydrogen=0;
		for(int i=0;i<index;i++) {
			hydrogen=initialDegrees.get(i)-sum(A[i]);
			for(int j=hIndex;j<(hIndex+hydrogen);j++) {
				A[i][j]=1;
				A[j][i]=1;
			}
			if(hydrogen==0) {
				hIndex++;
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
	
	public static void backwardBuild(int r, ArrayList<Integer> degrees,ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices, PrintWriter pWriter) throws IOException, CloneNotSupportedException, CDKException {
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
			pWriter.println("done");
			//if(connectivityTest(6,addHydrogens(A,4))){
				System.out.println("done");
				System.out.println(Arrays.deepToString(addHydrogens(A,4)));
				//depict(buildC(addHydrogens(A,4)),"C:\\Users\\mehme\\Desktop\\"+hoppa+".png");
				hoppa++;
			//}
			//System.out.println(Arrays.deepToString(addHydrogens(A,4)));
			pWriter.print(Arrays.deepToString(A)+"\n");
			//pWriter.print("test bu"+" "+y+" "+representatives.size()+" "+partitionList.size()+"\n");
			for(int h=0; h<representatives.size(); h++) {
				pWriter.print(h+" "+representatives.get(h)+"\n");
			}
			for(int k=(1); k<partitionList.size(); k++) {
				pWriter.print(k+" "+partitionList.get(k)+"\n");
			}
			
			for(int h=1; h<representatives.size(); h++) {
				pWriter.print(h+" "+representatives.size()+"\n");
			}
			
			for(int k=(2); k<partitionList.size(); k++) {
				pWriter.print(k+" "+partitionList.size()+"\n");
			}
			for(int h=0; h<representatives.size(); h++) {
				representatives.get(h).removeAll(representatives.get(h));
			}
			for(int k=(1); k<partitionList.size(); k++) {
				partitionList.get(k).removeAll(partitionList.get(k));
			}
			
			for(int h=0; h<representatives.size(); h++) {
				pWriter.print(h+" "+representatives.get(h)+"\n");
			}
			for(int k=(1); k<partitionList.size(); k++) {
				pWriter.print(k+" "+partitionList.get(k)+"\n");
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
					forwardBuild(r, degrees, partitionList.get(y), A, max, L, C, modified2, pWriter);
				}else {
					backwardBuild(r, degrees,partitionList.get(y), A, max, L, C, modified, pWriter);
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
	 
	 public static boolean blockTest(int index, int r, int[][] A, PrintWriter pWriter) {
		 pWriter.print("block test started"+" "+index+" "+r+" "+Arrays.deepToString(A)+"\n");
		 boolean check=true;
		 int y=findY(r);
		 int z=findZ(r);
		 pWriter.print("block"+"\n");
		 pWriter.print("index r y z"+" "+index+" "+r+" "+y+" "+z+"\n");
		 for(int s=0;s<partitionList.size();s++) {
			 pWriter.print("partition"+" "+s+" "+partitionList.get(s)+"\n");
		 }
		 for(int j=0;j<representatives.size();j++) {
			 pWriter.print("representatives"+" "+j+" "+representatives.get(j)+"\n");
		 }
		 for(int i=y;i<=z;i++) {
			 pWriter.print(i+" "+Arrays.deepToString(A)+" "+partitionList.get(i)+" "+canonicalPartition(i,partitionList.get(i))+"\n");
			 if(!blockDemo(i,r, A, partitionList.get(i),canonicalPartition(i,partitionList.get(i)),pWriter)){
				 check=false;
				 break;
			 }
		 }
		 clear(check, y,pWriter);
		 return check;
	 }
	 
	 /**
	  * If the blockTest is false, clear the reps and representatives 
	  * of the tested block. 
	  * @param check boolean blockTest
	  * @param y int beginning index of the block.
	  */
	 
	 public static void clear(boolean check, int y,PrintWriter pWriter) {
		 pWriter.print("clear check"+" "+y+" "+check+"\n");
		 
		 for(int i=y;i<representatives.size();i++) {
			 pWriter.print("reps"+" "+i+" "+representatives.get(i)+"\n");
		 }
		 
		 for(int i=y+1;i<partitionList.size();i++) {
			 pWriter.print("parts"+" "+i+" "+partitionList.get(i)+"\n");
		 }
		 
		 if(check==false) {
			 for(int i=y;i<representatives.size();i++) {
				 representatives.get(i).removeAll(representatives.get(i));
			 }
			 for(int i=y+1;i<partitionList.size();i++) {
				 partitionList.get(i).removeAll(partitionList.get(i));
			 }
		 }
		 
		 for(int i=y;i<representatives.size();i++) {
			 pWriter.print("sildi reps"+" "+i+" "+representatives.get(i)+"\n");
		 }
		 
		 for(int i=y+1;i<partitionList.size();i++) {
			 pWriter.print("sildi parts"+" "+i+" "+partitionList.get(i)+"\n");
		 }
		 
	 }
	 
	 public static int LValue(ArrayList<Integer> partEx, int degree) {
		 return (sum(partEx,(degree))-(degree)); //In Grund, the numeration starts with 1. Since we start with 0, it should be -(degree+1)+1, so just -degree
	 }
		
	 public static ArrayList<Permutation> cycleTranspositions(int index, ArrayList<Integer> partition) {
		 int total=sum(partition);
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 int lValue = LValue(partition,index);
		 //System.out.println("partitions and L"+" "+partition+" "+index+" "+lValue);
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
	 public static boolean blockDemo(int index, int r, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition,PrintWriter pWriter) {
		 pWriter.print("partitions"+" "+index+" "+partition+" "+newPartition);
		 for(int i=0;i<representatives.size();i++) {
			 pWriter.print("representatives"+" "+index+" "+representatives.get(i)+"\n");
		 }
		 pWriter.print(Arrays.deepToString(A)+"\n");		 
		 int y = findY(r);
		 boolean check= true;
		 pWriter.print("index r and y"+" "+index+" "+r+" "+y+"\n");
		 int total = sum(partition);
		 ArrayList<Permutation> cycleTrans= cycleTranspositions(index,partition);
		 for(Permutation cycle: cycleTrans) {
			 pWriter.print("cycle"+" "+cycle.toCycleString()+"\n");
		 }
		 for(Permutation cycle: cycleTrans) {
			 boolean formerCheck = formerPermutationsCheck(index,y, total, A, partition, newPartition, cycle,false,pWriter);
			 if(!formerCheck) {
				 pWriter.print("check is false"+"\n");
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
	 
	 public static ArrayList<Permutation> formerPermutations(boolean formerBlocks,int index, int y, int total, PrintWriter pWriter){
		 if(formerBlocks) {
			 return formerPermutationsFromFormerBlocks(y);
		 }else {
			 return formerPermutationsInBlock(index,y,total, pWriter);
		 }
	 }
	 
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
	 
	 public static ArrayList<Permutation> formerPermutationsInBlock(int index, int y, int total, PrintWriter pWriter) {
		 pWriter.print("index y total"+" "+index+" "+y+" "+total+"\n");
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
	 
	 /**
	  * To add a permutation to representatives.
	  * @param index representatives table index
	  * @param perm  Permutation
	  */
	 
	 public static void addRepresentatives(int index, Permutation perm, PrintWriter pWriter) {
		 pWriter.print("add rep"+" "+index+" "+perm.toCycleString()+" "+representatives.size());
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
	 
	 public static Permutation getCanonicalPermutatiom(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation former, Permutation cycle, PrintWriter pWriter) {
		 Permutation cycleM= former.multiply(cycle);
		 pWriter.print("getCanonicalPermutation"+"\n");
		 pWriter.print("permutation"+" "+cycleM+"\n");
		 Permutation canonical = idPermutation(total);
		 //if(descBlockwiseCheck2(index,y,A,newPartition)) {
		 if(descBlockwiseCheck2(index,y,A,newPartition)) {
			 if(!equalBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total), pWriter)) {
				 canonical=getEqualPermutation(cycleM, index, y, total,A, partition, newPartition, pWriter);
				 pWriter.print("equal perm"+" "+canonical.toCycleString()+"\n");
				 canonical=cycle.multiply(canonical);
				 /**if(canonical.isIdentity()) {
					 if(!descBlockCheck2(cycleM,newPartition,index,y,A,idPermutation(total),2,pWriter)) {
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
		 pWriter.print("getCanonical permutation"+" "+canonical.toCycleString());
		 return canonical;
	 }
	 
	 /**
	  * This function performs the canonical test with the former representatives.
	  * That can be the reps from the former representatives in the block. 
	  * If it outputs the canonical permutation as id, then not
	  * passed the canonical test and remove already added new entries from the
	  * list of representatives.
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
	 
	 public static boolean formerPermutationsCheck(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation perm, boolean formerBlocks, PrintWriter pWriter ) {
		 boolean check=true;
		 pWriter.print("index y and formerBlocks"+" "+index+" "+y+" "+formerBlocks+"\n");
		 ArrayList<Permutation> formerPerms= formerPermutations( formerBlocks,index, y, total, pWriter);
		 for(Permutation former: formerPerms) {
			 pWriter.print("formerPERMS"+" "+former+"\n");
		 }
		 Permutation canonical = idPermutation(total);
		 for(Permutation former: formerPerms) {
			 pWriter.print("perm nedir"+" "+perm.toCycleString()+"\n");
			 Permutation test=former.multiply(perm);
			 pWriter.print("cycle multiply by the former"+" "+test+"\n");
			 //if(getCanonicalTest(index, y, total, A, partition, newPartition, test, pWriter)) {
			 canonical=getCanonicalPermutatiom(index, y, total, A, partition, newPartition, former,perm, pWriter);
			 pWriter.print("getCanonical in former perms check"+" "+index+" "+canonical+" "+test+"\n");
			 pWriter.print(canonical.isIdentity()+" "+test.equals(canonical)+" "+descBlockwiseCheck2(index,y,A,newPartition)+" "+Arrays.toString(A)+" "+index+" "+y+" "+partition+" "+"\n");
			 if(canonical.isIdentity()) {
				 //if(descBlockwiseCheck2(index,y,A,newPartition) && descBlockCheck2(test,newPartition,index,y,A,idPermutation(total),2,pWriter)) { //TODO: Not enought for the example. Not des but id so if the array is desc but the comparison is not. Then we cant say continue.
				 if(descBlockwiseCheck2(index,y,A,newPartition)) {	 
					 addRepresentatives(index,idPermutation(total), pWriter);
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
					 addRepresentatives(index, idPermutation(total), pWriter); 
					 /**boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
					 	check=false;
						break;
				 	 }else {
						addRepresentatives(index, idPermutation(total), pWriter); 
				  	 }**/
				 }else {
					 addRepresentatives(index, canonical,pWriter);
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
				 
				 /**if(test.equals(canonical)) { 
					 boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
						 check=false;
						 break;
					 }else {
						 addRepresentatives(index,idPermutation(total)); 
					 //}
				 }else if(!canonical.isIdentity()) {
					 canonical=test.multiply(canonical);
					 addRepresentatives(index, perm.multiply(canonical));
					 boolean formerTest = formerBlocksRepresentatives(index,y, A, newPartition, canonical, pWriter);
					 if(!formerTest) {
						 check=false;
						 break;
					 }else {
						 addRepresentatives(index, canonical);
					 }
				 }else if(canonical.isIdentity() && descBlockwiseCheck2(index,y,A,partition)){
					 addRepresentatives(index,idPermutation(total));
				 }else {
					 check=false;
					 break;
				 }
			 }else {
				 check=false;
				 break;
			 }**/
		 }
		 //System.out.println("former permutation check"+" "+check);
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
	 
	 public static boolean formerBlocksRepresentatives(int index, int y, int[][] A, ArrayList<Integer> newPartition, Permutation perm, PrintWriter pWriter) {
		 boolean check=true;
		 ArrayList<ArrayList<Permutation>> removeList = describeRemovalList(y);
		 /**
		  * If we break within  the nested loop, need to put a label to 
		  * the outer loop also to break.
		  */
		 pWriter.print("formerBlockRepresentatives");
		 pWriter.print("y"+" "+y);
		 outer: for(int i=(y-1);i>=0;i--) {	
			 pWriter.print("former representations line i"+" "+i+"\n");
			 pWriter.print("former representations: "+"\n");
			 for(Permutation p: representatives.get(i)) {
				 pWriter.print(p.toCycleString()+"\n");
			 }
			 
			 /**
			  * No need to multiply them all since we test line by line
			  * in the list of representatives.
			  */
			 
			 for(Permutation rep: representatives.get(i)) {
				 pWriter.print("former in "+" "+i+" "+rep.toCycleString()+"\n");
				 pWriter.print("perm"+" "+perm.toCycleString()+"\n");
				 Permutation test= rep.multiply(perm);
				 if(!descBlockCheck2(test,newPartition,index,y,A,test,2,pWriter)) {
					 check=false;
					 break outer;
				 }else {
					 if(!equalBlockCheck(test,newPartition,index,y,A,test,pWriter)) {
						 if(!rep.isIdentity()) {
							 removeList.get(i).add(rep); 
						 }
					 }
				 }
			 }
		 }
		 pWriter.print("silmeden once");
		 for(int i=0;i<y;i++) {
			 pWriter.print(representatives.get(i)+"\n");
		 }
		 for(int i=0;i<removeList.size();i++) {
			 pWriter.print(removeList.get(i)+"\n");
		 }
		 if(check) {
			 removeFormerPermutations(y,removeList);
		 }
		 return check; 
	 }
	 
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
	 
	 public static boolean equalBlockCheck(Permutation cycleM,ArrayList<Integer> partition, int index, int y, int[][] A, Permutation perm, PrintWriter pWriter){
		 pWriter.print("equalblockcheck"+"\n");
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] original= A[index];
		 pWriter.print(cycleM.toCycleString()+" "+perm.toCycleString()+"\n");
		 Permutation mult=cycleM.multiply(perm);
		 original=actArray(A[mult.get(index)],mult);
		 pWriter.print("equalBlockCheck"+" "+Arrays.toString(canonical)+"\n");
		 pWriter.print("equalBlockCheck"+" "+Arrays.toString(original)+"\n");
		 //System.out.println("equal check"+" "+Arrays.toString(canonical)+" "+Arrays.toString(original));
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
	 public static boolean descBlockCheck2(Permutation cycleM,ArrayList<Integer> partition, int index, int y, int[][] A, Permutation perm,int order,PrintWriter pWriter){
		 boolean check=true;
		 int[] canonical= A[index];
		 int[] original= A[index];
		 if(order==1) {
			 //original=actArray(A[cycleM.get(index)],perm);
			 original=actArray(A[cycleM.get(index)],cycleM);
		 }else if(order==2) {
			 //original=actArray(A[cycleM.get(index)],cycleM.multiply(perm));
			 original=actArray(A[cycleM.get(index)],cycleM);
		 }
		 
		 int[] temp= new int[sum(partition)];
		 if(cycleM.get(index)<index) {
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
	 
	 public static boolean descBlockwiseCheck2(int index, int y, int[][] A, ArrayList<Integer> partition){
		 int[] array= A[index];
		 //System.out.println(Arrays.toString(array));
		 //System.out.println(partition);
		 boolean check=true;
		 int i=0;
		 for(Integer p:partition) {
			 //System.out.println(p);
			 Integer[] can= getBlocks(array,i,p+i);
			 //System.out.println(Arrays.toString(can));
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
	 
	 public static Permutation getCanonicalPermutation(int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, Permutation cycleM, PrintWriter pWriter) {
		 pWriter.print("getCanonicalPermutation"+"\n");
		 pWriter.print("permutation"+" "+cycleM+"\n");
		 //PermutationGroup group= getYoungGroup(newPartition,total);
		 Permutation canonical = idPermutation(total);
		 pWriter.print("des"+" "+descBlockCheck2(cycleM,newPartition,index,y,A,idPermutation(total),2,pWriter)+" "+"\n");
		 pWriter.print("equal"+" "+equalBlockCheck(cycleM,newPartition,index,y,A,idPermutation(total), pWriter)+" "+"\n");
		 if(descBlockCheck2(cycleM,newPartition,index,y,A,idPermutation(total),2,pWriter)){
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
				 canonical=cycleM.multiply(canonical);
			 }
		 }else {
			 
			 /**
			  * des already check equal or descending.
			  */
			 
			 /**canonical=getEqualPermutation(cycleM, index, y, total, A, partition, newPartition, pWriter);
			 if(!canonical.isIdentity()) {
				 canonical=cycleM.multiply(canonical);
			 }**/
		 }
		 return canonical;
	 }
	 
	 public static Permutation getEqualPermutation(Permutation cycleM, int index, int y, int total, int[][] A, ArrayList<Integer> partition, ArrayList<Integer> newPartition, PrintWriter pWriter) {
		 pWriter.print("getEqualPermutation"+"\n");
		 PermutationGroup group= getYoungGroup(newPartition,total);
		 Permutation canonical = idPermutation(total);
		 for(Permutation perm : group.all()) {
			 pWriter.print("equalda check"+" "+perm.toCycleString()+"\n");
			 if(!perm.isIdentity()) {
				 if(equalBlockCheck(cycleM,newPartition,index,y,A,perm,pWriter)){
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
			
			canonicalBlockbasedGenerator(degrees,partition);
			System.out.println(hoppa);
		 }
}
