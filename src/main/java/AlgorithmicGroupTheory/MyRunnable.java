package AlgorithmicGroupTheory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

public class MyRunnable implements Runnable {
	private final int[] degree;

	MyRunnable(int[] degree2) {
		this.degree = degree2;
	}

	public ThreadLocal<ArrayList<Integer>> getPartition(int[] degrees, ArrayList<Integer> partition){
		ThreadLocal<ArrayList<Integer>> initial= new ThreadLocal<ArrayList<Integer>>();
		ArrayList<Integer> newPartition = new ArrayList<Integer>();
		int i=0;
   	 	for(Integer p:partition) {
   	 		Integer[] subArray= trialAbstractClass.getBlocks(degrees,i,p+i);
   	 		newPartition.addAll(trialAbstractClass.getSubPartition(subArray));
   	 		i=i+p; 
   	 	}
   	    initial.set(newPartition);
   	 	return initial;
	}
	
	@Override
	public void run() {
			System.out.println(Arrays.toString(degree));
			ThreadLocal<List<ArrayList<Integer>>> partList= new ThreadLocal<List<ArrayList<Integer>>>();
			//trialAbstractClass.partitionList.clear();
			ThreadLocal<ThreadLocal<List<ArrayList<Integer>>>> formers= new ThreadLocal<ThreadLocal<List<ArrayList<Integer>>>>();
			//trialAbstractClass.formerPermutations.clear();
			trialAbstractClass.initialPartition=trialAbstractClass.getPartition(degree,trialAbstractClass.occurrences);
			//trialAbstractClass.partitionList.add(0,trialAbstractClass.initialPartition);
			List<ArrayList<Integer>> list= new ArrayList<ArrayList<Integer>>();
			list.add(0,trialAbstractClass.initialPartition);
			partList.set(list);
			try {
				trialAbstractClass.genStrip(degree);
			} catch (IOException | CloneNotSupportedException | CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	}
}