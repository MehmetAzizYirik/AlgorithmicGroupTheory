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
 * [1] Kerber, A., Laue, R., Meringer, M., RÃ¼cker, C. and Schymanski, E., 2013.
 * Mathematical chemistry and chemoinformatics: structure generation, elucidation
 * and quantitative structure-property relationships. Walter de Gruyter.
 * 
 * @author Mehmet Aziz Yirik
 */

package AlgorithmicGroupTheory;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.AtomContainerDiscretePartitionRefiner;
import org.openscience.cdk.group.Partition;
import org.openscience.cdk.group.PartitionRefinement;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

public class PermutationGroupFunctions {
	public static PermutationGroup group=PermutationGroup.makeSymN(4);
	public static IChemObjectBuilder builder =  SilentChemObjectBuilder.getInstance();
	public static int order=(int)group.order();
	public static int size=4;
	public static Map<String, Integer> valences; 
	public static boolean verbose = false;
	static String filedir = null;
	static String molecularFormula= null;
	static BufferedWriter fileWriter;
	static {
		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
			
		valences.put("C", 4);
		valences.put("N", 5);
		valences.put("O", 2);
		valences.put("S", 6);
		valences.put("P", 5);
		valences.put("F", 1);
		valences.put("I", 7);
		valences.put("Cl", 5);
		valences.put("Br", 5);
		valences.put("H", 1);
	}
	
	/**
	 * These are the basic functions used for different functions. 
	 */
	
	/**
	 * Calculates the power of a integer. 
	 * @param a a integer
	 * @param b power of the integer 
	 * @return the power the integer
	 */
	
	public static int power(int a, int b) {
		return (int)Math.pow(a, b);
	}
	
	/**
	 * Summation of the connected bond orders.
	 */
	
	public static int orderSummation(IAtomContainer molecule, int index){
		int orderCount=0;
		for(IBond bond: molecule.getAtom(index).bonds()){	
			orderCount=orderCount+bond.getOrder().numeric();
	    }
		return orderCount;
	}
	
	public static int[] buildZerosArray(int n) {
		int arr[] = new int[n];
		for(int i=0;i<arr.length;i++) {
		    arr[i] = 0;
		}
		return arr;
	}
	
	public static int[] cloneArray(int[] arr) {
		int[] result= new int[arr.length];
		for(int i=0;i<arr.length;i++) {
			result[i]=arr[i];
		}
		return result;
	}
	
	public static int[] normalizeArray(int[] arr) {
		for(int i=1;i<arr.length;i++) {
			arr[i]=arr[i]/i;
		}
		return arr;
	}
	
	public static void removeZeroArray(ArrayList<int[]> list) {
		if(zeroArray(list.get(0))){
			list.remove(list.get(0));
		}
	}
	public static ArrayList<Integer> disjoint2List(ArrayList<ArrayList<Integer>> list){
		ArrayList<Integer> result = new ArrayList<Integer>();
		for(ArrayList<Integer> r:list) {
			for(Integer i:r) {
				if(!result.contains(i)) {
					result.add(i);
				}
			}
		}
		return result;
	}
	
	public static void removeLastZeros(HashMap<Integer,Integer> map, int n) {
		for(int i=n;i>0;i--) {
			if(map.get(i)==0) {
				map.remove(i);
			}
			else {
				break;
			}
		}
	}
	
	public static void printListArray(ArrayList<int[]> list) {
		for(int i=0;i<list.size();i++) {
			System.out.println(Arrays.toString(list.get(i)));
		}
	}
	
	public static int findIncrement(int[] list) {
		int result=0;
		int index= notZero(list);
		for(int i=index;i<list.length;i++) {
			if(list[i+1]>list[i]) {
				result=result+i;
				break;
			}
		}
		return result;
	}
	public static boolean allEqual(int[] arr) {
		int index= notZero(arr);
		boolean equal= true;
		for(int i=index;i<arr.length-1;i++) {
			if(arr[i+1]!=arr[i]) {
				equal=false;
				break;
			}
		}
		return equal;
	}
	
	public static int notZero(int[] arr) {
		int index=0;
		for(int i=0;i<arr.length;i++) {
			if(arr[i]!=0) {
				index= index+i;
				break;
			}
		}
		return index;
	}
	
	public static boolean zeroArray(int[] bond) {
		boolean check= true;
		for(int i=0;i<bond.length;i++) {
			if(bond[i]!=0) {
				check=false;
				break;
			}
		}
		return check;
	}
	
	public static ArrayList<int[]> inc= new ArrayList<int[]>();
	public static boolean contains(ArrayList<int[]> list,int[] arr) {
		boolean check= false;
		for(int i=0;i<list.size();i++) {
			if(Arrays.equals(list.get(i), arr)) {
				check= true;
			}
		}
		return check;
	}
	
	public static void checkAdd(ArrayList<int[]> list, int[] arr) {
		if(!contains(list,arr)) {
			list.add(arr);
		}
	}
	
	public static void checkAddList(ArrayList<int[]> list, ArrayList<int[]> list2) {
		for(int[] l: list2) {
			checkAdd(list,l);
		}
	}
	
	public static int sum(ArrayList<Integer> list) {
		int sum=0;
		for(Integer i:list) {
			sum=sum+i;
		}
		return sum;
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
	
	public static int count(ArrayList<Integer> list, int i) {
		int count=0;
		for(Integer k:list) {
			if(k.equals(i)) {
				count++;
			}
		}
		return count;
	}
	
	/**
	 * Counting the frequency of integers in the lists.
	 * @param list	 the list of integer sets
	 * @return map of the frequencies.
	 */
	
	public static HashMap<Integer,Integer> countEntries(ArrayList<ArrayList<Integer>> list) {
		HashMap<Integer,Integer> map= new HashMap<Integer,Integer>();
		for(Integer s: getBase(group)) {
			map.put(s,0);
		}
		for(ArrayList<Integer> l:list) {
			int ind1= l.get(0);
			int ind2= l.get(1);
			map.put(ind1, map.get(ind1)+1);
			map.put(ind2, map.get(ind2)+1);
		}
		return map;
	}
	
	public static final Comparator<Integer> ASC_ORDER = new Comparator<Integer>() {
	    public int compare(Integer e1, Integer e2) { 
	        return -e2.compareTo(e1);
	    }
	};
	
	
	public static void arrAsc(int[] a) {
		int temp;
		for (int i = 0; i < a.length; i++) 
        {
            for (int j = i + 1; j < a.length; j++) 
            {
                if (a[i] > a[j]) 
                {
                    temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                }
            }
        }
	}
	public static final Comparator<ArrayList<Integer>> ASC_ORDER2 = new Comparator<ArrayList<Integer>>() {
		
		public int compare(ArrayList<Integer> l1, ArrayList<Integer> l2) { 
	    	return l1.toString().compareTo(l2.toString());
	    }
	};
	
	public static  Comparator<Integer> DES_ORDER = new Comparator<Integer>() {
	    public int compare(Integer e1, Integer e2) { 
	        return e2.compareTo(e1);
	    }
	};
	
	public static Comparator<String> DES_STRING_ORDER= new Comparator<String>() {
		public int compare(String o1, String o2) {
	        return -Integer.valueOf(valences.get(o1)).compareTo(Integer.valueOf(valences.get(o2)));
	    }
	};
	
	/**
	 * (Done)
	 * @param line
	 * @return
	 */
	public static boolean desOrderCheck(int[] line) {
		boolean check=true;
		int length=line.length;
		for(int i=0;i<length-1;i++) {
			if(line[i]<line[i+1]) {
				check=false;
				break;
			}
		}
		return check;
	}
	
	public static int[] addElement(int[] a, int e) {
        a  = Arrays.copyOf(a, a.length + 1);
        a[a.length - 1] = e;
        return a;
    }
	
	public static ArrayList<Integer> addElement2(ArrayList<Integer> a, int e) {
        a.add(e);
        return a;
    }
	
	public static int sum(int[] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	public static int sum(int[] array, int index) {
		int sum=0;
		for(int i=0;i<=index;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	/**
	 * SUm function for Gale Ryser Theorem. The right side of inequality. 
	 * @param array integer array
	 * @param k integer to check min value
	 * @return sum of the array with the minimality check. 
	 */
	public static int sumGale(int[] array, int k) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+Math.min(array[i], k);
		}
		return sum;
	}
	public static List<ArrayList<Integer>> buildArray2(int n,int d, int depth){
		List<ArrayList<Integer>> array= new ArrayList<ArrayList<Integer>>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(ArrayList<Integer> item: partition2(n-i,d,depth+1)) {
					item=addElement2(item,i);
			        if(item.size()==d) {
			        	if(sum(item)==n) {
			        		array.add(item);
			        	}
			        }else {
			        	array.add(item);
			        }
			}
		}
		return array;
	}
	
	public static List<int[]> buildArray(int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(int[] item: partition(n-i,d,depth+1)) {
					item=addElement(item,i);
			        if(item.length==d) {
			        	if(sum(item)==3) { 
			        		array.add(item);
			        	}
			        }else {
			        	array.add(item);
			        }
			}
		}
		return array;
	}
	
	/**
	 * Generating simple list. Just to put 1 into the ith index.
	 * @param var length of the list
	 * @param index where the entry is 1
	 * @return simple list with ones and zeros.
	 */
	
	public static ArrayList<Integer> simpleList(int var, int index){
		ArrayList<Integer> list= new ArrayList<Integer>();
		for(int i=0;i<var;i++) {
			if(i==index) {
				list.add(1);
			}else {
				list.add(0);
			}
		}
		return list;
	}
	
	public static int[] multiplyArray(int[] arr, int a) {
		int[] res= new int[arr.length];
		for(int i=0;i<arr.length;i++) {
			res[i]=arr[i]*a;
		}
		return res;
	}
	
	public static ArrayList<Integer> multiplyList(ArrayList<Integer> arr, int a) {
		ArrayList<Integer> list= new ArrayList<Integer>();
		for(int i=0;i<arr.size();i++) {
			list.add(arr.get(i)*a);
		}
		return list;
	}
	
	public static int[] sumArray(int[] arr1, int[] arr2) {
		int[] arr3= new int[arr1.length];
		for(int i=0;i<arr1.length;i++) {
			arr3[i]=arr1[i]+arr2[i];
		}
		return arr3;
	}
	
	public static ArrayList<Integer> sumList(ArrayList<Integer> arr1, ArrayList<Integer> arr2) {
		ArrayList<Integer> arr3= new ArrayList<Integer>();
		for(int i=0;i<arr1.size();i++) {
			arr3.add(arr1.get(i)+arr2.get(i));
		}
		return arr3;
	}
	
	public static HashMap<ArrayList<Integer>,Integer> multMultinomial(HashMap<int[],Integer> mult1,HashMap<int[],Integer>mult2) {
		HashMap<ArrayList<Integer>,Integer> result= new HashMap<ArrayList<Integer>,Integer>();
		for(int[] key:mult1.keySet()) {
			for(int[] key2:mult2.keySet()) {
				ArrayList<Integer> sumlist=arrayToList(sumArray(key,key2));
				if(result.containsKey(sumlist)) {
					result.put(sumlist, result.get(sumlist)+(mult1.get(key)*mult2.get(key2)));
				}else {
					result.put(arrayToList(sumArray(key,key2)),mult1.get(key)*mult2.get(key2));
				}
			}
		}
		return result;
	}
	
	public static HashMap<ArrayList<Integer>,Integer> multMultinomialList(HashMap<ArrayList<Integer>,Integer> mult1,HashMap<ArrayList<Integer>,Integer>mult2) {
		HashMap<ArrayList<Integer>,Integer> result= new HashMap<ArrayList<Integer>,Integer>();
		for(ArrayList<Integer> key:mult1.keySet()) {
			for(ArrayList<Integer> key2:mult2.keySet()) {
				ArrayList<Integer> sumlist=sumList(key,key2);
				if(result.containsKey(sumlist)) {
					result.put(sumlist, result.get(sumlist)+(mult1.get(key)*mult2.get(key2)));
				}else {
					result.put(sumlist,mult1.get(key)*mult2.get(key2));
				}
			}
		}
		return result;
	}
	
	public static ArrayList<Integer> arrayToList(int[] array) {
		ArrayList<Integer> list= new ArrayList<Integer>();
		for(int i=0;i<array.length;i++) {
			list.add(array[i]);
		}
		return list;
	}
	public static HashMap<ArrayList<Integer>,Integer> idMultinomial(int m,int classSize){
		HashMap<ArrayList<Integer>,Integer> mult= new HashMap<ArrayList<Integer>,Integer>();
		ArrayList<Integer> list= new ArrayList<Integer>();
		for(int i=0;i<m;i++) {
			list.add(0);
		}
		mult.put(list, classSize);
		return mult;
	}
	public static void divideByOrder(HashMap<ArrayList<Integer>,Integer> map) {
		for(ArrayList<Integer> key:map.keySet()) {
			map.replace(key, map.get(key), (Integer)(map.get(key)/order));
		}
	}
	
	public static HashMap<ArrayList<Integer>,Integer> sumLists(HashMap<ArrayList<Integer>,Integer> map1, HashMap<ArrayList<Integer>,Integer> map2){
		Set<ArrayList<Integer>> keySet = new HashSet<ArrayList<Integer>>();
        keySet.addAll(map1.keySet());
        keySet.addAll(map2.keySet());

        HashMap<ArrayList<Integer>, Integer> map3 = new HashMap<ArrayList<Integer>, Integer>();
        Integer val1, val2;
        for (ArrayList<Integer> key : keySet) {
            val1 = map1.get(key);
            val1 = (val1 == null ? 0 : val1);
            val2 = map2.get(key);
            val2 = (val2 == null ? 0 : val2);
            map3.put(key, val1 + val2);
        }
		return map3;
	}
	
	/**
	 * ***********************************************************************
	 */
	
	/**
	 * Acts a permutation p on a given set of integers. 
	 * The permutation permutates the integer set by the re-arrangement. 
	 * 
	 * @param set the set of integers 
	 * @param p a permutation from a permutation group
	 * @return modified set based on the acts of the permutation
	 */
	
	public static ArrayList<Integer> act(ArrayList<Integer> set, Permutation p) {
		ArrayList<Integer> modifiedSet = new ArrayList<Integer>();
		for(Integer element:set) {
			modifiedSet.add(p.get(element));
		}
		return modifiedSet;
	}
	
	public static int[] actArray(int[] strip, Permutation p) {
		int length= strip.length;
		int[] modified = new int[length];
		for(int i=0; i<length;i++) {
			modified[p.get(i)]=strip[i]; 
		}
		return modified;
	}
	
	public static Set<Integer> act(Set<Integer> set, Permutation p) {
		Set<Integer> modifiedSet = new HashSet<Integer>();
		for(Integer element:set) {
			modifiedSet.add(p.get(element-1)+1);
		}
		return modifiedSet;
	}
	
	
	
	/**
	 * Acts a group  on a given set of integers.  
	 * 
	 * @param set the set of integers 
	 * @param group a group acting on the set 
	 * @return modified set based on the acts of the group
	 */
	
	public static ArrayList<ArrayList<Integer>> actGroup(ArrayList<Integer> set, PermutationGroup group) {
		ArrayList<ArrayList<Integer>> list= new ArrayList<ArrayList<Integer>>();
		for(Permutation perm: group.all()) {
			ArrayList<Integer> modified = act(set, perm);
			if(!list.contains(modified)) {
				list.add(modified);
			}
		}
		return list;
	}
	
	/**
	 * Acts a product group  on a given set of integers.  
	 * 
	 * @param set the set of integers 
	 * @param productGroup a product group acting on the set 
	 * @return modified set based on the acts of the product group
	 */
	//TODO: This should be an action on a bijection. m over ( n 2) is the set from MOLGEN page 42.
	//TODO: This m maps are mapping from vertex pair to edge multiplicity. 
	public static ArrayList<ArrayList<Integer>> actProductGroup(ArrayList<Integer> set, ArrayList<ArrayList<Permutation>> productGroup) {
		ArrayList<ArrayList<Integer>> list= new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Permutation> perm: productGroup) {
			
		}
		return list;
	}
	
	/**
	 * All the mappings describing all the m-multigraphs on n nodes
	 * @param m maximum multiplicity in m-multigraph
	 * @param n number of nodes
	 */
	public static void mapping(int m, int n) {
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		
	}
	/**
	 * Group action for direct product groups
	 * @param product Direct product group
	 * @param group a permutation group
	 * @return action of direct product group on a permutation group
	 */
	
	public static ArrayList<Permutation> groupActionProduct(ArrayList<ArrayList<Permutation>> product, PermutationGroup group) {
		ArrayList<Permutation> result= new ArrayList<Permutation>();
		for(Permutation perm: group.all()) {
			for(ArrayList<Permutation> list: product) {
				result.add(list.get(0).multiply(perm.multiply(list.get(1).invert())));
			}
		}
		return result;
	}
	

	/**
	 * Group action of direct product groups on vertex pairs.
	 * @param product Direct product group
	 * @param group a permutation group
	 * @return action of direct product group on a permutation group
	 */
	
	public static ArrayList<Permutation> groupActionProduct(ArrayList<ArrayList<Permutation>> product, int k) {
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(k,getBase(group));
		ArrayList<Permutation> result= new ArrayList<Permutation>();
		for(Permutation perm: group.all()) {
			for(ArrayList<Permutation> list: product) {
				result.add(list.get(0).multiply(perm.multiply(list.get(1).invert())));
			}
		}
		return result;
	}
	/**
	 * Gets the size of group's base
	 * @param group a permutation group
	 */
	
	public static void setSize(PermutationGroup group) {
		PermutationGroupFunctions.size= group.getSize();
	}
    
	/**
	 * Generates the subsets of a integer set. 
	 * @param set
	 * @return list of unique subsets
	 */
	
	public static Set<Set<Integer>> subSet( Set<Integer> set ) {
		Integer[] elements = set.toArray(new Integer [set.size()]);
        int power = 1 << elements.length;
        Set<Set<Integer>> subsets = new HashSet<Set<Integer>>();
        for( int count = 0; count < power; count++ ) {
            Set<Integer> subset = new HashSet<Integer>();
            for( int bit = 0; bit < elements.length; bit++ ) {
            	if( (count & 1 << bit) != 0 ) {
                    subset.add( elements[bit] );
                }
            }
            subsets.add( subset );
        }
        return subsets;
    }
	
	/**
	 * Generates all the k-subsets whose size is k
	 * @param k size of subset
	 * @param set a integer set
	 * @return list of k-subsets
	 */
	
	public static Set<Set<Integer>> ksubSet(int k,ArrayList<Integer> set ) {
		Integer[] elements = set.toArray(new Integer [set.size()]);
        int power = 1 << elements.length;
        Set<Set<Integer>> subsets = new HashSet<Set<Integer>>();
        for( int count = 0; count < power; count++ ) {
            Set<Integer> subset = new HashSet<Integer>();
            for( int bit = 0; bit < elements.length; bit++ ) {
                if( (count & 1 << bit) != 0 ) {
                    subset.add( elements[bit] );
                }
            }
            if(subset.size()==k) {
                subsets.add( subset );
            }
        }
        return subsets;
    }
	
	/**
	 * It is the same as ksubSet function. The input and output types are ArrayList rather
	 * than Set<Integer>
	 * 
	 * Generates all the k-subsets whose size is k
	 * @param k size of subset
	 * @param set a integer set
	 * @return list of k-subsets
	 */
	
	public static ArrayList<ArrayList<Integer>> ksubSet2(int k,ArrayList<Integer> set ) {
		Integer[] elements = set.toArray(new Integer [set.size()]);
        int power = 1 << elements.length; 
        ArrayList<ArrayList<Integer>> subsets = new ArrayList<ArrayList<Integer>>();
        for( int count = 0; count < power; count++ ) {
        	ArrayList<Integer> subset = new ArrayList<Integer>();
            for( int bit = 0; bit < elements.length; bit++ ) {
                if( (count & 1 << bit) != 0 ) {
                    subset.add( elements[bit] );
                }
            }
            if(subset.size()==k) {
                subsets.add( subset );
            }
        }
        return subsets;
    }
	
	/**
	 * Checks the given set and the list of subsets are disjoint or not.
	 * @param list list of subsets
	 * @param set the set to check with the list
	 * @return boolean disjoint or not.
	 */
	
	public static boolean disjointCheck(ArrayList<ArrayList<Integer>> list, ArrayList<Integer> set) {
		boolean check=true;
		for(ArrayList<Integer> l:list) {
			if(!Collections.disjoint(l,set)) {
				check=false;
				break;
			}
		}
		return check;
	}
	
	/**
	 * Generating the k-disjoint subsets.
	 * The disjoint k subsets are needed for Young subgroup representation.
	 * 
	 * @param k the subset size
	 * @param set integer set 
	 * @return the list of disjoint k-subsets 
	 */
	
	public static ArrayList<ArrayList<ArrayList<Integer>>> disjoints= new ArrayList<ArrayList<ArrayList<Integer>>>();
	public static ArrayList<ArrayList<ArrayList<Integer>>> disjoint(int k,ArrayList<Integer> set ){
		ArrayList<ArrayList<Integer>> list= ksubSet2(k,set);
		ArrayList<ArrayList<Integer>> firstSets= new ArrayList<ArrayList<Integer>>(); // first element is zero
		for(ArrayList<Integer> l: list) {
			if(l.get(0)==0) {
				firstSets.add(l);
			}
		}
		list.removeAll(firstSets);
		for(ArrayList<Integer> s:firstSets) {
			ArrayList<ArrayList<Integer>> dis= new ArrayList<ArrayList<Integer>>();
			dis.add(s);
			disjointBuild(dis,list);
		}
		return disjoints;
	}
	
	/**
	 * Cloning the input list, ArrayList<ArrayList<Integer>>
	 * @param list the list of subsets
	 * @return clone of the list
	 */
	
	public static ArrayList<ArrayList<Integer>> cloneList(ArrayList<ArrayList<Integer>> list){
		ArrayList<ArrayList<Integer>> result= new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> i:list) {
			result.add(i);
		}
		return result;
	}
	
	/**
	 * Checks two subset lists are disjoint or not. 
	 * @param list1 the list of subsets
	 * @param list2 the other list of subsets
	 * @return boolean disjoint or not.
	 */
	
	public static boolean disjointCheckAll(ArrayList<ArrayList<Integer>> list1, ArrayList<ArrayList<Integer>>list2) {
		boolean check= false;
		for(ArrayList<Integer> l2: list2) {
			if(disjointCheck(list1,l2)) {
				check=true;
				break;
			}
		}
		return check;
	}
	
	/**
	 * For the list of subsets, generating the disjoint subsets.
	 * In the disjoint function, the subsets having 0 are listed.
	 * Then these elements are used for the disjoint comparison and
	 * generation of disjoint subsets for each of these subsets with 0.
	 * 
	 * @param set
	 * @param list
	 */
	
	public static void disjointBuild(ArrayList<ArrayList<Integer>> set, ArrayList<ArrayList<Integer>> list ) {
		for(ArrayList<Integer> l:list) {
			if(disjointCheck(set,l)) {
				ArrayList<ArrayList<Integer>> set2=cloneList(set);
				set2.add(l);
				if(disjointCheckAll(set2,list)) {
					disjointBuild(set2,list);
				}else {
					disjoints.add(set2);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param group a permutation group
	 * @param set a set of integers
	 * @return the list of permutations stabilizing the input set
	 */
	
	public static ArrayList<Permutation> getStabilizers(PermutationGroup group,ArrayList<Integer> set){
		ArrayList<Permutation> stabilizers = new ArrayList<Permutation>();
		for(Permutation p:group.all()){
			ArrayList<Integer> stabilizer = new ArrayList<Integer>();
			for(Integer i :set) {
				stabilizer.add(p.get(i));
			}
			if(stabilizer.equals(set)) {
				stabilizers.add(p);
			}
		}
		return stabilizers;
	}
	
	/**
	 * Rather than the group action on the integer set, the action
	 * of a permutation is considered.
	 * 
	 * @param g a permutation
	 * @param sets list of subsets
	 * @return the list of stabilized subsets under the permutation action.
	 */
	
	public static int getStabilizers2(Permutation g, ArrayList<ArrayList<Integer>> sets){
		ArrayList<ArrayList<Integer>> stabilizedSets = new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> set: sets) {
			ArrayList<Integer> stabilize = new ArrayList<Integer>();
			for(Integer i :set) {
				stabilize.add(g.get(i));
			}
			if(stabilize.equals(set)) {
				stabilizedSets.add(set);
			}
		}
		return stabilizedSets.size();
	}
	
	/**
	 * Checks the integer set is stable or not under the permutation action
	 * 
	 * @param g a permutation
	 * @param set a set of integers
	 * @return 1 or 0 meaning the set is stable or not.
	 */
	
	public static int checkStab(Permutation g, ArrayList<Integer> set){
		int no=0;
        if(set.equals(act(set,g))) {
			no++;
		}
		return no;
	}
	
	/**
	 * Returns the cycle type of a permutation. In the permutation,
	 * cycles are counted based on their lengths and the counting return
	 * as the cycle type array of the permutation.
	 * 
	 * ( The algorithm is from C.A.G.E.S. Book. )
	 * 
	 * @param n the size of the return array (group size +1)
	 * @param g
	 * @return cycle type array of the permutation.
	 */
	
	public static int[] type(int n,Permutation g) {
		boolean[] p = new boolean[size];
		int[] t = new int[n];
		for(int i=0;i<n-1;i++) {
			p[i]=true;
			t[i+1]=0;
		}
		for(int i=0;i<n-1;i++) {
			if(p[i]) {
				int l=1;
				int j=i;
				p[j]=false;
				while(p[g.get(j)]) {
					l=l+1;
					j=g.get(j);
					p[j]=false;
				}
				t[l]=t[l]+1;
			}
		}
		return t;
	}
	
	/**
	 * Counting the n combinatoion of m.
	 * @param m integer
	 * @param n integer
	 * @return the combination value of two integers.
	 */
	
	public static int combination(int m, int n) {
		return factorial(m) / (factorial(n) * factorial(m - n));
	}
	
	/**
	 * Calculating factorial of an integer.
	 * @param i integer
	 * @return factorial of the integer.
	 */
	
	public static int factorial(int i){
		if (i==0){
			return 1;
		}
		return i * factorial(i - 1);
	}
	
	/**
	 * Calculating the size of a Youngsubgroup. 
	 * If a young subgroup is the direct product of its subgroups build
	 * with respect to the integer partition, then the multiplication of 
	 * the factorial of each entry in the partition gives the Young group
	 * order.
	 * @param partition ArrayList<Integer> partition to build a Young subgroup.
	 * @return int multiplication of the factorial of each entry of the partition.
	 */
	
	public static int youngGroupOrder(ArrayList<Integer> partition){
		int result=1;
		for(Integer i:partition) {
			if(i!=1) {
				result*=factorial(i);
			}
		}
		return result;
	}
	
	/**
	 * Building the orbits of a integer set under the group action of
	 * a permutation group
	 * 
	 * @param group a permutation group
	 * @param set a integer set
	 * @return the list of orbits under the group action.
	 */
	
	public static ArrayList<ArrayList<Integer>> getOrbits(PermutationGroup group, ArrayList<Integer> set){
		ArrayList<ArrayList<Integer>> orbits = new ArrayList<ArrayList<Integer>>();
		for(Permutation p:group.all()){
			ArrayList<Integer> orbit = new ArrayList<Integer>();
			for(Integer i :set) {
				orbit.add(p.get(i));
			}
			orbit.sort(ASC_ORDER);

			if(!orbits.contains(orbit)) {
				orbits.add(orbit);
			}
		}
		orbits.sort(ASC_ORDER2);
		return orbits;
	}
	
	/**
	 * Orbits of a group by the group action.
	 * @param group a permutation group
	 * @param group2 a permutation group
	 * @return the list of orbits. 
	 */
	
	public static ArrayList<ArrayList<Permutation>> getOrbits(PermutationGroup group, PermutationGroup group2){
		ArrayList<ArrayList<Permutation>> orbits = new ArrayList<ArrayList<Permutation>>();
		for(Permutation p:group.all()){
			ArrayList<Permutation> orbit = new ArrayList<Permutation>();
			for(Permutation p2 :group2.all()) {
				orbit.add(p.multiply(p2));
			}
			if(!orbits.contains(orbit)) {
				orbits.add(orbit);
			}
		}
		return orbits;
	}
	
	/**
	 * From the list of orbits, getting the minimum ones from the lexicographical order.
	 * @param group a permutation group 
	 * @param set an integer set
	 * @return set of minimal orbits. 
	 */
	
	public static ArrayList<ArrayList<Integer>> getMinOrbits(PermutationGroup group, ArrayList<Integer> set){
		ArrayList<ArrayList<Integer>> orbits = new ArrayList<ArrayList<Integer>>();
		for(Permutation p:group.all()){
			ArrayList<Integer> orbit = new ArrayList<Integer>();
			for(Integer i :set) {
				orbit.add(p.get(i));
			}
			orbit.sort(ASC_ORDER);
			if(!orbits.contains(orbit)) {
				orbits.add(orbit);
			}
		}
		return orbits;
	}
	
	/**
	 * Gets list of permutation under the orbits.
	 * @param group a permutation group
	 * @param set a integer set
	 * @return the list of permutations.
	 */
	
	public static ArrayList<ArrayList<Permutation>> getOrbitsPerm(PermutationGroup group, ArrayList<Permutation> set){
		ArrayList<ArrayList<Permutation>> orbits = new ArrayList<ArrayList<Permutation>>();
		for(Permutation p:group.all()){
			ArrayList<Permutation> orbit = new ArrayList<Permutation>();
			for(Permutation i :set) {
				orbit.add(p.multiply(i));
			}
			if(!orbits.contains(orbit)) {
				orbits.add(orbit);
			}
		}
		return orbits;
	}
	

	public static ArrayList<ArrayList<Permutation>> getOrbitsPermList(PermutationGroup group, ArrayList<ArrayList<Permutation>> set){
		ArrayList<ArrayList<Permutation>> orbits = new ArrayList<ArrayList<Permutation>>();
		for(ArrayList<Permutation> s:set) {
			orbits.addAll(getOrbitsPerm(group,s));
		}
		return orbits;
	}
	
	/**
	 * Gets all the orbits without putting them in ascending order and without checking duplicates. 
	 * @param group a permutation group
	 * @param set a integer set
	 * @return all the orbits under the group action
	 */
	
	public static ArrayList<ArrayList<Integer>> getOrbits2(PermutationGroup group, ArrayList<ArrayList<Integer>> set){
		ArrayList<ArrayList<Integer>> orbits = new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> s:set) {
			orbits.addAll(getOrbits(group,s));
		}
		return orbits;
	}
	
	/**
	 * Gets the base of a permutation group.
	 * @param group a permutation group
	 * @return the base set of the group
	 */
	
	public static ArrayList<Integer> getBase(PermutationGroup group) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		for(int i=0;i<group.getSize();i++) {
			base.add(i);
		}
		return base;
	}
	
	/**
	 * 
	 * Cauchy Frobenious Theorem - Number of Orbits
	 *
	 **/
	
	/**
	 * Counts the number of orbits for k-subsets.
	 * @param k size of subsets
	 * @return number of orbits
	 */
	
	public static int numberOfOrbits(int k) {
		int count=0;
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(k,getBase(group));
		for(Permutation perm:group.all()) {
			count+=getStabilizers2(perm,subsets);
		}
		return count/group.all().size();
	}
	
	public static int numberOfOrbit(PermutationGroup group, ArrayList<ArrayList<Integer>> subset) {
		int count=0;
		for(Permutation perm:group.all()) {
			for(ArrayList<Integer> set:subset) {
				count=count+checkStab(perm,set);
			}
		}
	return count/group.all().size();
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
	
	public static void depictHighlight(IAtomContainer molecule,Iterable<IChemObject> object, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		//depiction.withSize(200, 200).withZoom(4).depict(molecule).writeTo(path);
		depiction.withHighlight(object, Color.GREEN).withCarbonSymbols().withAtomValues().withAtomMapNumbers().withAtomColors().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
		//depiction.withAtomColors().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path);
	}
	
	public static HashMap<Integer, Integer> setCycleList(int size){
		HashMap<Integer, Integer> list = new HashMap<Integer, Integer>(size);
		for(int i=1;i<=size;i++) {
			list.put(i, 0);
		}
		return list;
	}
	
	/**
	 * Builds a map of cycle types in the permutation group
	 * @param group
	 * @return map of cycle types in the group
	 */
	
	public static HashMap<Permutation, HashMap<Integer, Integer>> cycleType(PermutationGroup group) {
		HashMap<Permutation, HashMap<Integer, Integer>> map= new HashMap<Permutation, HashMap<Integer, Integer>>();
		for(Permutation perm: group.all()) {
			HashMap<Integer, Integer> list = setCycleList(group.getSize());
			String [] cycles = perm.toCycleString().split("\\)"); 
			for(String l:cycles) { 
				int count=0;
				String[] length= l.split(",");
				for(String d:length) {
					count=count+d.length()-1;
				}
				list.replace(count, list.get(count), (list.get(count)+1));
			}
			map.put(perm, list);
		}
		return map;
	}
	
	/**
	 * Building the conjugacy classes of the permutation group.
	 * @param group a permutation group
	 * @return
	 */
	
	public static Multimap<HashMap<Integer, Integer>, Permutation> conjugacyClasses(PermutationGroup group) {
		Multimap<HashMap<Integer, Integer>, Permutation> multiMap = ArrayListMultimap.create();
		HashMap<Permutation, HashMap<Integer, Integer>> map=cycleType(group);
		for(Permutation entry:map.keySet()) {
			multiMap.put(map.get(entry),entry);
		}
		return multiMap;
	}
	
	/**
	 * Generating a group from a permutation
	 * @param perm a permutation	
	 * @return permutation group generated by a permutation
	 */
	
	public static PermutationGroup generateGroup(Permutation perm) {
		List<Permutation> generators = new ArrayList<Permutation>();
        generators.add(perm);
        return new PermutationGroup(size, generators);
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
	 * Grund 1.1.16
	 * Get generators for a set of elements to build their permutation group.
	 * Due to the Permutation class input should be int[].
	 * @return
	 */
	
	public static List<Permutation> getGenerators(int[] values,int total){
		List<Permutation> generators = new ArrayList<Permutation>();
		int length = values.length;
		if(length==1) {
			generators.add(new Permutation(total));
		}else if (length==2) {
			generators.add(new Permutation(values));
		}else {
			int[] firstValues= new int[] {values[0],values[1]};
			Permutation perm1= new Permutation(firstValues);
			Permutation perm2= new Permutation(values);
			generators.add(perm1);
			generators.add(perm2);
		}
		return generators;
	}
	
	/**
	 * This is the Younggroup for a given partition.
	 * Here, without generators, we just go through all the permutations in the main group.
	 * We check for each permutation, whether they keep the partition the same or modify it.
	 * getYounggroupL indicates list. So output is a List<Permutation>
	 * 
	 * @param partition occurences of the atoms as the partition.
	 * @param total number of atoms
	 * @return
	 */
	
	public static List<Permutation> getYoungGroupL(ArrayList<Integer> partition,int total){
		PermutationGroup group  = PermutationGroup.makeSymN(total);
		ArrayList<ArrayList<Integer>> parts=subPartitionsList(partition);
		List<Permutation> perms= new ArrayList<Permutation>();
		for(Permutation perm:group.all()) {
			if(actOnPartition(parts,perm)) {
				perms.add(perm);
			}
		}
		return perms;
	}
	
	public static boolean actOnPartition(ArrayList<ArrayList<Integer>> subParts, Permutation p){
		ArrayList<Integer> mod= new ArrayList<Integer>();
		boolean check=true;
		for(ArrayList<Integer> part: subParts) {
			mod=act(part,p);
			mod.sort(DES_ORDER);
			if(!listEqual(part,mod)) {
				check=false;
			}
		}
		return check;
	}
	
	public static ArrayList<Integer> listToArr(int[] arr){
		ArrayList<Integer> list= new ArrayList<Integer>();
		for(int i=0;i<arr.length;i++) {
			list.add(arr[i]);
		}
		list.sort(DES_ORDER);
		return list;
	}
	
	public static boolean listEqual(ArrayList<Integer> list, ArrayList<Integer> list2) {
		boolean check=true;
		for(int i=0;i<list.size();i++) {
			if(list.get(i)!=list2.get(i)) {
				check=false;
				break;
			}
		}
		return check;
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
	 * For 3.3.3, rather than generating a Youngsubgroup, we take the representatives.
	 * These representatives are used in the canonical test not all the group members.
	 * @param reps Representative permutations built by the previous partition
	 * @param cycles Cycles are generated based on the current partition and the previous one.
	 * @return
	 */
	
	public static void updatePermTree(ArrayList<Permutation> cycles){
		ArrayList<Permutation> newReps = new ArrayList<Permutation>();	
		for(Permutation cycle: cycles) {
			for(Permutation rep: permTree) {
				Permutation n=rep.multiply(cycle);
				if(!newReps.contains(n)) {
					newReps.add(n);
				}
			}
		}
		permTree=newReps;
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
	 * Counts the number of graphs by the acting of the permutation group 
	 * @param group a permutation group
	 * @param m multiplicity
	 * @return number of graphs
	 */
	
	public static int numberOfGraphs(PermutationGroup group, int m) {
		int count=0;
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		Multimap<HashMap<Integer, Integer>,Permutation> classes= conjugacyClasses(group);
		for(HashMap<Integer, Integer> key:classes.keySet()) {
			Permutation perm = (Permutation) classes.get(key).toArray()[0];
			System.out.println(perm.toCycleString());
			PermutationGroup group2=generateGroup(perm);
			count=count+classes.get(key).size()*power(m,numberOfOrbit(group2,subsets));
		}
		return count/group.all().size();
	}
	
	/**
	 * Polya's enumeration method:
	 * 
	 * Group reduction method  for which cycle index polynomial is built.
	 */
	
	/**
	 * Building a map of cycle index coefficients for a permutation group. 
	 * Cylce index polynomial relies on the the cycle type function.
	 * @param group
	 * @return cycle index coefficients of the conjugacy classes representatives
	 */
	
	public static HashMap<String, Integer> cycleIndexCoefficient(PermutationGroup group) {
		Multimap<String, String> map = ArrayListMultimap.create();
		HashMap<String, Integer> coef = new HashMap<String, Integer>();
		for(Permutation perm: group.all()) {
			map.put(Arrays.toString(type((group.getSize()+1),perm)),Arrays.toString(type((group.getSize()+1),perm)));
		}
		for(String key:map.keySet()) {
			coef.put(key, map.get(key).size());
		}
		return coef;
	}
	
	/**
	 * Set of orbits of length i of <g> on X
	 * 
	 * @param group permutation group
	 * @param set a integer set
	 * @param l length of orbits
	 * @return the set of orbits wth length l
	 */
	
	public static List<Set<Integer>> getLengthBasedOrbits(PermutationGroup group, Set<Set<Integer>> set, int l){
		List<Set<Integer>> orbits = new ArrayList<Set<Integer>>();
		for(Permutation p:group.all()){
			for(Set<Integer> s:set) {
				Set<Integer> orbit = new HashSet<Integer>();
				for(Integer i :s) {
					orbit.add(p.get(i));
				}
				if(orbit.size()==l && !orbits.contains(orbit)) {
					orbits.add(orbit);
				}
			}
		}
		return orbits;
	}
	
	/**
	 * Numbe rof orbits with the specified orbit length.
	 * @param perm a permutation
	 * @param set set of subsets
	 * @param i orbits length
	 * @return the number of specified orbits.
	 */
	
	public static int numberOfLengthOrbits(Permutation perm, Set<Set<Integer>> set,int i) {
		PermutationGroup group= generateGroup(perm);
		return getLengthBasedOrbits(group,set,i).size();
	}
	
	/**
	 * Building the array of cycle type for a permutation. But the induced
	 * form of the permutation is considered. For more details;
	 * please look at MOLGEN Book. 
	 * 
	 * @param perm a permutation
	 * @return cycle type array
	 */
	
	public static int[] cycleTypeInduced(Permutation perm) {
		ArrayList<ArrayList<Integer>> subset= ksubSet2(2,getBase(group));
		int size = group.getSize();
		int[] cycle= buildZerosArray(size+1);
        List<Permutation> generators = new ArrayList<Permutation>();
        generators.add(perm);
        PermutationGroup group2 = new PermutationGroup(size, generators);
		for(ArrayList<Integer> set:subset) {
			cycle[getOrbits(group2,set).size()]=cycle[getOrbits(group2,set).size()]+1;
		}
		return normalizeArray(cycle);
	}
	
	/**
	 * Permutation Types like explain in MOLGEN Book.
	 * Ordering the cycle types in descending order. 
	 */
	
	public static ArrayList<Integer> cyclicFactors(Permutation perm) {
		int[] cyclic=type(group.getSize()+1,perm);
		ArrayList<Integer> lengths= new ArrayList<Integer>();
		for(int i=1;i<cyclic.length;i++) {
			if(cyclic[i]>0) {
				lengths.add(cyclic[i]);
			}
		}
		Collections.sort(lengths,DES_ORDER);
		return lengths;
	}
	/**
	 * Like given in MOLGEN Book, the n-tuple representation of cycles
	 * @param perm a permutation
	 * @return map of ntuples of a permutation
	 */
	
	public static HashMap<Integer,Integer> nTuple(Permutation perm) {
		ArrayList<Integer> list= cyclicFactors(perm);
		Collections.sort(list, DES_ORDER);
		System.out.println("cyclic"+" "+list);
		int n=sum(list);
		HashMap<Integer,Integer> map= new HashMap<Integer,Integer>();
		for(int i=1;i<=n;i++) {
			map.put(i,count(list,i));
		}
		removeLastZeros(map,n);
		return map;
	}

	/**
	 * Partitioning integer n into d parts. 
	 * @param n integer to partition
	 * @param d number of parts o partition
	 * @param depth
	 * @return The list of integer partition 
	 */
	
	public static List<ArrayList<Integer>> partition2(int n, int d,int depth) {
		if(d==depth) {
			List<ArrayList<Integer>> array= new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer> take=new ArrayList<Integer>();
			array.add(take);
			return array;
		}
		return buildArray2(n,d,depth);
		
	}
	
	public static List<int[]> partition(int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray(n,d,depth);
		
	}
	
	/**
	 * Calculating the multinomial coefficient for the given power array
	 * 
	 * @param n main power over the multinomial.
	 * @param array power array, the powers of each element in the multinomial
	 * @return specified moltinomial coefficient
	 */
	
	public static int multiCoef(int n, int[] array) {
		int mult=1;
		for(int i=0;i<array.length;i++) {
			mult=mult*factorial(array[i]);
		}
		return factorial(n)/mult;
	}
	
	public static int multiCoefList(int n, ArrayList<Integer> array) {
		int mult=1;
		for(int i=0;i<array.size();i++) {
			mult=mult*factorial(array.get(i));
		}
		return factorial(n)/mult;
	}
	
	/**
	 * Building a monomial for number of variable with the power.
	 * @param var number of variables
	 * @param inPower the power of the variables
	 * @return map representing the multinomials. 
	 */
	
	public static HashMap<ArrayList<Integer>,Integer> monomial(int var, int inPower){
		HashMap<ArrayList<Integer>,Integer> map= new HashMap<ArrayList<Integer>,Integer>();
		for(int i=0;i<var;i++) {
			map.put(multiplyList(simpleList(var,i),inPower), 1);
		}
		return map;
	}
	
	/**
	 * Construction of a multinomial for the given number of variable and power information 
	 * @param power power of the multinomial
	 * @param var number of variables
	 * @param inPower power of each variable
	 * @return a multinomial.
	 */
	
	public static HashMap<ArrayList<Integer>,Integer> multinomialList(int power, int var , int inPower ) {
		if(power==1) {
			return monomial(var,inPower);
		}
		HashMap<ArrayList<Integer>,Integer> map= new HashMap<ArrayList<Integer>, Integer>();
		for(ArrayList<Integer> arr:partition2(power,var,0)) {
			map.put(multiplyList(arr,inPower),multiCoefList(power,arr));
		}
		return map;
	}
	
	public static HashMap<int[],Integer> multinomial(int power, int var , int inPower ) {
		HashMap<int[],Integer> map= new HashMap<int[], Integer>();
		for(int[] arr:partition(power,var,0)) {
			map.put(multiplyArray(arr,inPower),multiCoef(power,arr));
		}
		return map;
	}
	
	
	/**
	 * Group Reduction Function. The formula given in MOLGEN Book pg 77.
	 * @param group permutation group 
	 * @param m multiplicity
	 * @return value of GRF for the multiplicity and permutation group
	 */
	
	public static HashMap<ArrayList<Integer>, Integer> GRF(PermutationGroup group, int m) {
		Multimap<HashMap<Integer, Integer>, Permutation> map=conjugacyClasses(group);
		HashMap<ArrayList<Integer>, Integer> sum= new HashMap<ArrayList<Integer>, Integer>();
		for(HashMap<Integer, Integer> key:map.keySet()) {
			Permutation perm = (Permutation) map.get(key).toArray()[0];
			HashMap<ArrayList<Integer>,Integer> mult= idMultinomial(m,map.get(key).size());
			for(int i=1;i<cycleTypeInduced(perm).length;i++) {
				if(cycleTypeInduced(perm)[i]!=0) {
					int power= cycleTypeInduced(perm)[i];
					mult=multMultinomialList(mult,multinomialList(power,m,i));
				}
			}
			sum=sumLists(sum,mult);
		}
		divideByOrder(sum);
		return sum;
	}
	
	
	
	/**
	 * Coset and Double Cosets on Groups - MOLGEN 1.36
	 *
	 */
	
	public static ArrayList<Permutation> leftCoset(Permutation permutation, PermutationGroup group) {
		ArrayList<Permutation> list= new ArrayList<Permutation>();
		for(Permutation perm: group.all()) {
			list.add(permutation.multiply(perm));
		}
		return list;
	}
	
	public static ArrayList<Permutation> rightCoset(Permutation permutation, PermutationGroup group) {
		ArrayList<Permutation> list= new ArrayList<Permutation>();
		for(Permutation perm: group.all()) {
			list.add(perm.multiply(permutation));
		}
		return list;
	}
	
	public static HashMap<Permutation, ArrayList<Permutation>> leftCosets(PermutationGroup group, PermutationGroup subGroup) {
		HashMap<Permutation, ArrayList<Permutation>> map= new HashMap<Permutation, ArrayList<Permutation>>(); 
		for(Permutation perm: group.all()) {
			ArrayList<Permutation> list= new ArrayList<Permutation>();
			for(Permutation perm2: subGroup.all()) {
				list.add(perm.multiply(perm2));
			}
			map.put(perm, list);
		}
		return map;
	}
	
	public static HashMap<Permutation, ArrayList<Permutation>> rightCosets(PermutationGroup group, PermutationGroup subGroup) {
		HashMap<Permutation, ArrayList<Permutation>> map= new HashMap<Permutation, ArrayList<Permutation>>(); 
		for(Permutation perm: group.all()) {
			ArrayList<Permutation> list= new ArrayList<Permutation>();
			for(Permutation perm2: subGroup.all()) {
				list.add(perm2.multiply(perm));
			}
			map.put(perm, list);
		}
		return map;
	}
	
	/**
	 * Direct product of two permutation groups is also a group.
	 * 
	 * @param group  a permutation group
	 * @param group2 another permutation group
	 * @return direct prodcut group
	 */
	
	public static ArrayList<ArrayList<Permutation>> groupDirectProduct(PermutationGroup group, PermutationGroup group2) {
		ArrayList<ArrayList<Permutation>> product = new ArrayList<ArrayList<Permutation>>(); 
		for(Permutation perm: group.all()) {
			for(Permutation perm2: group2.all()){
				ArrayList<Permutation> list= new ArrayList<Permutation>();
				list.add(perm);
				list.add(perm2);
				product.add(list);
			}
		}
		return product;
	}
	
	
	
	
	/**
	 * Direct product of a group with a list of stabilizers. The stabilizers set is also a permutation group.
	 * @param group a permutation group
	 * @param stab a stabilizier group
	 * @return direct product of permutation group with the stabilizer group.
	 */
	
	public static ArrayList<ArrayList<Permutation>> groupDirectProduct2(PermutationGroup group, ArrayList<Permutation> stab) {
		ArrayList<ArrayList<Permutation>> product = new ArrayList<ArrayList<Permutation>>(); 
		for(Permutation perm: group.all()) {
			for(Permutation perm2: stab){
				ArrayList<Permutation> list= new ArrayList<Permutation>();
				list.add(perm);
				list.add(perm2);
				product.add(list);
			}
		}
		return product;
	}
	
	
	/**
	 * 
	 * Double Cosets: Orbits for the group build by the group action of a product group
	 *
	 * @param product direct product group
	 * @param group a permutation group
	 * @return lis of double cosets of the permutation group
	 */
	
	public static ArrayList<ArrayList<Permutation>> doubleCosets(ArrayList<ArrayList<Permutation>> product, PermutationGroup group) {
		ArrayList<ArrayList<Permutation>> result= new ArrayList<ArrayList<Permutation>>();
		for(Permutation perm: group.all()) {
			result.add(doubleCoset(product,perm));
		}
		return result;
	}
	
	/**
	 * Double coset of a permutation
	 * 
	 * @param product direct product group
	 * @param perm a permutation
	 * @return double cosets of a permutation
	 */
	
	public static ArrayList<Permutation> doubleCoset(ArrayList<ArrayList<Permutation>> product, Permutation perm) {
		ArrayList<Permutation> result= new ArrayList<Permutation>();
		for(ArrayList<Permutation> list: product) {
			result.add(list.get(0).multiply(perm.multiply(list.get(1).invert())));
		}
		return result;
	}
	
	/**
	 * Double coset of an integer. For gluing lemma
	 * @param product direct product group
	 * @param x an integer
	 * @return double cosets of a integer
	 */
	public static ArrayList<Permutation> doubleCosetInt(ArrayList<ArrayList<Permutation>> product, Permutation perm) {
		ArrayList<Permutation> result= new ArrayList<Permutation>();
		for(ArrayList<Permutation> list: product) {
			result.add(list.get(0).multiply(perm.multiply(list.get(1).invert())));
		}
		return result;
	}
	
	
	/**
	 * Fundamental Lemma: orbits, cosets and double cosets
	 */
	
	/**
	 * G/Gx , the orbits of the stabilizer,Gx, by the group action of G.
	 * @param group permutation group
	 * @param subGroup permutation subgroup
	 * @param set integer set
	 * @return the list of permutations in the orbit of the stabilizer.
	 */
	
	public static ArrayList<Permutation> stabilizerOrbits(PermutationGroup group, PermutationGroup subGroup, ArrayList<Integer> set) {
		ArrayList<Permutation> stabilizers= getStabilizers(group,set);
		ArrayList<Permutation> result= new ArrayList<Permutation>();
		for(Permutation perm: subGroup.all()) {
			for(Permutation perm2: group.all()) {
				for(Permutation stab: stabilizers) {
					result.add(perm2.multiply(perm.multiply(stab)));
				}
			}
		}
		return result;
	}
	
	/**
	 * G/Gx , the orbits of the stabilizer,Gx, by the group action of G.
	 * @param group 	a permutation group
	 * @param subGroup 	a permutation subgroup
	 * @param set 		list of subsets.
	 * @return the list of permutations in the orbits of the stabilizer
	 */
	
	public static ArrayList<ArrayList<Permutation>> stabilizerOrbitsList(PermutationGroup group, PermutationGroup subGroup, ArrayList<ArrayList<Integer>> sets) {
		ArrayList<ArrayList<Permutation>> result= new ArrayList<ArrayList<Permutation>>();
		for(ArrayList<Integer> set: sets) {
			System.out.println("set"+" "+set);
			result.add(stabilizerOrbits(group,subGroup,set));
		}
		return result;
	}
	
	
	/**
	public static void constructTrans(PermutationGroup group, PermutationGroup subGroup) {
		ArrayList<ArrayList<Integer>> subSets= ksubSet2(2,getBase(group));
		ArrayList<ArrayList<ArrayList<Permutation>>> lis= new ArrayList<ArrayList<ArrayList<Permutation>>>();
		ArrayList<ArrayList<ArrayList<Permutation>>> lis2= new ArrayList<ArrayList<ArrayList<Permutation>>>();
		ArrayList<ArrayList<Permutation>> stabOrbits = stabilizerOrbitsList(group,subGroup,subSets);
		for(ArrayList<Integer> set: subSets) {
			ArrayList<ArrayList<Permutation>> dCosets    = doubleCosets(groupDirectProduct2(subGroup, getStabilizers(group,set)),group);
			for(ArrayList<Permutation> l: dCosets) {
				System.out.println("dcoset"+" "+dCosets.size()+" "+l.size());
			}
		}
		for(ArrayList<Permutation> list: stabOrbits) {
			System.out.println("stab"+" "+list.size()+" "+stabOrbits.size());
		}
	}**/
	

	/**
	 * MOLGEN - minimal orbit representation
	 * @param group a permutation group 
	 * @param subSets the list of subsets
	 * @return the minimal orbit representative of the subsets
	 */
	
	public static ArrayList<Integer> minRep(PermutationGroup group, ArrayList<ArrayList<Integer>> subSets) {
		ArrayList<Integer> min= new ArrayList<Integer>();
		for(ArrayList<Integer> set: subSets) {
			ArrayList<ArrayList<Integer>> orbit=getOrbits(group, set);
			if(orbit.get(0).toString().compareTo(min.toString())<0) {
				min=orbit.get(0);
			}
		}
		return min;
	}
	
	/**
	 * MOLGEN - Algorithm 5.1: Bond adding function P
	 */
	
	/**
	 * Returns the minimal source index for the bond adding.
	 * @param list list of vertex labels.
	 * @return source index for bond adding.
	 */
	
	public static int minSourceIndex(ArrayList<ArrayList<Integer>> list) {
		ArrayList<Integer> atoms = getBase(group);
		HashMap<Integer,Integer> map= countEntries(list);
		for(Integer key:map.keySet()) {
			if(map.get(key)>=2) {
				atoms.remove(key);
			}
		}
		atoms.sort(ASC_ORDER);
		return atoms.get(0);
	}
	
	/**
	 * Returns the minimal target index for the bond adding.
	 * @param list list of vertex labels.
	 * @return target index for bond adding.
	 */
	
	public static int minTargetIndex(ArrayList<ArrayList<Integer>> list) {
		ArrayList<Integer> atoms = getBase(group);
		for(ArrayList<Integer> l: list) {
			atoms.removeAll(l);
		}
		atoms.sort(ASC_ORDER);
		return atoms.get(0);
	}
	
	/**
	 * Extends the bond list my adding the minimal bond to the list.
	 * @param list the list of vertex labels.
	 * @return adding a new vertex pair to the list meaning a new bond.
	 */
	
	public static ArrayList<ArrayList<Integer>> minBondAdd(ArrayList<ArrayList<Integer>> list) {
		ArrayList<Integer> pair= new ArrayList<Integer>();
		pair.add(minSourceIndex(list));
		pair.add(minTargetIndex(list));
		list.add(pair);
		return list;
	}
	
	
	/**
	 * MOLGEN: Homomorphism principle
	 * 
	 * Mapping m multigraph to (m-1) multigraph
	 */
	
	/**
	 * The algorithm is from molgen book. 
	 * @param pair vertex pair 
	 * @return
	 */
	public static int bond(ArrayList<Integer> pair) {
		return 1;
	}
	
	/**
	 * Returns an empty array in the size of multiplicity. 
	 * @param n multiplicity
	 * @return an n-length empty array
	 */
	
	public static int[] bondMultiplicity(int n) {
		ArrayList<ArrayList<Integer>> sets=ksubSet2(2,getBase(group));
		int[] mult= new int[sets.size()];
		for(int i=0;i<sets.size();i++) {
			mult[i]=0;
		}
		return mult;
	}
	
	/**
	 * Homomorphism Principle
	 * @param pair vertex pair 
	 * @param m bond multiplicity
	 * @return  mapping the m-multigraph to (m-1)
	 */
	
	public static int homomorphism(ArrayList<Integer> pair, int m) {
		if(bond(pair)<(m-1)) {
			return bond(pair);
		}else {
			return (m-2); 
		}
	}
	
	/**
	 * From molgen book, the inverse homomorphism for m-multiplicity. 
	 */
	
	/**
	 * 
	 * Increases the bond multiplcity and used in the inverse homomorphism function. 
	 * 
	 * @param bond the array of bond information.
	 * @param index index of the bond in the bond array
	 */
	
	public static void bondIncrease(int[] bond, int index) {
		if(index!=0) {
			for(int i=0;i<=(index-1);i++) {
					bond[bond.length-1-i]=bond[bond.length-1-i]+1;
			}	
		}
	}
	
	/**
	 * As given in MOLGEN book, incrementing the bond multiplicity by one
	 * in the bond array.
	 * @param bond the array of bond multiplicities
	 * @return the modified bond array.
	 */
	
	public static ArrayList<int[]> inverseHomomorphism(int[] bond) {
		ArrayList<int[]> result= new ArrayList<int[]>();
		int size=bond.length+1; 
		for(int i=0;i<size;i++) {
			result.add(cloneArray(bond));
		}
		for(int i=0;i<size;i++) {
			bondIncrease(result.get(i),i);
		}
		return result;
	}
	
	
	/** 
	 * Besides the inverHomomorphism function, another way of 
	 * increasing the bond multiplicities in m-multigraphs. 
	 * @param bond bond multiplicity array
	 * @param m maximum multiplicity in m-multigraph
	 */
	
	public static void increaseBond(int[] bond, int m) {
		if(zeroArray(bond)) {
			for(int[] arr:increaseToM(bond)) {
				increaseBond(arr,m);
			}
		}else if(allEqual(bond) && bond[0]!=m ) {
			increaseToM(bond,m);
		}else {
			int index= findIncrement(bond);
			int zero= notZero(bond);
			checkAdd(inc, cloneArray(bond));
			//inc.add(cloneArray(bond));
			for(int i=index;i>=zero;i--) {
				if(bond[i]!=0) {
					bond[i]=bond[i]+1;
				}
				checkAdd(inc, cloneArray(bond));
				//inc.add(cloneArray(bond));
			}
			if(!allEqual(bond)) {
				increaseBond(bond, m);
			}else {
				if(bond[bond.length-1]!=m) {
					increaseToM(bond,m);
				}
			}
		}
	}
	
	/**
	 * Generation of m-multigraphs. For a given set of numbers increasing
	 * each bond multiplicity by one.  
	 * 
	 * @param bond bond multiplicity array
	 * @return modified list of bond arrays
	 */
	
	public static ArrayList<int[]> increaseToM(int[] bond) {
		int iter;
		if(zeroArray(bond)) {
			iter = bond.length+1;
		}else {
			iter = bond.length-notZero(bond)+1;
		}
		ArrayList<int[]> result= new ArrayList<int[]>();
		for(int i=0;i<iter;i++) {
			result.add(cloneArray(bond));
		}
		for(int i=0;i<iter;i++) {
			bondIncrease(result.get(i),i);
		}
		removeZeroArray(result);
		checkAddList(inc, result);
		//inc.addAll(result);
		return result;
	}
	
	/**
	 * Generation of m-multigraphs. For a given set of numbers increasing them 
	 * in an order until they reach the upper limit, m.
	 * @param bond bond multiplicity array
	 * @param m    maximum multiplicity to reach.
	 */
	
	public static void increaseToM(int[] bond,int m) {
		int iter;
		if(zeroArray(bond)) {
			iter = bond.length+1;
		}else {
			iter = bond.length-notZero(bond)+1;
		}
		ArrayList<int[]> result= new ArrayList<int[]>();
		for(int i=0;i<iter;i++) {
			result.add(cloneArray(bond));
		}
		for(int i=0;i<iter;i++) {
			bondIncrease(result.get(i),i);
		}
		
		removeZeroArray(result);
		checkAddList(inc, result);
		//inc.addAll(result);
		if((result.get(result.size()-1)[bond.length-1])!=m) {
			increaseToM(result.get(result.size()-1),m);
		}
	}

	
	/**
	 * Young Subgroup construction
	 */
	
	/**
	 * Build the Young subgroup for the given disjoint set
	 * @param set   a disjoint set
	 * @param group a permutation group
	 * @return Young subgroups for the disjoint set
	 */
	
	public static ArrayList<Permutation> youngSubgroup(ArrayList<ArrayList<Integer>> set,PermutationGroup group) {
		ArrayList<Permutation> permutations = new ArrayList<Permutation>();
		ArrayList<Integer> elements=getBase(group);
		elements.removeAll(disjoint2List(set));
		for(Permutation perm: group.all()) {
			for(Integer i: elements) {
				if(perm.get(i)==i) {
					permutations.add(perm);
				}
			}
		}
		return permutations;
	}
	
	/**
	 * For all the disjoint sets, generating the Young subgroups.
	 * @param sets a set of disjoint subsets
	 * @param group a permutation group
	 * @return Young subgroups of the disjoint sets
	 */
	
	public static ArrayList<ArrayList<Permutation>> youngSubgroups(ArrayList<ArrayList<ArrayList<Integer>>> sets, PermutationGroup group) {
		ArrayList<ArrayList<Permutation>> subGroups=  new ArrayList<ArrayList<Permutation>>();
		for(ArrayList<ArrayList<Integer>> set: sets) {
			subGroups.add(youngSubgroup(set,group));
		}
		return subGroups;
	}
	

	/**
	 * For the given size,k, generating all the k-disjoint subsets first; then
	 * generating the Young subgroups for the Permutation group.
	 * @param k size of subsets
	 * @param group a permutation group
	 * @return list of young subgroups for the k-subsets.
	 */
	
	public static ArrayList<ArrayList<Permutation>> youngSubgroupRepresentation(int k, PermutationGroup group) {
		ArrayList<ArrayList<ArrayList<Integer>>> kDisjointSets= disjoint(k,getBase(group));
		return youngSubgroups(kDisjointSets,group);
	}
	
	/**
	 * Degree Partition- Gruner Thesis Chapter 7
	 */
	
	/**
	 * For all the permutations of symmetry group n; the permutations in lexicographical
	 * order are the set partitions. cdk group package has partition also but not 
	 * generating all the partitions of a set. Still might be useful.  
	 */
	
	//TODO: Think about degree partition vs permutations over the degree set.
	
	/**
	 * Gale-Ryser Criteria is used in the degree partition step.	
	 * @param a int array for degree sequences
	 * @param b int array for degree sequences
	 * @return true if satisfies the criteria.
	 */
	
	public static boolean galeRyserCriteria(int[] a, int[] b) {
		boolean check= true;
		a = Arrays.stream(a).boxed().sorted((k, l) -> l.compareTo(k)).mapToInt(i -> i).toArray();
		if(sum(a)==sum(b)) {
			for(int i=0;i<a.length;i++) {
				if(!(sum(a,i)<=sumGale(b, i))) {
					check=false;
					break;
				}
			}
		}else {
			check= false;
		}
		return check;
	}
	
	
	/**
	 * 
	 * Gluing Lemma: this method of joining partial structures
	 * without producing isomorphic duplicates is known as gluing
	 * lemma. Here, partial structures are the substructures
	 * called macroatoms.the fixed map is a fixed m-multigraph. 
	  * This means the list of vertex pair with bond multiplicity. 
	 * 
	 **/
	
	 //TODO: How should we decide for a fixed m-multigraph ?
	 public static void gluingLemma(ArrayList<Integer> Y) {
		 ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		 PermutationGroup groupX = PermutationGroup.makeSymN(subsets.size());
		 PermutationGroup groupY = PermutationGroup.makeSymN(Y.size());
		 ArrayList<ArrayList<Permutation>> productGroup = groupDirectProduct(groupX, groupY);
		 ArrayList<Integer> multiplicities=groupActionMaps(productGroup, subsets);
	 }
	 
	 /**
	  * For substructures, automorphism group of the graphs
	  * should be calculated. 
	  */
	 
	public static PermutationGroup automorphismGroup(String formula) {
		IAtomContainer ac= MolecularFormulaManipulator.getAtomContainer(formula, builder);
		AtomContainerDiscretePartitionRefiner refiner = PartitionRefinement.forAtoms().create();
	    PermutationGroup autoGroup = refiner.getAutomorphismGroup(ac);
	    return autoGroup;
	 }
	
	/** 
	 * 
	 * Like in dioxine example, MOLGEN, the atom symbols should
	 * be classified. The gluing lemma is the for the bijection between
	 * the symmetry groups of both sets.
	 * 
	 */
	 //TODO: Direct product can be a 2*n matrix 
	 
	/**
	 * Returns the equivalence classes of the atoms.
	 * @param ac atom container
	 * @return equivalence class partition of atoms.
	 */
	
	public static Partition automorphismPartition(IAtomContainer ac){
		AtomContainerDiscretePartitionRefiner refiner = PartitionRefinement.forAtoms().create();
		Partition autoPart = refiner.getAutomorphismPartition(ac);
		return autoPart;
	}
	
	 /**
	  * The group action on set of mappings. Molgen pg 52.
	  * The group action is on the vertex pairs not on the
	  * edge multiplicity of the mappings.
	  * 
	  * The group action is part of gluing lemma.
	  */
	 
	 public static ArrayList<Integer> groupActionMaps( ArrayList<ArrayList<Permutation>> group, ArrayList<ArrayList<Integer>> set) {
		 ArrayList<Integer> list = new ArrayList<Integer>();
		 for(ArrayList<Permutation> perm: group) {
			 for(ArrayList<Integer> s: set) {
				 //list.add(replace(s,perm));
			 }
		 }
		 return list;
	 }	
	 
	 /**
	  * The double coset generation for gluing lemma. From the 1994 first paper. 
	  * @param G permutation group
	  * @param H permutation group
	  */
	 
	 public static ArrayList<Integer> doubleCosetsGluing(ArrayList<ArrayList<Permutation>> group, ArrayList<ArrayList<Integer>> set) {
		 ArrayList<Integer> list = new ArrayList<Integer>();
		 for(ArrayList<Permutation> perm: group) {
			 for(ArrayList<Integer> s: set) {
				 list.add(gluingLeft(s,perm));
			 }
		 }
		 return list;
	 }
	 
	 public static int fixedMap(ArrayList<Integer> pair) {
		 return 1;
	 }
	 
	 
	 /**
	  * The left and the right side of the gluing lemma. 
	  * @param s
	  * @param perm
	  * @return
	  */
	 
	 public static int gluingLeft(ArrayList<Integer> s, ArrayList<Permutation> perm) {
		 ArrayList<Integer> set= new ArrayList<Integer>(); 
		 Permutation inverse = perm.get(0).invert();
		 int multiplicity= perm.get(1).get(fixedMap(act(s,inverse)));
		 return multiplicity;
	 }
	 
	 //TODO: This fixed map cannot map everything to the same vertex pair or the same
	 //TODO: Multiplicity. Otherwise why do we need gluing lemma if the right side always
	 //TODO: Mapped to the same vertex pair for each f elements of maps. Understand the maps.
	 public static int gluingRight(ArrayList<Integer> s, ArrayList<Permutation> perm) {
		 ArrayList<Integer> set= new ArrayList<Integer>(); 
		 Permutation inverse = perm.get(0).invert();
		 int multiplicity= perm.get(1).get(fixedMap(act(s,inverse)));
		 return multiplicity;
	 }
	 
	 /**
	  * To be able to understand the gluing lemma, we need to find a fixed map 
	  * in the set of m-multigraphs on n nodes. In MOLGEN Book, all the maps
	  * are considered as a multigraph. 
	  */
	
	 /**
	  * This is the functions finding the fixed vertex pair under the
	  * group action. Here, the group is generated by a permutation.
	  * The definition of fixed map is given on page 81. 
	  * @param p permutation
	  * @return fixed vertex pair.
	  */
	 
	 public static ArrayList<ArrayList<Integer>> findFixedMap(Permutation p) {
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		ArrayList<ArrayList<Integer>> subsets2= ksubSet2(2,getBase(group));
		for(Permutation perm:generateGroup(p).all()) {
			for(ArrayList<Integer> m:subsets) {
				if(!m.equals(act(m,perm))) {
					subsets2.remove(m);
				}
			}
		}
		return subsets2;
	 }
	 
	 /**
	  * Check atom is saturated or not.
	  */
	 
		
	public static boolean saturationChecker(IAtomContainer molecule, int index) throws CloneNotSupportedException, CDKException, IOException{
		if ((orderSummation(molecule,index))>= (int)valences.get(molecule.getAtom(index).getSymbol())){ 
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * Get the list of atom indices with the free valences.
	 * @param ac IAtomContainer
	 * @return the list of unsaturated atom indices.
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
	 */
	public static SaturationChecker saturation;
	 public static ArrayList<Integer> getFreeValences(IAtomContainer ac) throws CloneNotSupportedException, CDKException, IOException {
		 ArrayList<Integer> freeV= new ArrayList<Integer>();
		 for(int i=0;i<ac.getAtomCount();i++) {
			 if(!saturationChecker(ac,i)) {
				 freeV.add(i);
			 }
		 }
		 return freeV;
	 }
	 
	 /**
	  * For the list of atom symbol, grouping the same symbols
	  * and assigning the acting permutation group. 
	  * @param list
	  * @return map of symbols and acting permutation group
	  */
	 
	 public static HashMap<Character,PermutationGroup> symbolClassification(ArrayList<Character> list) {
		 ArrayList<Character> unique = new ArrayList<Character>();
		 HashMap<Character,PermutationGroup> map= new HashMap<Character,PermutationGroup>();
		 for(Character symbol: list) {
			 if(!unique.contains(symbol)) {
				 unique.add(symbol);
			 }
		 }
		 for(Character sym: unique) {
			 map.put(sym, PermutationGroup.makeSymN(countCharacter(list,sym)));
		 }
		 return map;
	 }
	 
	 /**
	  * Counting the number of occurences of symbols in the character list
	  * @param list list of symbols
	  * @param symbol atom symbol
	  * @return number of occurence
	  */
	 
	 public static int countCharacter(ArrayList<Character> list, Character symbol) {
		 int count=0;
		 for(Character sym: list) {
			 if(sym.equals(symbol)) {
				 count++;
			 }
		 }
		 return count;
	 }
	 
	 
	 //TODO: This direct product group is needed for the group action on symbol list.
	 /**
	  * Generating the direct product group of a list of permutation groups
	  * @param list list of permutation groups
	  */
	 
	 public static void multiDirectProduct(ArrayList<PermutationGroup> list) {
		 ArrayList<ArrayList<Permutation>> grp= groupDirectProduct(list.get(0),list.get(1));
		 for(int i=2;i<list.size();i++) {
			 grp=directProductList(grp,list.get(i));
		 }
	 }
	 
	 /**
	  * 
	  * Direct group product of more than 2 groups.
	  * @param perm list of permutation pairs
	  * @param grp permutation group
	  * @return direct product
	  */
	 public static ArrayList<ArrayList<Permutation>> directProductList(ArrayList<ArrayList<Permutation>> perm, PermutationGroup grp) {
		 for(Permutation prm: grp.all()) {
			 for(ArrayList<Permutation> lst: perm) {
				 lst.add(prm);
			 }
		 }
		 return perm;
	 }
	 
	 /**
	  * Direct Product of List of Permutation Groups
	  * @param groups list of permutation groups
	  * @return the direct produtct of the groups
	  */
	 
	 public static ArrayList<ArrayList<Permutation>> directProductGroups(ArrayList<PermutationGroup> groups) {
		 groups.iterator();
		 ArrayList<ArrayList<Permutation>> product=groupDirectProduct(groups.get(0),groups.get(1)); 
		 for(int i=2;i<groups.size();i++) {
			 product=productListGroup(product,groups.get(i));
		 }
		 return product;
	 }
	 
	 
	 public static ArrayList<ArrayList<Permutation>> productListGroup(ArrayList<ArrayList<Permutation>> list, PermutationGroup group) {
		 ArrayList<ArrayList<Permutation>> out= new ArrayList<ArrayList<Permutation>>();
		 for(Permutation perm: group.all()) {
			 for(ArrayList<Permutation> pair: list) {
				 ArrayList<Permutation> l= new ArrayList<Permutation>();
				 l.addAll(pair);
				 l.add(perm);
				 out.add(l);
			 }	 
		 }
		 return out;
	 }
	 
	 
	 public static String fixedMapDioxine(int index) {
		 String atom="";
		 if(0<=index && index<=3) { //How did we numerated these free valences?
			 atom=atom+"Cl";
		 }else {
			 atom=atom+"H";
		 }
		 return atom;
	 }
	 
	 /**
	  * Group action on fixedMapDioxine.
	  * @param group an acting group
	  * @param index input for the fixed map.
	  * @return map of new indices and the return values.
	  */
	 
	 public static HashMap<Integer, String> actionOnMap(PermutationGroup group, int index) {
		 HashMap<Integer, String> map = new HashMap<Integer, String>();
		 for(Permutation perm: group.all()) {
			 int index2=perm.invert().get(index);
			 map.put(index2, fixedMapDioxine(index2));
		 }
		 return map;
	 }
	 
	 /**
	  * Right group action on string values.
	  * @param group
	  * @param group2
	  * @param index
	  */
	 
	 public static void rightAction(PermutationGroup group, PermutationGroup group2, int index) {
		 HashMap<Integer, String> map= actionOnMap(group, index);		 
		 for(String str: map.values()) {
			 for(Permutation perm: group2.all()) {
				 for(Permutation perm2: group2.all()) {
					 perm.multiply(perm2); // Group action on strings ?	 
				 }
			 } 
		 }
	 }
	 
	 /**
	  * Generating the groups acting on the sublists
	  * @param set list of strings
	  * @return map of sublist and their acting groups
	  */
	 
	 public static HashMap<List<String>, PermutationGroup> actionGroupsOnStrings(ArrayList<String> set) {
		 HashMap<List<String>, PermutationGroup> map= new HashMap<List<String>, PermutationGroup>();
		 ArrayList<List<String>> sublists=subListGenerator(set);
		 for(List<String> list: sublists) {
			 map.put(list, PermutationGroup.makeSymN(list.size()));
		 }
		 return map;
	 }
	 
	 /**
	  * Group action on the sublist of a list of strings
	  * @param set list of strings
	  * @param group 
	  */
	 public static void actionOnStrings(ArrayList<String> set) {
		 HashMap<List<String>, PermutationGroup> map= actionGroupsOnStrings(set);
		 ArrayList<PermutationGroup> groups = new ArrayList<>(map.values());
		 ArrayList<ArrayList<Permutation>>  directProduct=directProductGroups(groups);
		 ArrayList<List<String>> list= subListGenerator(set);
		 for(ArrayList<Permutation> perm: directProduct) {
			 for(List<String> l:list) {
				 
			 }
		 }
		 
	 }
	 
	 /**
	  * Check at which indices the string element is changed. 
	  * @param set list of String elements
	  * @return the list of indices where the next element is changed. 
	  */
	 
	 public static ArrayList<Integer> splitIndices(ArrayList<String> set) {
		 ArrayList<Integer> indices= new ArrayList<Integer>();
		 for(int i=0;i<set.size()-1;i++) {
			 if(set.get(i)!=set.get(i+1)) {
				 indices.add(i+1);
			 }
		 }
		 return indices;
	 }
	 
	 /**
	  * Groups the equivalent strings and generates the sublists.
	  * @param list list of strings
	  * @return list of sublists
	  */
	 
	 public static ArrayList<List<String>> subListGenerator(ArrayList<String> list){
		 Collections.sort(list);
		 ArrayList<List<String>> subsets= new ArrayList<List<String>>();
		 ArrayList<Integer> indices= splitIndices(list);
		 int start=0;
		 for(Integer index: indices) {
			 List<String> set= list.subList(start, index);
			 subsets.add(set);
			 start=index;
		 }
		 List<String> set= list.subList(start, list.size());
		 subsets.add(set);
		 return subsets;
	 }
	 
	 
	 /**
	  * Generating the double cosets in gluing lemma. 
	  * @param group a permutation group
	  * @return
	  */
	 public static String gluingLemmaDoubleCosets(int index){
		 String out="";
		 //TODO: Direct product of more than 2 groups should be described and its action.
		 return out;
		 
	 }
	 
	 public static ArrayList<Permutation> asc_order(Permutation first, Permutation second) {
		 ArrayList<Permutation>  perms= new ArrayList<Permutation>();
		 int size1= first.getValues().length;
		 int size2= second.getValues().length;
		 if(size1<size2) {
			 perms.add(first);
			 perms.add(second);
		 }else{
			 perms.add(second);
			 perms.add(first);
		 }
		 return perms;
	 }
	 
	 /**
	  * Multiplygen is the general version of permutation multiplication.
	  * In the CDK Group package, only the same size permutations can be multiplied.
	  * That is why we neeed this function also. 
	  * 
	  * @param first a permutation
	  * @param second a permutation
	  * @return multiplication of these permutations. 
	  */
	 
	 public static Permutation multiplyGen(Permutation first, Permutation second) {
		 ArrayList<Permutation> perms= asc_order(first,second);
		 int[] firstValues= perms.get(0).getValues();
		 int[] secondValues= perms.get(1).getValues();
		 Permutation newPerm = new Permutation(perms.get(1));
	     for (int i = 0; i < firstValues.length; i++) {
	    	 newPerm.getValues()[i] = secondValues[firstValues[i]];
	     }
	     return newPerm;
	 }
	 
	 public static boolean ascCheck(ArrayList<Integer> list) {
		 boolean check= true;
		 for(int i=0;i<list.size()-1;i++) {
			 if(list.get(i)>list.get(i+1)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static int[] getTabloid(Permutation perm) {
		 int[] arr= new int[perm.size()/2];
	     for(int i=0;i<(perm.size()/2);i++){
		     arr[i]=perm.get(i);
		 }
		 return arr;
	 }
	 
	 public static boolean ascCheck(Permutation perm) {
		 boolean check=true;
		 for(int i=0;i<(perm.size()/2)-1;i++) {
			 if(perm.get(i)>perm.get(i+1)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static boolean inTheList(ArrayList<int[]> list, int[] arr) {
		 boolean check= false;
		 for(int i=0;i<list.size();i++) {
			 if(Arrays.equals(list.get(i), arr)) {
				 check=true;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static boolean inTheList(HashSet<ArrayList<Integer>> list, ArrayList<Integer> arr) {
		 boolean check= false;
		 for(ArrayList<Integer> l:list) {
			 if(l.equals(arr)) {
				 check=true;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static boolean inList(HashSet<HashSet<ArrayList<Integer>>> list, HashSet<ArrayList<Integer>> arr) {
		 boolean check= false;
		 for(HashSet<ArrayList<Integer>> l: list) {
			 if(l.equals(arr)) {
				 check=true;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Implementation of Orbit Fundamental Lemma. Molgen Book Example page 168.
	  * @param group a permutation group
	  * @param group2 a permutation group
	  * @return truncated tabloids for the orbits.
	  */
	 
	 public static  ArrayList<int[]> truncatedTabloids(PermutationGroup group, PermutationGroup group2) {
		 ArrayList<int[]> arrl= new ArrayList<int[]>();
	     //ArrayList<Permutation> lp= new ArrayList<Permutation>();
	     for(Permutation permutation: group.all()) {
	    	 for(Permutation permutation2:group2.all()) {
	    		 Permutation p=permutation.multiply(permutation2);
	    		 if(ascCheck(permutation)) {
	    			 int[] h= getTabloid(permutation);
	    			 if(!inTheList(arrl,h)) {
	    				 arrl.add(h);
	    				 //lp.add(permutation);
	    			 } 
	    		 }
	    	 } 
	     }
	     return arrl;
	 }
	 
	 
	 /**
	  * Gneerating list of unique orbits for the list of truncated tabloids
	  * and the acting group.
	  * 
	  * @param truncated truncated tabloids for the acting group
	  * @param group a permutation group 
	  * @return list of unique orbits
	  */
	 
	 public static HashSet<HashSet<ArrayList<Integer>>> fundamentalLemma(PermutationGroup group, PermutationGroup group2, ArrayList<Permutation> group3) {
		 HashSet<HashSet<ArrayList<Integer>>> orbits= new HashSet<HashSet<ArrayList<Integer>>>(); 
		 ArrayList<int[]> truncated= truncatedTabloids(group, group2);
		 for(int j=0;j<truncated.size();j++) {
		    HashSet<ArrayList<Integer>> orbit= new HashSet<ArrayList<Integer>>();
		    for(Permutation perm: group3) {
		    	ArrayList<Integer> l= new ArrayList<Integer>();
		    	for(int k=0;k<truncated.get(j).length;k++) {
		    		l.add(perm.get(truncated.get(j)[k]));
		    	}
		    	l.sort(ASC_ORDER);
		    	if(!inTheList(orbit,l)) {
		    		orbit.add(l);
		    	}
		    }
		    System.out.println(orbit);
		    if(!inList(orbits,orbit)) {
		    	orbits.add(orbit);
		    }
		  }
		 return orbits;
	 }
	 
	 /**
	  * Existence criteria for graphs from Grund's Thesis.
	  */
	 
	 /**
	  * In the existence criteria we need the degree list in descending order. 
	  * @param degrees list of degrees of a graphs
	  */
	 
	 public static void descendingOrder(ArrayList<Integer> degrees) {
		 degrees.sort(DES_ORDER);
	 }
	 
	 /**
	  * Checks whether the given degree list satisfies the existence criteria of graphs
	  * @param edges the number of edges.
	  * @param degrees the list of vertex degrees.
	  * @return true if satisfies the existence criteria.
	  */
	 public static boolean existenceCriteria(int edges, ArrayList<Integer> degrees) {
		 boolean check= false;
		 descendingOrder(degrees);
		 if(sum(degrees)==(2*edges) && degrees.get(0)<=edges){
			 check=true;
		 }
		 return check;
	 }
	 
	 
	 /**
	  * Graph connectivity
	  * @param vertices number of vertices
	  * @param edges number of edges
	  * @return true if connected; satisfying the criteria.
	  */
	 
	 public static boolean graphConnectivity(int vertices, int edges) {
		 boolean check= false;
		 if(edges>=vertices-1) {
			 check=true;
		 }
		 return check;
	 }
	 
	 /**
	  * Component Based Graph Existence Criteria from Grund's Thesis 2.1.6.
	  * @param components number of components
	  * @param vertices number of vertices
	  * @param edges number of edges
	  * @return true the graph exists with the degree distribution and (components-1) components
	  */
	 
	 public static boolean componentCriteria(int components, int vertices, int edges) {
		 boolean check=false;
		 if(components>1 && edges> vertices- components) {
			 check=true;
		 }
		 return check;
	 }
	 
	 /**
	  * Multigraph existence criteria
	  * @param degrees the list of degrees
	  * @return true if satisfies the multigraph existence criteria 
	  */
	 public static boolean multiGraphExistence(ArrayList<Integer> degrees) {
		 boolean check=false;
		 if(sum(degrees)>=2*(degrees.size()-1)) {
			 check=true;
		 }
		 return check;
	 }
	 
	 /**
	  * Simple graph existence criteria
	  * @param degrees the list of degrees
	  * @param edges number of edges
	  * @return true if satisfies the criteria.
	  */
	 
	 //TODO: Why if any they said in that theorem. Grund thesis 2.1.9. 
	 public static boolean simpleGraphExistence(ArrayList<Integer> degrees, int edges) {
		 boolean check= false;
		 if(sum(degrees)==(2*edges) && criteria4simple(degrees)) {
			 check=true;
		 }
		 return check;
	 }
	 
	 public static boolean criteria4simple(ArrayList<Integer> degrees) {
		 boolean check=true;
		 int sum=sum(degrees);
		 for(int i=0;i<degrees.size();i++) {
			 if(sum>i*(i-1)+minimalSum(degrees,i)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 public static int minimalSum(ArrayList<Integer> degrees, int j) {
		 int sum=0;
		 for(int i=j+1;i<degrees.size();i++) {
			 sum=sum+Math.min(j, degrees.get(i));
		 }
		 return sum;
	 }
	 
	 //Grund Thesis 2.1.10; existence of a simple graph. 
	 
	 public static boolean existenceCriteriaSimple(ArrayList<Integer> degrees) {
		 boolean check= false;
		 int sum=sum(degrees);
		 if(sum>=2*(degrees.size()-1)) {
			 check=true;
		 }
		 return check;
	 }
	 
	 //Grund thesis 2.2.6, the second criteria is coded. 
	 //TODO: We need also 2.2.8 satz from chapter 2. 
	 //TODO: Can be better, understand why we need V1 and other for this summation.
	 public static int multVectorSum(ArrayList<int[]> mults) {
		 int[] v1= buildMultArray(mults,1);
		 int[] v2= buildMultArray(mults,2);
		 int[] v3= buildMultArray(mults,3);
		 int sum=0;
		 for(int i=0;i<mults.size();i++) {
			 sum=sum+v1[i]+v2[i]+v3[i];
		 }
		 return sum;
	 }
	 
	 public static boolean multVectorCriteria1(ArrayList<int[]> mults) {
		 boolean check=false;
		 if(multVectorSum(mults)>=2*(mults.size()-1)) {
			 check=true;
		 }
		 return check;
	 }
	 
	 /**
	  * This array types are build for Lemma 2.2.6. An example is given
	  * in 2.2.4 
	  * @param arr list of multiplicities of vertices
	  * @param index the index requested for this array types. 
	  * @return multiplicity based arrays
	  */
	 public static int[] buildMultArray(ArrayList<int[]> arr, int index) {
		 int [] res= new int[arr.size()];
		 for(int i=0; i<arr.size();i++) {
			 res[i]=arr.get(i)[index];
		 }
		 return res;
	 }
	 
	 //TODO: How can we check the simple graphic criteria for these generated arrays ?
	 //TODO: 2.2.5 is not clear. Understand it.
	 public static void multVectorCriteria2(ArrayList<int[]> arr) {
		 int[] V1= buildMultArray(arr,1);
		 int[] V2= buildMultArray(arr,2);
		 int[] V3= buildMultArray(arr,3);
		 int[] V1V2 = sumArray(V1,V2);
		 int[] V1V3 = sumArray(V1,V3);
	 }
	 
	 public static boolean multVectorCriteria3(ArrayList<int[]> arr) {
		 boolean check=true;
		 int[] V1= buildMultArray(arr,1);
		 int[] V3= buildMultArray(arr,3);
		 for(int i=0;i<arr.size();i++) {
			 if(!(V1[i]==0 && V3[i]==0)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 
	 //Grund Thesis 1.1.16
	 /**
	  * Building the partitioning of the atom symbol occurences.
	  * @param list occurences of the atom symbols in the formula.
	  */
	 
	 public static ArrayList<Set<Integer>> blockSets(ArrayList<Integer> list) {
		 ArrayList<Set<Integer>> sets= new ArrayList<Set<Integer>>(); 
		 for(int i=0;i<list.size();i++) {
			 Set<Integer> set= new HashSet<Integer>();
			 set.add(blockSum1(list,i));
			 set.add(blockSum2(list,i));
			 sets.add(set);
		 }
		 return sets;
	 }
	 
	 public static int blockSum1(ArrayList<Integer> list, int index) {
		 int sum=1;
		 
		 if((index-1)==-1) {
			return sum;				 
		 }else {
			 for(int i=0;i<=index-1;i++) {
				 sum=sum+list.get(i);
			 }
		 }
		 return sum;
	 }
	 
	 public static int blockSum2(ArrayList<Integer> list, int index) {
		 int sum=0;
		 for(int i=0;i<=index;i++) {
			 sum=sum+list.get(i);
		 }
		 return sum;
	 }
	 
	 public static boolean stabParts(ArrayList<Set<Integer>> parts, Permutation perm) {
		 boolean check= true;
		 for(Set<Integer> part: parts) {
			 if(!act(part,perm).equals(part)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Grund 1.1.16
	  * @param list
	  * @return
	  */
	 public static ArrayList<Permutation> partitionYoungGroup(ArrayList<Integer> list) {
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 PermutationGroup grp= PermutationGroup.makeSymN(sum(list));
		 ArrayList<Set<Integer>> parts= blockSets(list);
		 for(Permutation perm: grp.all()) {
			 if(stabParts(parts,perm)) {
				 perms.add(perm);
			 }
		 }
		 return perms;
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
	 
	 
	 public static int Lsum(int[][] max, int i, int j, int size) {
		 int sum=0;
		 for(int k=j;k<size;k++) {
			 sum=sum+max[i][k];
		 }
		 return sum;
	 }
	 
	 
	 //Page 35, the bottom upper triangular matrices.
	 public static int Lsum2(int[][] A, int i, int j) {
		 int sum=0;
		 for(int k=0;k<j;k++) {
			 sum=sum+A[i][k];
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
	 
	 
	//Page 35, the bottom upper triangular matrices.
	 public static int Csum2(int[][] A, int i, int j) {
		 int sum=0;
		 for(int k=0;k<i;k++) {
			 sum=sum+A[k][j];
		 }
		 return sum;
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
	 
	 /**
	  * 3.2.7. Multiplicity matrix
	  * @param A  adjacency matrix
	  * @return
	  */
	 
	 public static int[][] multiplicityMat(int[][] A, int m){
		 int size = A.length; 
		 int[][] mult= new int[size][m];
		 for(int i=0;i<size;i++) {
			 for(int j=0;j<m;j++) {
				 mult[i][j]= countMult(A,i,j);
			 }
		 }
		 return mult;
	 }
	 
	 public static int countMult(int[][] A, int i, int j) {
		 int count=0;
		 int size= A.length;
		 for(int k=0;k<size;k++) {
			 if(A[i][k]==j) {
				 count++;
			 }
		 }
		 return count;
	 }
	 
	//Set diagonal to zeros.
	public static void zeroDiagonal(ArrayList[][] mat,int m) {
		for(int i=0; i<mat.length;i++) {
			for(int j=0;j<mat.length;j++) {
				if(i==j) {
					for(int k=0;k<m;k++) {
						mat[i][j].add(0);
					}
				}
			}
		}
	}
	 
	 /**
	  * 3.2.7 E matrix
	  * @param A adjacency matrix
	  * @param m maximal multiplicity
	  * @return e matrix as described in 3.2.7
	  */
	 public static ArrayList[][] eMatrix(int[][] mult, int m){
		 int size= mult.length;
		 ArrayList[][] eMatrix=TwoDimArrayList(size);
		 zeroDiagonal(eMatrix,m);
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) {
				 for(int k=0;k<m;k++) {
					 eMatrix[i][j].add(eMatrixEntry(mult,i,j,k));
					 eMatrix[j][i].add(eMatrixEntry(mult,i,j,k));
				 }
			 }
		 }
		 return eMatrix;
		 
	 }
	 
	 /**
	  * 3.2.7. E matrix entry
	  */
	 
	 public static int eMatrixEntry(int[][] mult,int i, int j, int k) {
		 if(mult[i][k]>0 && mult[j][k]>0) {
			 return 1;
		 }else {
			 return 0;
		 }
	 }
	 
	 
	 /** 
	  * Initialize 2D ArrayList.
	  */
	 
	 public static ArrayList[][] TwoDimArrayList(int size){
		 ArrayList[][] e= new ArrayList[size][size];
		 for(int i=0;i<size;i++) {
			 for(int j=0;j<size;j++) {
				 e[i][j]=new ArrayList<Integer>();
			 }
		 }
		 return e;
	 }
	 
	 /**
	  * Print 2D ArrayList
	  */
	 
	 public static void print2DArrayList(ArrayList[][] list) {
		 int size= list[0].length;
		 for(int i=0;i<size;i++) {
			 for(int j=0;j<size;j++) {
				 System.out.println(list[i][j]);
			 }
		 }
	 }
	 
	 /**
	  * Check a multiplicity arraylist has any non zero entry
	  */
	 
	 public static boolean zeroCheck(ArrayList<Integer> list) {
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
	  * Maximal Multiplicity
	  */
	 public static int maxMult(ArrayList<Integer> list) {
		 int index=0;
		if(zeroCheck(list)) {
			return index;
		}else {
			for(int i=(list.size())-1;i>=0;i--) {
				 if(list.get(i)!=0) {
					 index=index+i;
					 break;
				 }
			 }
			 return index+1;
		}
	 }
	 
	 /**
	  * Maximal matrix based on E matrix entries for multiplicities. 3.2.7
	  */
	 
	 public static int[][] maximalMultiplicityMatrix(int[][] mult, int m){
		 int size= mult.length;
		 ArrayList[][] eMat= eMatrix(mult,m);
		 int[][] max= new int[size][size];
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) {
				 max[i][j]=max[j][i]=maxMult(eMat[i][j]);
			 }
		 }
		 return max;
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
	 
	 /**
	  * 3.2.7 L matrix entry.
	  */
	 
	 public static int multLEntry(int i, int j, int k, ArrayList[][] eMatrix, int[][] mult) {
		 int count=0;
		 int size= eMatrix.length;
		 for(int s=j+1;s<size;s++) {
			 count=count+(int)eMatrix[i][s].get(k);
		 }
		 return Math.min(mult[i][k], count);
	 }
	 

	 /**
	  * 3.2.7 L matrix entry.
	  */
	 
	 public static int multCEntry(int i, int j, int k, ArrayList[][] eMatrix, int[][] mult) {
		 int count=0;
		 int size= eMatrix.length;
		 for(int s=i+1;s<size;s++) {
			 count=count+(int)eMatrix[s][j].get(k);
		 }
		 return Math.min(mult[j][k], count);
	 }
	 
	 /**
	  * 3.2.7 L matrix
	  */
	 public static ArrayList[][] multLMatrix(int[][]mult, int m){
		 int size= mult.length;
		 ArrayList[][] eMatrix=eMatrix(mult,m);
		 ArrayList[][] LMatrix=TwoDimArrayList(size);
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) {
				 for(int k=0;k<m;k++) {
					 LMatrix[i][j].add(multLEntry(i, j, k, eMatrix, mult));
					 LMatrix[j][i].add(multLEntry(i, j, k, eMatrix, mult));
				 }
			 }
		 }
		 return LMatrix;
	 }
	 
	 /**
	  * 3.2.7 Linverse entry
	  */
	 
	 public static int multL2Entry(int i, int j, int k , int[][]A, int[][]mult) {
		 int count=0;
		 for(int s=0; s<j;s++) {
			 count=count+indicatorFunction(A[i][s],k);
		 }
		 return mult[i][k]-count;
	 }
	 
	 /**
	  * 3.2.7 L Inverse Matrix 
	  */
	 
	 public static ArrayList[][] multL2Matrix(int[][]A, int m){
		 int size= A.length;
		 int[][] mult=multiplicityMat(A,m);
		 ArrayList[][] L2Matrix=TwoDimArrayList(size);
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) { //TODO: i+1 ?
				 for(int k=0;k<m;k++) {
					 L2Matrix[i][j].add(multL2Entry(i, j, k, A, mult));
				 }
			 }
		 }
		 return L2Matrix;
	 }
	 
	 /**
	  * 3.2.7 C matrix
	  */
	 public static ArrayList[][] multCMatrix(int[][] mult, int m){
		 int size= mult.length;
		 ArrayList[][] eMatrix= eMatrix(mult,m);
		 ArrayList[][] CMatrix= TwoDimArrayList(size);
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) { //TODO: i+1 ?
				 for(int k=0;k<m;k++) {
					 CMatrix[i][j].add(multCEntry(i, j, k, eMatrix, mult));
				 }
			 }
		 }
		 return CMatrix;
	 }
	 
	 /**
	  * 3.2.7 Cinverse entry
	  */
	 
	 public static int multC2Entry(int i, int j, int k , int[][]A, int[][]mult) {
		 int count=0;
		 for(int s=0; s<i;s++) {
			 count=count+indicatorFunction(A[s][j],k);
		 }
		 return mult[j][k]-count;
	 }
	 
	 /**
	  * 3.2.7 C2 matrix
	  */
	 
	 public static ArrayList[][] multC2Matrix(int[][]A, int m){
		 int size= A.length;
		 int[][] mult=multiplicityMat(A,m);
		 ArrayList[][] C2Matrix=TwoDimArrayList(size);
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) { //TODO: i+1 ?
				 for(int k=0;k<m;k++) {
					 C2Matrix[i][j].add(multC2Entry(i, j, k, A, mult));
				 }
			 }
		 }
		 return C2Matrix;
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
	 
	 //Filling the adjacency matrices. The criterias are given in 3.2.1.
	 public static int entry(int[][] A,int[][] max, int[][] L, int[][] C, int[][]L2, int[][]C2, int i, int j) {
		 int number=0;
		 int minimal= Math.min(max[i][j],Math.min(L2[i][j],C2[i][j]));
		 for(int k=minimal;k>=minimal;k--) {
			 if((L2[i][j]-k<=L[i][j]) && (C2[i][j]-k<=C[i][j])) {
				 number+=k;
				 A[i][j]=number;
				 break;
			 }
		 }
		 return number;
	 }
	 
	 
	 
	 
	// Algorithm 3.2.9
    public static void canonicalMatrix2(int[][] mult,int m) throws IOException{
    	int size    = mult.length;
    	int[][] A   = new int[size][size];
    	ArrayList[][] eMat= eMatrix(mult,m);
    	int[][] max = maximalMultiplicityMatrix(mult,m);
		ArrayList[][] L   = multLMatrix(mult,m);
		ArrayList[][] L2   = multLMatrix(A,m);
		ArrayList[][] C   = multCMatrix(mult,m);
		ArrayList[][] C2   = multCMatrix(A,m);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forward2(mult,A,max,L,C,indices,m);
	}
		
    /**
     * 3.2.9 Forward criteria 
     */
    
    public static boolean forBackCriteria(int h,int l, int c, int l2, int c2, int m) {
    	boolean check=true;
    	if(!(l2>0 && c2>0)) {
    		check=false;
    	}else {
    		for(int k=1;k<m;k++) {
        		if(!(l2-indicatorFunction(h,k)<=l && c2-indicatorFunction(h,k)<=c)) {
        			check=false;
        			break;
        		}
        	}
    	}
    	return check;
    }
    
	//3.2.9. Step forward
	public static void forward2(int[][] mult,int[][] A, int[][]max, ArrayList[][]L, ArrayList[][]C, ArrayList<Integer> indices, int m) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int maximum= max[i][j];
	 	for(int h=maximum;h>0;h--) {
	 		int l2= multL2Entry(i,j,h,A,mult);
			int c2= multC2Entry(i,j,h,A,mult);
			int l = (int)L[i][j].get(h);
			int c = (int)C[i][j].get(h);
	 		if(forBackCriteria(h,l,c,l2,c2,m)) {
	 			A[i][j]=A[j][i]=h;
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backward2(mult,A, max,L, C, indices,m);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				forward2(mult,A, max, L, C, modified,m);
	 			}
	 		}else {
	 			backward2(mult,A, max, L, C, indices,m);
	 		}
	 	}
	}
		 
	//3.2.9 Step backward
	public static void backward2(int[][] mult,int[][] A, int[][]max, ArrayList[][]L, ArrayList[][]C, ArrayList<Integer> indices, int m) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		if(i==max.length-2 && j==max.length-1) {
			count++;
			writeMatrix(A);
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int maximal= A[i][j];
				if(maximal>0) {
					for(int h=maximal;h>0;h--) {
						int l2= multL2Entry(i,j,h,A,mult);
						int c2= multC2Entry(i,j,h,A,mult);
						int l = (int)L[i][j].get(h);
						int c = (int)C[i][j].get(h);
						if(forBackCriteria(h,l,c,l2,c2,m)) {
							ArrayList<Integer> modified2=successor(modified,max.length);
							forward2(mult,A, max, L, C, modified2,m);
						}else {
							backward2(mult, A, max, L, C, modified,m);
						}
					}
				}
			}else {
				backward2(mult,A,max,L,C,indices,m);
			}
		}
	}
	  
	 public static void printMatrix(int[][] mat) {
		for (int[] row : mat) {	 
			System.out.println(Arrays.toString(row)); 
	    } 
	 }
	 
	 
	 /**
	  * Algorithm 3.1.6: Canonical Generation Algorithm
	  */
	 
	 
	 
	 public static void writeMatrix(int[][] mat) throws IOException {
		 fileWriter.write(String.format("Adjacency matrix - %d", output.size()));
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
	 
	 public static void writeLine(int index, String s, int[][] mat) throws IOException {
		 fileWriter.write(String.format("Matrix row - %d %s", index, s));
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
	 
	 public static void writePartition(int index, ArrayList<Integer> old, ArrayList<Integer> newPart) throws IOException {
		 fileWriter.write(String.format("New and old partition - %d", index));
		 fileWriter.newLine();
		 fileWriter.write(old +" "+ newPart);
		 fileWriter.newLine();
		 fileWriter.write("----------");
		 fileWriter.newLine();
		 fileWriter.flush();
	 }
	 public static ArrayList<String> atomList(String molecularFormula) {
		 ArrayList<String> symbols= new ArrayList<String>();
		 String[] atoms = molecularFormula.split("(?=[A-Z])");
	     for (int i=0;i<atoms.length;i++) {
	    	 String[] info = atoms[i].split("(?=[0-9])", 2);
	         for(int k=0;k<(info.length > 1 ? Integer.parseInt(info[1]):1);k++) {
	        	 symbols.add(info[0]);
	         }
	     }
	     symbols.sort(DES_STRING_ORDER);
	     return symbols;
	 }
	 
	 
	 public static ArrayList<Integer> degrees(ArrayList<String> atomList) {
		 ArrayList<Integer> degrees= new ArrayList<Integer>();
		 for(String atom: atomList) {
			 degrees.add(valences.get(atom));
		 }
	     return degrees;
	 }
	 
	 public static void generate (String molecularFormula, String filedir) throws IOException {
		 ArrayList<Integer> degrees= degrees(atomList(PermutationGroupFunctions.molecularFormula));
		 if(verbose) System.out.println("Generating adjacency matrices for"+" "+PermutationGroupFunctions.molecularFormula+" "+"...");
		 canonicalMatrix(degrees);
		 if(verbose) System.out.println(count+" "+"adjacency matrices were generated and written to the given file.");
	 }
	 
	 
	 //Definition 3.2.7
	 public static int indicatorFunction(int x, int k) {
		 if(x==k) {
			return 1; 
		 }else {
			return 0;
		 }
	 }
	 
	 
	 
	 /**
	  * Canonization functions Chapter 3.3
	  */
	 
	 
	 /**
	  * Example 3.3.9. Permutation action on a vector
	  */
	 public static int[] permutationAction(int[] vector, Permutation perm) {
		 int[] v2= new int[vector.length];
		 for(int i=0;i<vector.length;i++) {
			 v2[i]=vector[perm.get(i)];
		 }
		 return v2;
	 }
	 
	 /**
	  * Definition 3.3.2
	  * Partition generation
	  */
	 public static ArrayList<Integer> partitionRule(ArrayList<Integer> partEx, int degree){
		 ArrayList<Integer> partNew = new ArrayList<Integer>();
		 if(partEx.get(degree-1)>1) { //Indices starts with zero that is why -1
			 addOnes(partNew,degree);
			 partNew.add(partEx.get(degree-1)-1);
			 for(int k=degree;k<partEx.size();k++) {
				 partNew.add(partEx.get(k));
			 }
		 }else if(partEx.get(degree-1)==1) {
			 addOnes(partNew,degree);
			 for(int k=degree+1;k<partEx.size();k++) {
				 partNew.add(partEx.get(k));
			 }
		 }
		 return partNew;
	 }
	 
	 
	 /**
	  * Clean and tested functions for Grund Thesis 
	 * @throws CloneNotSupportedException 
	  */
	
	 public static IAtomContainer mat2AC(int[][] matrix) throws CloneNotSupportedException {
		 IAtomContainer ac = new AtomContainer();
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));

		 int size=matrix.length;
		 for(int i=0;i<size;i++) {
			 for(int j=i+1;j<size;j++) {
				 if(matrix[i][j]!=0) {
					 ac=addBond(ac,i,j,matrix[i][j]);
				 }
			 }
		 }
		 return ac;
	 }
	 
	 public static IAtomContainer addBond(IAtomContainer ac, int index, int index2, int order) throws CloneNotSupportedException {
		 IAtomContainer ac2= ac.clone();
		 if(order==1) {
			 ac2.addBond(index,index2, Order.SINGLE);
		 }else if(order==2) {
			 ac2.addBond(index,index2, Order.DOUBLE);
		 }else if(order==3) {
			 ac2.addBond(index,index2, Order.TRIPLE);
		 }
		 return ac2;
	 }
	 
	/** 
	 * Grund Thesis - Algorithm 3.2.3 (DONE)
	 * Simple multigraph generation
	 * @param degrees atom valences  
	 * @throws IOException
	 */
	 
	public static void matrixGenerator(ArrayList<Integer> degrees) throws IOException{
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
		forward(degrees,A,max,L,C,indices); //It was originally without partition entry.
	}
	
	public static void canonicalMatrixGenerator(ArrayList<Integer> degrees, ArrayList<Integer> partition) throws IOException{
		int size    = degrees.size();
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] C   = upperTriangularC(degrees);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forwardCanonical(degrees,partition,A,max,L,C,indices); //It was originally without partition entry.
	}
	/**
	 * This algorithm performs canonical test with permutation tree of Youngsubgroup chain.
	 * @param degrees degree array for atoms
	 * @param partition partition of atom types
	 * @throws IOException
	 */
	
	public static void canMatGen(ArrayList<Integer> degrees, ArrayList<Integer> partition) throws IOException{
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
		Permutation id= new Permutation(size);
		permTree.add(id);
		setFile();
		forwardCanRep(degrees,partition,A,max,L,C,indices); //It was originally without partition entry.
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
	 			System.out.println(Arrays.deepToString(A)+" "+i+" "+j);
	 			System.out.println("forward"+" "+i+" "+j);
	 			//TODO: I need to check two things:
	 			/**
	 			 * First, I need to check whether I am in the same line or move to next one, sometimes it goes back and no need to check for that case.
	 			 * Second, this partition check is when we reach the last column and need to move to next row. 
	 			 */
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backward(degrees,A, max,L, C, indices);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				System.out.println(modified.get(0)+" "+modified.get(1)+" "+i+" "+j);
	 				//TODO: 
	 				/**
	 				 * if i new > i old and j=the last column; so we finish one row then moved to next one. Then we need to check canonicity.
	 				 */
	 				//System.out.println(Arrays.deepToString(A));
	 				//System.out.println(i+" "+j+" "+modified+" "+"indices forward successor");
	 				forward(degrees,A, max, L, C, modified);
	 			}
	 		}else {
	 			backward(degrees,A, max, L, C, indices);
	 		}
	 	}
	}
	
	public static void forward(ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
	 	int minimal= Math.min(max[i][j],Math.min(l2,c2));
	 	for(int h=minimal;h>=0;h--) {
	 		if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
	 			System.out.println("general"+" "+i+" "+j+" "+h+" "+Arrays.deepToString(A));
	 			A[i][j]=A[j][i]=h;
	 			System.out.println("generalAfter"+" "+i+" "+j+" "+h+" "+Arrays.deepToString(A));
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backwardDemo(degrees,partition,A, max, L, C, indices);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				System.out.println("modified"+" "+modified.get(0)+" "+modified.get(1)+" "+i+" "+j+" "+A.length);
	 				if(modified.get(0)>i && j==A.length-1) { 
	 					System.out.println("forward"+" "+modified.get(0)+" "+modified.get(1)+" "+i+" "+j);
	 					System.out.println("partition"+" "+partition);
	 					partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					System.out.println("canonical partition"+" "+partition);
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					System.out.println("row clean"+" "+Arrays.toString(A[i]));
	 					A[i]=canonicalRow(group,partition,A[i]);
	 					System.out.println("row canonical"+" "+Arrays.toString(A[i]));
	 					System.out.println("canonicalCheck"+" "+canonicalRowCheck(A[i],group,partition,A.length));
	 					System.out.println("A[i] can"+" "+Arrays.toString(A[i]));
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						//A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						System.out.println("refined partition"+" "+partition);
	 						forward(degrees, partition, A, max, L, C, modified);
	 					}
	 				}else {
	 					forward(degrees, partition, A, max, L, C, modified);
	 				}
	 					
		 			/**if(j==A.length-1) {
	 					/**partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					//System.out.println("before canonical"+" "+partition+" "+Arrays.toString(A[i]));
	 					//System.out.println(Arrays.toString(A[i])+" "+"canBefore"+" "+partition);
	 					System.out.println("partition"+" "+partition);
	 					System.out.println(i+" "+Arrays.deepToString(A));
	 					A[i]=canonicalRow(group,partition,A[i]);
	 					System.out.println(i+" "+Arrays.deepToString(A)+" "+"after");
	 					//System.out.println(Arrays.toString(A[i])+" "+"canAfter");
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						//System.out.println("canonical");
	 						//A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						//System.out.println("refined"+" "+partition);
	 						ArrayList<Integer> modified=successor(indices,max.length);
	 		 				forward(degrees,partition,A, max, L, C, modified);
	 					}
	 				}else {
	 					ArrayList<Integer> modified=successor(indices,max.length);
		 				forward(degrees,partition,A, max, L, C, modified);
	 				}**/
	 			}
	 		}else {
	 			backwardDemo(degrees, partition,A, max, L, C, indices);
	 		}
	 	}
	}
	
	public static void forwardCanonical(ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
	 	int minimal= Math.min(max[i][j],Math.min(l2,c2));
	 	for(int h=minimal;h>=0;h--) {
	 		if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
	 			System.out.println("general"+" "+i+" "+j+" "+h+" "+Arrays.deepToString(A));
	 			A[i][j]=A[j][i]=h;
	 			System.out.println("generalAfter"+" "+i+" "+j+" "+h+" "+Arrays.deepToString(A));
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backwardCanonical(degrees,partition,A, max, L, C, indices);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				System.out.println("modified"+" "+modified.get(0)+" "+modified.get(1)+" "+i+" "+j+" "+A.length);
	 				if(modified.get(0)>i && j==A.length-1) { 
	 					System.out.println("forward"+" "+modified.get(0)+" "+modified.get(1)+" "+i+" "+j);
	 					System.out.println("row clean"+" "+Arrays.toString(A[i]));
	 					System.out.println("normal"+" "+Arrays.toString(A[i]));
	 					//A[i]=canonicalRow(group,partition,A[i]);
	 					System.out.println("row canonical"+" "+Arrays.toString(A[i]));
	 					//System.out.println("canonicalCheck"+" "+canonicalRowCheck(A[i],group,partition,A.length));
	 					System.out.println("A[i] can"+" "+Arrays.toString(A[i]));
	 					System.out.println("partition"+" "+partition);
	 					partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					System.out.println("canonical partition"+" "+partition);
	 					if(canonicalRowCheck(A[i],i,partition,A.length)) {
	 						//A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						System.out.println("refined once part"+" "+partition);
	 						partition=refinedPartitioning(partition,A[i]);
	 						System.out.println("refined Array"+" "+Arrays.toString(A[i])+" "+i);
	 						System.out.println("refined partition"+" "+partition);
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
	public static BufferedWriter file;
	public static void setFile() throws IOException {
		 BufferedWriter f = new BufferedWriter(new FileWriter("C:\\Users\\mehme\\Desktop\\output.txt"));
		 file=f;
	}
	public static ArrayList<Permutation> permTree= new ArrayList<Permutation>(); 
	public static void forwardCanRep(ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		ArrayList<Integer> cyclePart= new ArrayList<Integer>();
		for(Integer i: partition) {
			cyclePart.add(i);
		}
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
	 	int minimal= Math.min(max[i][j],Math.min(l2,c2));
	 	for(int h=minimal;h>=0;h--) {
	 		if((l2-h<=L[i][j]) && (c2-h<=C[i][j])) {
	 			A[i][j]=A[j][i]=h;
	 			if(i==(max.length-2) && j==(max.length-1)) {
	 				backwardCanRep(degrees,partition,A, max, L, C, indices);
	 			}else {
	 				ArrayList<Integer> modified=successor(indices,max.length);
	 				if(modified.get(0)>i && j==A.length-1) {
	 					file.write("part before"+" "+partition+"\n");
	 					//System.out.println("part before"+" "+partition);
	 					partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					file.write(i+" "+"part canonical"+" "+partition+" "+Arrays.toString(A[i])+"\n");
	 					//System.out.println(i+" "+"part canonical"+" "+partition+" "+Arrays.toString(A[i]));
	 					file.write("permTree"+" "+permTree+"\n");
	 					if(canonicalRepRowCheck(A[i],i,partition,A.length)) {
	 						partition=refinedPartitioning(partition,A[i]);
	 						file.write("refined"+" "+partition+"\n");
	 						//System.out.println("refined"+" "+partition);
	 						ArrayList<Permutation> cycles= permTreeCycles(cyclePart,partition);
	 						file.write("cycles"+" "+cycles+"\n");
	 						//System.out.println("cycles"+" "+cycles);
	 						permTree=cycles;
	 						//updatePermTree(cycles);
	 						//System.out.println("permTree"+" "+permTree);
	 						forwardCanRep(degrees, partition, A, max, L, C, modified);
	 					}
	 				}else {
	 					forwardCanRep(degrees, partition, A, max, L, C, modified);
	 				}
	 			}
	 		}else {
	 			backwardCanRep(degrees, partition,A, max, L, C, indices);
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
	
	public static void backward(ArrayList<Integer> degrees,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
		int i=indices.get(0);
		int j=indices.get(1);
		int l2= LInverse(degrees,i,j,A);
		int c2= CInverse(degrees,i,j,A);
		if(i==max.length-2 && j==max.length-1) {
			
			/**
			 * When I just add the matrix to the list it did not do anything.
			 * This is not a proper solution but that is ok for now to keep it like this.
			 */
			
			int[][] mat2= new int[A.length][A.length]; 
			for(int k=0;k<A.length;k++) {
				for(int l=0;l<A.length;l++) {
					mat2[k][l]=A[k][l];
				}
			}
			//writeLine(111,"b",mat2);
			System.out.println(Arrays.deepToString(mat2));
			//output.add(mat2);
			//System.out.println("bu"+" "+Arrays.deepToString(mat2));
			//System.out.println("bu"+" "+Arrays.deepToString(A));
			//writeMatrix(A);
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
					//writeLine(i,"b",A);
					backward(degrees, A, max, L, C, modified);
				}
			}
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
			//writeLine(111,"b",mat2);
			System.out.println("done");
			System.out.println(Arrays.deepToString(A));
			//writeMatrix(A);
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					System.out.println(Arrays.deepToString(A)+" "+"back"+" "+i);
					ArrayList<Integer> modified2=successor(modified,max.length);
	 				forwardCanonical(degrees, partition, A, max, L, C, modified2);
	 				
					//forward(degrees,partition,A, max, L, C, modified2);
					/**if(j==A.length-1) {
	 					partition=canonicalPartition(i,partition);
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						ArrayList<Integer> modified2=successor(indices,max.length);
	 		 				forward(degrees,partition,A, max, L, C, modified2);
	 					}
	 				}else {
	 					ArrayList<Integer> modified2=successor(modified,max.length);
						forward(degrees,partition,A, max, L, C, modified2);
	 				}**/
				}else {
					backwardCanonical(degrees,partition, A, max, L, C, modified);
				}
			}
		}
	}
	
	public static void backwardCanRep(ArrayList<Integer> degrees,ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
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
			//writeLine(111,"b",mat2);
			System.out.println("done");
			System.out.println(Arrays.deepToString(A));
			//writeMatrix(A);
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					ArrayList<Integer> modified2=successor(modified,max.length);
	 				forwardCanRep(degrees, partition, A, max, L, C, modified2);
	 				
					//forward(degrees,partition,A, max, L, C, modified2);
					/**if(j==A.length-1) {
	 					partition=canonicalPartition(i,partition);
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						ArrayList<Integer> modified2=successor(indices,max.length);
	 		 				forward(degrees,partition,A, max, L, C, modified2);
	 					}
	 				}else {
	 					ArrayList<Integer> modified2=successor(modified,max.length);
						forward(degrees,partition,A, max, L, C, modified2);
	 				}**/
				}else {
					backwardCanRep(degrees,partition, A, max, L, C, modified);
				}
			}
		}
	}
	
	public static ArrayList<int[][]> output= new ArrayList<int[][]>();
	public static void backward(ArrayList<Integer> degrees,ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices) throws IOException {
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
			//writeLine(111,"b",mat2);
			System.out.println("done");
			System.out.println(Arrays.deepToString(A));
			//writeMatrix(A);
		}else{
			ArrayList<Integer> modified=predecessor(indices, max.length);
			i= modified.get(0);
			j= modified.get(1);
			if(i>0 && j-i==1) {
				int x= A[i][j];
				if(x>0 && (l2-(x-1)<=L[i][j]) && (c2-(x-1)<=C[i][j])) {
					A[i][j]=A[j][i]=x-1;
					System.out.println(Arrays.deepToString(A)+" "+"back"+" "+i);
					ArrayList<Integer> modified2=successor(modified,max.length);
					if(modified2.get(1)>j) {
						System.out.println("backward"+" "+modified2.get(1)+" "+i+" "+j);
	 					partition=canonicalPartition(i,partition); //TODO: Might need to test again
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					A[i]=canonicalRow(group,partition,A[i]);
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						//A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						forward(degrees, partition, A, max, L, C, modified);
	 					}
	 				}else {
	 					forward(degrees, partition, A, max, L, C, modified2);
	 				}
					//forward(degrees,partition,A, max, L, C, modified2);
					/**if(j==A.length-1) {
	 					partition=canonicalPartition(i,partition);
	 					PermutationGroup group=getYoungGroup(partition,A.length);
	 					if(canonicalRowCheck(A[i],group,partition,A.length)) {
	 						A[i]=canonicalRow(A[i],group,partition); //TODO: Need to check whether it modifies the A matrix really.
	 						partition=refinedPartitioning(partition,A[i]);
	 						ArrayList<Integer> modified2=successor(indices,max.length);
	 		 				forward(degrees,partition,A, max, L, C, modified2);
	 					}
	 				}else {
	 					ArrayList<Integer> modified2=successor(modified,max.length);
						forward(degrees,partition,A, max, L, C, modified2);
	 				}**/
				}else {
					backward(degrees,partition, A, max, L, C, modified);
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
	  * Inverse partition functions.
	  */
	 /**
	  * Line criteria partitioning to implement on refined partition (DONE)
	  * @param refinedPartition refined partition based on matrix row.
	  * @return
	  */
	 
	 public static ArrayList<Integer> partitionLineCriteria(ArrayList<Integer> refinedPartition){
		 ArrayList<Integer> partNew = new ArrayList<Integer>();
		 addOnes(partNew,2);
		 if(refinedPartition.get(1)>1) {
			 partNew.add(refinedPartition.get(1)-1);
			 for(int k=2;k<refinedPartition.size();k++) {
				 partNew.add(refinedPartition.get(k));
			 }
		 }else if(refinedPartition.get(1)==1){
			 for(int k=2;k<refinedPartition.size();k++) {
				 partNew.add(refinedPartition.get(k));
			 }
		 }
		 return partNew;
	 }
	 
	 /**
	  * This desCheck is used to check whether the blocks of a row is in
	  * descending order or not. Rather than creating another function
	  * this function also returns the descending order check and terminates
	  * if the block is not in descending order. If the desCheck is not true
	  * then the matrix is not canonical.
	  */
	 
	 public static boolean desCheck=true; 
	 
	 /**
	  * Refined Partition - Grund 3.3.11
	  * Create the refined partition from the new row. The occurences of
	  * the numbers in the blocks of the row. 
	  * @param partition Partition
	  * @param row int[] row of a matrix to check canonicality. 
	  * @return
	  */
	 
	 /**public static ArrayList<Integer> refinedPartitioning(ArrayList<Integer> partition, int[] row){
		 System.out.println("refined ici"+" "+partition+" "+Arrays.toString(row));
		 ArrayList<Integer> refined= new ArrayList<Integer>();
		 int index=0;
		 int count=1;
		 for(Integer p:partition) {
			 System.out.println("subpart"+" "+p);
			 System.out.println("sub indices"+" "+index+" "+(p+index));
			 for(int i=index;i<p+index;i++) {
				 System.out.println("row[i] and row[i+1] "+" "+row[i]+" "+row[i+1]);
				 if(i+1<p+index) {
					 if(row[i]==row[i+1]) {
						 count++;
					 }else if(row[i]>row[i+1]){
						 refined.add(count);
						 count=1;
					 }else{
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
	 }**/
	 
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
	  * (DONE) 3.3.8. Canonical Permutation 
	  * The action should act on a matrix line based on partitioning
	  * Should return also a canonical permutation; the one puting the
	  * row entries in descending order with respoect to the blocks. 
	  * @param group PermutationGroup created by the partition info
	  * @param partition Partition
	  * @param array row of a matrix
	  * @return
	  */
	 
	 public static ArrayList<Permutation> canonicalPermutation(PermutationGroup group, ArrayList<Integer> partition, int[] array){
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 int[] canCheck = new int[array.length];
		 for(Permutation perm: group.all()) {
			 canCheck = actArray(array,perm);
			 //System.out.println("canCheck"+" "+Arrays.toString(canCheck));
			 if(blockwiseDesOrderCheck(canCheck,partition)) {
				 perms.add(perm);
				 System.out.println(Arrays.toString(canCheck)+" "+perm);
			 }
			 
		 }
		 return perms;
	 }
	 
	 
	 public static int[] canonicalRow(PermutationGroup group, ArrayList<Integer> partition, int[] array){
		 int[] canCheck = new int[array.length];
		 for(Permutation perm: group.all()) {
			 canCheck = actArray(array,perm);
			 if(blockwiseDesOrderCheck(canCheck,partition)) {
				 break;
			 }
			 
		 }
		 return canCheck;
	 }
	 
	 /**public static Permutation canonicalPermutation(PermutationGroup group, ArrayList<Integer> partition, int[] array){
		 Permutation canonical= null;
		 System.out.println("first"+" "+Arrays.toString(array));
		 int[] canCheck = new int[array.length];
		 for(Permutation perm: group.all()) {
			 canCheck = act(array,perm);
			 System.out.println("canonical Act"+" "+Arrays.toString(canCheck));
			 if(blockBasedEqual(array,canCheck,partition)) {
				 System.out.println("original"+" "+Arrays.toString(array));
				 System.out.println("canonical"+" "+Arrays.toString(canCheck));
				 System.out.println("perm"+" "+perm.toCycleString());
				 canonical=perm;
			 }
			 break;
		 }
		 return canonical;
	 }**/
	 
	 public static int[] canonicalRow(int[] array,PermutationGroup group, ArrayList<Integer> partition){
		 int[] canRow = new int[array.length];
		 for(Permutation perm: group.all()) {
			 canRow = actArray(array,perm);
			 if(desBlockCheck(partition,canRow) && blockBasedEqual(array,canRow,partition)) {
				 break;
			 }
		 }
		 return canRow;
	 }
	 
	 /**
	  * To check whether the blocks of a row array is in descending order.
	  * @param partition refined partition 
	  * @param row row of a matrix
	  * @return
	  */
	 
	 public static boolean desBlockCheck(ArrayList<Integer> partition, int[] row){
		 boolean check=true;
		 int index=0;
		 for(Integer p:partition) {
			 for(int i=index;i<p+index;i++) {
				 if(i+1<p+index) {
					 if(row[i]<row[i+1]) {
						 check=false;
						 break;
					 }
				 }
			 }
			 index=index+p;
		 }
		 return check;
	 }
	 
	 //TODO: DesBlockCheck is wrong. Not one by one all the members of the block but blockwise check should be.
	 public static boolean desBlockCheck(ArrayList<Integer> partition, int[] row, int[] modified){
		 boolean check=true;
		 int index=0;
		 System.out.println("desBlockCheck");
		 for(Integer p:partition) {
			 System.out.println("range"+" "+index+" "+(p+index));
			 for(int i=index;i<p+index;i++) {
				 System.out.println("org and modified"+" "+row[i]+" "+modified[i]);
				 if(row[i]<modified[i]) {
					 check=false;
					 break;
				 }
			 }
			 index=index+p;
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
			 System.out.println("indices"+" "+index+" "+(p+index));
			 System.out.println("canonical");
			 System.out.println("getBlock"+" "+Arrays.toString(getBlocks(canonical,index,p+index)));
			 System.out.println("Int"+" "+toInt(getBlocks(canonical,index,p+index)));
			 System.out.println("original");
			 System.out.println("getBlock"+" "+Arrays.toString(getBlocks(original,index,p+index)));
			 System.out.println("IntCan"+" "+toInt(getBlocks(canonical,index,p+index)));
			 System.out.println("IntOrg"+" "+toInt(getBlocks(original,index,p+index)));
			 if(toInt(getBlocks(canonical,index,p+index))<toInt(getBlocks(original,index,p+index))) {
				 check=false;
				 break;
			 }
			 index=index+p;
		 }
		 return check;
	 }
	 
	 public static boolean descBlockCheck(ArrayList<Integer> partition, int[] canonical, int[] original){
		 boolean check=true;
		 int index=0;
		 for(Integer p:partition) {
			 Integer[] can= getBlocks(canonical,index,p+index);
			 Integer[] org= getBlocks(original,index,p+index);
			 if(!Arrays.equals(can,org)) {
				 if(toInt(can)>toInt(org)) {
					 check=true;
					 break;
				 }else if(toInt(can)<toInt(org)) {
					 check=false;
					 break;
				 }
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
	 
	 /**
	  * Grund 1.1.16 (Done)
	  * Generating the permutation group for a partition.
	  * Here, the generators of all the symmetry groups are taken. 
	  * This function returns us the Young subgroup of a partition. 
	  * @param partition int partition.
	  * @return PermutationGroup Young subgroup
	  */
		
	public static PermutationGroup getPermutationGroup(ArrayList<Integer> partition, int total) {
		int[] part=listToArray(partition);
		List<int[]> parts=subPartitions(part);
		for(int[] i:parts) {
			//System.out.println(Arrays.toString(i));
		}
		List<Permutation> generators= new ArrayList<Permutation>();
		for(int[] p: parts) {
			//System.out.println(getGenerators(p,total));
			generators.addAll(getGenerators(p,total));
		}
		return generateGroup(generators,total);
	}
	
	public static int[] listToArray(ArrayList<Integer> list) {
		int[] array = new int[list.size()];
		int count = 0;
		for(Integer l:list) {
			array[count]=l.intValue();
			count++;
		}
		return array;
	}
	
	/**
	* 3.3.11 Canonical Line Test (Done)
	* @param A int matrix
	* @param partition given degree partition
	* @param total sum of degree distribution
	* @return
	*/
	 
	/**public static boolean canonicalLineTest(int[][] A, ArrayList<Integer> partition, int total) {
		boolean check=true;
		for(int i=0; i<A[0].length;i++) {
			PermutationGroup group=getPermutationGroup(partition,total);
			ArrayList<Permutation> canonicalPerm=canonicalPermutation(group,A[i]);
			if(canonicalPerm.size()!=0) {
				A[i]=act(A[i],canonicalPerm.get(0));
				partition=partition(refinedPartitioning(partition,A[i]),2);
				A[i]=act(A[i],canonicalPerm.get(0));
			}else {
				check=false;
				break;
			}
		}
		return check;
	}**/
	
	
	//TODO: There should be the canonical permutation not descending but equal
	// As given in 3.3.9; the vectors are equal with respect to block elements.
	//In 3.3.11 the canonization permutation is defined with respect to the partition,
	// I did not understand where I can use that canonical perm equal to the original row.
	
	public static boolean canonicalLineTest(int[][] A, ArrayList<Integer> partition, int total) {
		boolean check=true;
		partition=partitionWDegree(partition,2);
		for(int i=0; i<A[0].length;i++) {
			PermutationGroup group=getPermutationGroup(partition,total);
			/**
			 * We prefered to have ArrayList return since there is no canonical
			 * permutation for some cases. We check with the list size.
			 */
			if(canonicalRowCheck(A[0],group,partition)) {
				partition=refinedPartitioning(partition,A[i]);
			}else {
				check=false;
				break;
			}
		}
		return check;
	}
	
	/**
	 * Grund 3.3.9 
	 * To check two arrays are equal blockwise.
	 * We need that check after a canonical perm act on the row.
	 * @param row row of a matrix
	 * @param row2 row of a matrix
	 * @param partition Degree partition
	 * @return
	 */
	
	public static boolean blockBasedEqual(int[] row, int[] row2, ArrayList<Integer> partition) {
		boolean check=true;
		int index=0;
		for(Integer p:partition) {
			if(!Arrays.equals(desOrderBlock(row,index,index+p), desOrderBlock(row2,index,index+p))) {
				check=false;
				break;
			}
			index=index+p;
		}
		return check;
	}
	
	public static boolean blockwiseDesOrderCheck(int[] row, ArrayList<Integer> partition) {
		boolean check=true;
		int index=0;
		for(Integer p: partition) {
			//System.out.println("Partition" +" "+p);
			//System.out.println("Indices" +" "+(index)+" "+(index+p));
			if(!desOrderCheck(row,index,index+p-1)) {
				check=false;
				break;
			}
			index=(index+p);
			//System.out.println("index changed"+" "+index);
		}
		return check;
	}
	
	public static boolean desOrderCheck(int[] row, int begin, int end) {
		boolean check=true;
		for(int i=begin;i<end;i++) {
			if(row[i]<row[i+1]) {
				check=false;
				break;
			}
		}
		return check;
	}
	
	public static Integer [] desOrderBlock(int[] row, int begin, int end) {
		Integer[] sub=subArray(row,begin,end);
		Arrays.sort(sub, DES_ORDER);
		return sub;
	}
	
	public static Set<Integer> subSet(int[] array, int begin, int end) {
		Set<Integer> sub= new HashSet<Integer>();
		for(int i=begin;i<end;i++) {
			sub.add(array[i]);
		}
		return sub;
	}
	
	public static Integer[] subArray(int[] array, int begin,int end) {
		return IntStream.range(begin, end).mapToObj(i->array[i]).toArray(Integer[]::new);
	}
	/**
	 * Grund 3.1.6. Canonical Test Step 5
	 * 
	 * As explained in 3.1.6, I should check with all the permutations 
	 * of the automorphism group. Descending order check is based on partition
	 * as given in 3.3.8.
	 * @param row row of a matrix
	 * @param group Permutation group built with respect to the partition
	 * @return
	 */
	
	/**
	 * Int i added just because of 
	 * @throws IOException 
	 **/
	
	public static boolean canonicalRowCheck(int[] row,int i, ArrayList<Integer> partition, int total) throws IOException {
		boolean check=true;
		PermutationGroup group=getYoungGroup(partition,total);
		for(Permutation perm: group.all()) {
			//file.write("canRowCheck"+" "+perm+"\n");
			//System.out.println("canRowCheck"+" "+perm);
			//if(!desBlockCheck(partition,row,actArray(row,perm))) {
			//if(!desBlockwiseCheck(partition,row,actArray(row,perm))) {
			if(!descBlockCheck(partition,row,actArray(row,perm))) {	
				check=false;
				break;
			}
		}
		return check;
	}
	
	/**
	 * Canonical Row Check is performed by checking the descending order with every permutation
	 * of a Young subgroup built with respect to a partition. Rather than dealing with that many
	 * permutations, we get the Young subgroup representatives. These representatives are used
	 * in the canonical row check. 
	 * @param row row of connectivity matrix
	 * @param i row index
	 * @param reps Young group representatives from the previous partition.
	 * @param cycles cycles are generated with respect to the previous and current partition at the ith row.
	 * @param partition current partition at the ith row.
	 * @param total number of atoms. 
	 * @return
	 * @throws IOException 
	 */
	
	public static boolean canonicalRepRowCheck(int[] row, int i, ArrayList<Integer> partition, int total) throws IOException {
		boolean check=true;
		if(i==0) { //In the first row, descending order check is performed with its Young subgroup
			canonicalRowCheck(row,i,partition,total);
		}else {  
			for(Permutation perm: permTree) {
				file.write("permTr"+" "+perm.toCycleString()+"\n");
			}
			for(Permutation perm: permTree) {
				if(!descBlockCheck(partition,row,actArray(row,perm))) {
					check=false;
					break;
				}
			}
		}
		return check;
	}
	
	public static ArrayList<Integer> canonicalPartition(int i, ArrayList<Integer> partition){
		if(i==0) {
			return partitionWDegree(partition,2);
		}else {
			return partitionWDegree(partition,1);
		}
	}
	
	/**
	* 3.3.11 Canonical Line Test (Done)
	* @param A int array
	* @param partition given degree partition
	* @param total
	* @return
	*/
	 
	public static boolean canonicalLineTest(int[] A, ArrayList<Integer> partition, int total) {	
		boolean check=true;
		PermutationGroup group=getPermutationGroup(partition,total);
		if(!canonicalRowCheck(A,group,partition)) {
			check=false;
		}
		return check;
	}
	
	
	/**
	 * ****************************************************
	 */
	
	
	 public static ArrayList<Integer> partition(ArrayList<Integer> part, int degree){
		 for(int i=1;i<=degree;i++) {
			 ArrayList<Integer> part2 = partitionRule(part,i);
			 part=part2;
		 }
		 return part;
	 }
	 /**
	  * add number of 1s into a ArrayList
	  */
	 
	 public static void addOnes(ArrayList<Integer> list, int number) {
		 for(int i=0;i<number;i++) {
			 list.add(1);
		 }
	 }
	 
	 public static ArrayList<Integer> decreaseOne(ArrayList<Integer> list) {
		 ArrayList<Integer> list2= new ArrayList<Integer>();
		 for(int i=0;i<list.size();i++) {
			 list2.add(list.get(i)-1);
		 }
		 return list2;
	 }
	 	 
	 /**
	  * Definition 1.1.16 (DONE)
	  * 
	  */
	 
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
	  * Definition 1.1.16 (DONE)
	  */
	 
	 public static int[] subPartition(int[] partition, int index){
		 int[] sub= new int[0];
		 int former = sum(partition,index-1)+1;
		 System.out.println(former+" "+"former");
		 int latter = sum(partition,index);
		 System.out.println(latter+" "+"later");
		 if(former==latter) {
			 sub=addElement(sub,former);
		 }else {
			 for(int i=former;i<=latter;i++) {
				 sub=addElement(sub,i);
			 }
		 }
		 return sub;
	 }
	 
	 public static List<int[]> subPartitions(int[] partition){
		 List<int[]> parts = new ArrayList<int[]>();
		 for(int i=0;i<partition.length;i++) {
			parts.add(subPartition(partition,i)); 
		 }
		 return parts;
	 }
	 
	 public static ArrayList<ArrayList<Integer>> subPartitionsList(ArrayList<Integer> partition){
		 ArrayList<ArrayList<Integer>> parts = new ArrayList<ArrayList<Integer>>();
		 for(int i=0;i<partition.size();i++) {
			parts.add(subPartition(partition,i)); 
		 }
		 return parts;
	 }
	 
	 public static int[] toIntArray(ArrayList<Integer> set) {
		 int[] array= new int[set.size()];
		 for(int i=0;i<set.size();i++) {
			 array[i]=set.get(i);
		 }
		 return array;
	 }
	 
	 
	 public static boolean automorphismPermutationCheck(Permutation perm, ArrayList<Integer> part) {
		 boolean check=true;
		 for(int[] list: subPartitions(part)) {
			 list=decreaseOne(list);
			 if(!list.containsAll(act(list,perm))) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	  	 
	 /**
	  * Descending order check for two arrays
	  */
	 
	 public static boolean desOrderCheck(int[] arr1, int[] arr2) {
		 boolean check=true;
		 int size= arr1.length;
		 for(int i=0;i<size;i++) {
			 if(arr1[i]<arr2[i]) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 
	 /**
	  * 3.3.11. Getting second degree partitioning of the second degree partition (DONE)
	  */
	 
	 public static ArrayList<Integer> lineCriteria(ArrayList<Integer> part, int degree){
		 ArrayList<Integer> newPart= subPartition(part,degree-1);
		 return partition(partition(newPart,1), 2);
	 }
	 
	 /**
	  * Canonical Test 3.3.11 (This is when we have group elements are arrayList.
	  * But if we keep it as PermutationGroup, it will be easier for us to access
	  * to the leftTraversals.
	  */
	 
	 public static boolean canonicalTest(int[][] A, ArrayList<Integer> part) {
		 boolean check=true;
		 for(int i=1;i<A.length;i++) {
			 ArrayList<Permutation> perms=group4Partition(lineCriteria(part, i-1));
			 if(!groupActionVectors(A[i],perms)) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * 3.3.8 Canonizor Permutations
	  */
	 
	 public static ArrayList<Permutation> group4Partition(ArrayList<Integer> part){
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 PermutationGroup group= PermutationGroup.makeSymN(sum(part));
		 for(Permutation perm: group.all()) {
			 if(automorphismPermutationCheck(perm,part)) {
				 perms.add(perm);
			 }
		 }
		 return perms;
	 }
	 
	 /**
	  * 3.3.11 Permutation Group Action on a Matrix Row
	  */
	 public static boolean groupActionVectors(int[] a, ArrayList<Permutation> perms) {
		 boolean check=true;
		 for(Permutation perm:perms) {
			 if(!desOrderCheck(a,permutationAction(a,perm))) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * 3.1.1. A matrix strip 
	  */
	 
	 public static int[][] matrixStrip(int[][] A , int index, ArrayList<Integer> parts){
		 int size= A.length;
		 int[][] mat= new int[size][size];
		 ArrayList<Integer> subPart= subPartition(parts,index);
		 for(Integer i: subPart) {
			 mat[i]= A[i];
		 }
		 return A;
	 }
	 
	 /**
	  * 3.1.5 Automorphism group of a matrix strip
	  */
	 	 
	 public static ArrayList<Permutation> stripAuto(int[][] A, int index, ArrayList<Integer> part){
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 for(int deg=1; deg<=index;deg++) {
			 int[][] strip= matrixStrip(A,deg,part);
			 for(Permutation perm:group4Partition(part)) {
				 if(autoCheckMatStrip(strip,perm)) {
					 perms.add(perm);
				 }
			 } 
		 }
		 return perms;
	 }
	 
	 /**
	  * Actin g a permutation on matrix strips.
	  * @param A
	  * @param index
	  * @param perm
	  * @return
	  */
	 
	 public static boolean autoCheckMatStrip(int[][] A,Permutation perm) {
		 boolean check=true;
		 int size=A.length;
		 for(int i=0;i<size;i++) {
			 if(!A[i].equals(permutationAction(A[i],perm))) {
				 check=false;
				 break;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Centralizer group chain and other functions between 3.3.1 and 3.3.11
	  */
	 
	 /**
	  * Definition 3.3.3
	  * In the definition, l(i) is not clearly explained. With my assumption;
	  * we need to take the former partition and index as inputs.
	  * Sum the entries of the partition until the index (including the index value) 
	  */
	 
	 public static int LValue(ArrayList<Integer> partEx, int degree) {
		 return (sum(partEx,degree-1)-degree+1); 
	 }
	 
	 /**
	  * Definition 3.3.3.
	  * The number of left cosets of a subgroup in a group.
	  * This number is called the index of the subgroup in the bigger group.
	  * This index relies on Lagrange's Index Theorem.
	  */
	 
	 public static int LIndex(ArrayList<Integer> subPart, ArrayList<Integer> part) {
		 return youngGroupOrder(part)/youngGroupOrder(subPart);
	 }
	 
	 
	 /**
	  * Grund 3.3.3 Group tree cycles for cosets
	  * @param index i power of the partition
	  * @param partition partition
	  * @return List<int[]> list of components as described in Grund Thesis 3.3.3
	  */
	 
	 public static List<int[]> cycleComponents(int index, ArrayList<Integer> partition){
		 List<int[]> components = new ArrayList<int[]>();
		 int lValue= LValue(partition,index);
		 for(int j=1;j<=lValue;j++) {
			 int[] component=new int[2];
			 component[0]=index;
			 component[1]=j+index-1;
			 components.add(component);
		 }
		 return components;
	 }
	 
	 /**
	  * In the case of 3.3.19 example, they take the cycle representatives
	  * one by one for each row, that is why we have this table representation.
	  * ArrayList<Permutation> are stored as ArrayList
	  */
	 
	 public static ArrayList<ArrayList<Permutation>> representatives = new ArrayList<ArrayList<Permutation>>();
	 
	 public static ArrayList<ArrayList<Integer>> partition2 = new ArrayList<ArrayList<Integer>>();
	 /**
	  * (After Grund Answers) First it builds cycle representatives as given 
	  * in 3.3.3. Then the Aut group (the perm tree) is updated. For the table
	  * representation like in 3.3.19, representatives of all the rows are
	  * stored in representatives ArrayList
	  * @param index row index
	  * @param formerPartition former partition
	  * @param indexChanges from which index the partition is divided to new partition
	  * @param size number of atoms
	  * @return
	  */
	 
	 public static ArrayList<Permutation> cycleRepresentatives(int index, ArrayList<Integer> indexChanges, int size) {
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 perms.add(idPermutation(size));
		 int lValue= LValue(refinedPartitions.get(index),index);
	     for(int i=1;i<=lValue;i++) {
	    	 int[] values= idValues(size);
			 for(int changes:indexChanges) {
				 values[changes]=changes+i;
				 values[changes+i]=changes;
			 }
			 Permutation p = new Permutation(values);
			 perms.add(p);
		 }
	     representatives.add(index, perms);
	     updatePermTree(perms);
		 return perms;
	 }
	 
	 public static ArrayList<Integer> findChanges(ArrayList<Integer> extPart, ArrayList<Integer> part) {
		 ArrayList<Integer> changes= new ArrayList<Integer>();
		 int extS= extPart.size();
		 int newS= part.size();
		 if(extS!=newS) {
			 outer:
			 for(int i=0; i<extS;i++) {
				 for(int j=0;j<newS;j++) {
					 if(extPart.get(i)!=part.get(j)) {
						 if(i==(extS-1)) {
							 changes.add(sum(extPart));
							 break outer;
						 }else {
							 changes.add(sum(extPart,i));
							 j++;
						 }
					 }
				 }
			 }
		 }
		 return changes;
	 }
	 
	 
	 /**
	  * (After Grund Answers) First it builds cycle representatives as given 
	  * in 3.3.3. Then the Aut group (the perm tree) is updated.For the table
	  * representation like in 3.3.19, representatives of all the rows are
	  * stored in representatives ArrayList.
	  * 
	  * Different from the former version, here, we generate all the cycles 
	  * for the index as in 3.3.19. The last partition is the one does not
	  * have the indexth entry so we need to build all the cycles like in
	  * 3.3.3 and 3.3.19. 
	  * @param index
	  * @param lastartition
	  * @param size
	  * @return
	  */
	 
	 public static ArrayList<Permutation> cycleRepresentatives(int index, int size) {
		 ArrayList<Permutation> perms= new ArrayList<Permutation>();
		 perms.add(idPermutation(size));
		 int lValue= LValue(refinedPartitions.get(index),index);
	     for(int i=1;i<=lValue;i++) {
	    	 int[] values= idValues(size);
	    	 values[index]=index+i;
	    	 values[index+i]=index;
			 Permutation p = new Permutation(values);
			 perms.add(p);
		 }
	     representatives.add(index, perms);
	     updatePermTree(perms);
		 return perms;
	 }
	
	 /**
	  * Build id permutation
	  */
	 
	 public static Permutation idPermutation(int size) {
		 Permutation perm= new Permutation(size);
		 return perm;
	 }
	 
	 /**
	  * (After Grund's Answer)
	  * I need the strip ( block) based canonical test.
	  * So taking a strip to check canonicity not line by line.
	  */
	 
	 //TODO: The xth partition or all the partitions can be stored in an ArrayList
	 public static boolean canonicalBlockTest(int[] row,int r, int rowIndex) {
		 boolean check=true;
		 if(rowIndex==findY(r)) {
			 check=yCanonicalRowCheck(row,rowIndex);
		 }else {
			 ArrayList<Permutation> perms= blockPermutations((r+1),rowIndex); //TODO: Without using perms, just arrays ?
			 ArrayList<Integer> formerPartition= partition2.get(rowIndex-1);
			 for(Permutation perm : perms) {
					if(!desBlockwiseCheck(formerPartition,row,actArray(row,perm))) {
						check=false;
						break;
					}
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * Refined Partitions a input partition eklenmeli.
	  * @param A
	  * @throws IOException
	  */
	 public static boolean firstBlock(int[][]A) throws IOException {
		 boolean check=true;
		 ArrayList<Integer> partition= inputPartition;
		 int len= inputPartition.get(0);
		 int matSize= A.length;
		 for(int i=0;i<len;i++) {
			 partition=canonicalPartition(i,refinedPartitions.get(i));
			 if(canonicalRowCheck(A[i],i,partition,matSize)) {
				 refinedPartitions.add(refinedPartitioning(partition,A[i]));
				 representatives.add(cycleRepresentatives(i,findChanges(refinedPartitions.get(i),refinedPartitions.get(i+1)),matSize));
			 }else {
				 check=false;
				 break;
			 }
		 } 
		 if(check) {
			 fillRepresentatives(matSize);
		 }
		 return check;
	 }
	 
	 //TODO: For the other strips, I need to update these representatives and refinedPartitions.
	 public static void fillRepresentatives(int total) {
		 int z=findZ(0);
		 for(int i=z;i<total;i++) {
			 refinedPartitions.add(partitionWDegree(refinedPartitions.get(i),1));
			 representatives.add(cycleRepresentatives(i,total));
		 }
	 }
	
	 //TODO: fOR THE Y ROW, I will check with the former ones but then How can I decide what the canonical perm ?
	 public static boolean yCanonicalRowCheck(int[] rowY, int yValue) {
		 boolean check=true;
		 ArrayList<Permutation> reps= representatives.get(yValue); //TODO: Without using perms, just arrays ?
		 ArrayList<Integer> formerPartition= partition2.get(yValue-1);
		 for(Permutation perm: reps) {
				if(!desBlockwiseCheck(formerPartition,rowY,actArray(rowY,perm))) {
					check=false;
					break;
				}
			}
		 return check;
	 }
	
	 public static int partitionSize;
	 public static ArrayList<Integer> inputPartition= new ArrayList<Integer>();

	 /**
	  * According to 3.3.1 in Grund Thesis, we need to fill a matrix 
	  * strip by strip. For a strip r, first fill the matrix upto z value. 
	  * Filling the rth strip first. Then test the strip whether canonical
	  * or not based on the representatives. 
	  */
	 
	public static void canonicalBlockGenerator(ArrayList<Integer> degrees, ArrayList<Integer> partition) throws IOException{
		int size    = degrees.size();
		partitionSize=partition.size();
		inputPartition=partition;
		int r=1;
		int[][] A   = new int[size][size];
		int[][] max = maximalMatrix(degrees);
		int[][] L   = upperTriangularL(degrees);
		int[][] C   = upperTriangularC(degrees);
		ArrayList<Integer> indices= new ArrayList<Integer>();
		indices.add(0);
		indices.add(1);
		zeroDiagonal(A);
		forwardCanonicalBlock(degrees,partition,A,max,L,C,indices,r); //It was originally without partition entry.
	}
	
	public static ArrayList<ArrayList<Integer>> refinedPartitions= new ArrayList<ArrayList<Integer>>();
	public static void forwardCanonicalBlock(ArrayList<Integer> degrees, ArrayList<Integer> partition,int[][] A, int[][]max, int[][]L, int[][]C, ArrayList<Integer> indices, int r) throws IOException {
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
	 					if(canonicalBlockTest(A[i],r,i)) { //Based on former perms, check canonical or not then add new perms if it is canonical.
	 						partition=refinedPartitioning(partition,A[i]);
	 						refinedPartitions.add(partition);
	 						ArrayList<Integer> changes=findChanges(refinedPartitions.get(i-1),partition); //TODO: In case i=1 or 0, the first partition is used.
	 						cycleRepresentatives(i,changes,max.length);
	 						if(i==findZ(r)) {
	 							r++;
	 							forwardCanonicalBlock(degrees, partition, A, max, L, C, modified,r);
	 						}else {
	 							forwardCanonicalBlock(degrees, partition, A, max, L, C, modified,r);
	 						}
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
	  * 3.3.18 From the yth line, to ith representatives,
	  * multiply them all to find the permutation set.
	  */
	 
	 public static ArrayList<Permutation> blockPermutations(int yValue, int rowIndex){
		 ArrayList<Permutation> yReps = representatives.get(yValue);	
		 for(int i=(yValue+1);i<=rowIndex;i++) {
			 ArrayList<Permutation> newReps= new ArrayList<Permutation>();
			 for(Permutation perm : representatives.get(i)) {
				 for(Permutation yPerm:yReps) {
					 Permutation n=yPerm.multiply(perm);
					 if(!newReps.contains(n)) {
						 newReps.add(n);
					 } 
				 }
			 }
			 yReps=newReps;
		 }
		 return yReps;
	 }
	 
	 /**
	  * 3.3.3. definition is not clear. Based on Example 3.3.19, I defined this function.
	  * The cycles are build based on the changes between the ext and new refined partitions. 
	  */
	 
	 //TODO: Clean the code
	 public static ArrayList<Permutation> permTreeCycles(int level, ArrayList<Integer> extPart, ArrayList<Integer> part){
		 int total=sum(part);
		 int LValue=LValue(extPart,level);
		 ArrayList<Permutation> cycles= new ArrayList<Permutation>();
		 cycles.add(new Permutation(idValues(total)));
		 if(lastPartition(extPart) && lastPartition(part)) {
			 if(LValue!=1) {
				 if(extPart.size()!=part.size()) {
					 int index=sum1s(part,(part.size()-1));
					 for(int i=index+1;i<LValue+index;i++) {
						 cycles.add(buildCycles(index,i,total));
					 } 
				 }
			 }
		 }else {
			 int index=0;
			 List<Integer> indices= new ArrayList<Integer>();
			 if(LValue!=1) {
				 for(Integer i: extPart) {
					 if(i!=part.get(index)) {
						 indices.add(sumPrev(part,index));
						 index++;
						 if(indices.size()==LValue) {
							 break;
						 }
					 }
					 index++;
				 } 
				 cycles.add(buildPermutation(indices,total));
			 }
			 
		 }
		 return cycles;
	 }
	 
	 public static Permutation buildCycles(int former, int latter, int total) {
		 int[] entries=idValues(total);
		 entries[former]=latter;
		 entries[latter]=former;
		 Permutation perm=new Permutation(entries); 
		 return perm;
	 }
	 
	 public static Permutation buildPermutation(List<Integer> indices, int total) {
		 Permutation cycle= new Permutation(buildEntriesArray(indices,total));
		 return cycle;
	 }
	 
	 public static int[] buildEntriesArray(List<Integer> indices,int total) {
		 int[] entries=idValues(total);
		 for(Integer i:indices) {
			 entries[i]=i+1;
			 entries[i+1]=i;
		 }
		 return entries;
	 }
	 public static boolean lastPartition(ArrayList<Integer> partition) {
		 boolean check=true;
		 int count=0;
		 for(Integer i:partition) {
			 if(i!=1) {
				 if(count==1) {
					 check=false;
					 break;
				 }else {
					 count++;
				 }
			 }
		 }
		 return check;
	 }
	 
	 public static int sum1s(ArrayList<Integer> part, int index) {
		 int result=0;
		 for(int i=0;i<index;i++) {
			 result=result+part.get(i);
		 }
		 return result-1;
	 }
	 
	 public static int sumPrev(ArrayList<Integer> part, int index) {
		 int result=0;
		 for(int i=0;i<=index;i++) {
			 result=result+part.get(i);
		 }
		 return result-1;
	 }
	 
	 /**
	  * Grund 3.3.5. Coset components as described for Ui-1/Ui
	  * Here, the Ui is ith centralizer of the automorphism.
	  * It returns the list of permutation components of centralizers.
	  * @param Ui The ith centralizer group
	  * @param partition degree partition
	  * @param index centralizer index 
	  * @return List<Permutation> list of permutation components.
	  */
	 
	 public static List<Permutation> cosetComponents(List<Permutation> Ui,ArrayList<Integer> partition, int index){
		 List<Permutation> perms= new ArrayList<Permutation>();
		 int lValue= LValue(partition,index);
		 for(int j=2;j<lValue;j++) {
			 perms.add(Ui.get(j));
		 }
		 return perms;
	 }
	 
	 /**
	  * 3.3.1. Permutation tree
	  */
	 
	 public static ArrayList<Permutation> permTreeRep(ArrayList<Integer> part, int degree){
		 ArrayList<Permutation> out= new ArrayList<Permutation>();
		 int l=LValue(partition(part,degree-1),degree);
		 for(int j=1;j<=l;j++) {
			 int[] perm= new int[2];
			 perm[0]=degree;
			 perm[1]=j+degree-1;
			 out.add(new Permutation(perm));
		 }
		 return out;
	 }
	 
	 /**
	  * Centralizer group of a given index as described in 3.3.1. Grund Thesis.
	  * @param group PermutationGroup U
	  * @param index index of the subgroup
	  */
	 
	 //TODO: The action of U groups on rows and matrix stips must be analized. 
	 public static ArrayList<Permutation> centralizer(PermutationGroup group, int index) {
		ArrayList<Permutation> centralizer= new ArrayList<Permutation>();
		//TODO: This index might be changed later. 0 to index or index-1. 
		for(Permutation perm: group.all()) {
			if(stabilizerCheck(index,perm)) {
				centralizer.add(perm);
			}
		}
		return centralizer;
	 }
	 
	 public static boolean stabilizerCheck(int index, Permutation p) {
		 boolean check=true;
		 for(int i=0; i<index;i++) {
			 if(i!=p.get(i)) {
				 check=false;
			 }
		 }
		 return check;
	 }
	 
	 /**
	  * 3.3.7 Left Coset indices J(i) for vij values. 
	  */
	 
	 /**
	  * Understand and test this one to build V groups.
	  */
	 
	 public static ArrayList<Integer> leftClassRepIndices(PermutationGroup group, int i, int[][]A, int index, ArrayList<Integer> part) {
		 ArrayList<Integer> vReps= new ArrayList<Integer>();
		 List<Permutation> perms=group.getLeftTransversal(i);
		 int[][] Ar= matrixStrip(A,index,part);
		 ArrayList<Permutation> reps=permTreeRep(part,index);
		 for(int j=0;j<perms.size();j++) {
			 for(int k=0;k<reps.size();k++) {
				 if(autoCheckMatStrip(Ar,perms.get(j).multiply(reps.get(k)))) {
					 vReps.add(j);
				 }
			 }
			 
		 }
		 return vReps;
	 }
	 
	 /**
	  * 3.3.7 Must take a matrix strip as the input. The matrix strip function is not needed.
	  */
	 
	 public static ArrayList<Permutation> leftClassRepIndices(PermutationGroup group, int  index, int[] strip) {
		 ArrayList<Permutation> vReps= new ArrayList<Permutation>();
		 List<Permutation> perms=group.getLeftTransversal(index);
		 for(Permutation perm:perms) {
			for(Permutation perm2:perms) {
				if(Arrays.equals(strip, act(strip,perm2.multiply(perm)))){
					vReps.add(perm2);
				}
			}
		 }
		 return vReps;
	 }
	 
	 /**
	  * 3.3.7 Left Coset indices vij permutations. 
	  */
	 
	 public static ArrayList<Permutation> leftClassRepPerms(PermutationGroup group, int i, int[][]A, int index, ArrayList<Integer> part) {
		 ArrayList<Permutation> vReps= new ArrayList<Permutation>();
		 List<Permutation> perms=group.getLeftTransversal(i);
		 int[][] Ar= matrixStrip(A,index,part);
		 ArrayList<Permutation> reps=permTreeRep(part,index);
		 for(int j=0;j<perms.size();j++) {
			 for(int k=0;k<reps.size();k++) {
				 Permutation v= perms.get(j).multiply(reps.get(k));
				 if(autoCheckMatStrip(Ar,v)) {
					 vReps.add(v);
				 }
			 }
		 }
		 return vReps;
	 }
	 
	 /**
	  * Canonical strips - 3.3.14
	  */
	 
	 /**
	  * I made a mistake in the former version. r is the block index,
	  * so the input for the function.
	  */
	 
	 public static int xValue(ArrayList<Integer> part, int rValue) {
		 int sum=0;
		 for(int i=0;i<=(rValue-1);i++) { //TODO = or not; check again
			 sum=sum+part.get(i);
		 }
		 return sum;
	 }
	 
	 public static int yValue(ArrayList<Integer> part, int rValue) {
		 return xValue(part,rValue)+1; // the xValue can be fixed, no need to calculate again.
	 }
	 
	 public static int zValue(ArrayList<Integer> part, int rValue) {
		 return xValue(part,rValue)+part.get(rValue);
	 }
	 
	 public static void setFileWriter() throws IOException {
		 PermutationGroupFunctions.fileWriter = new BufferedWriter(new FileWriter(filedir));
	 }
	 
	 //TODO: Adjacency matrix to atomcontainer.
	 private void parseArgs(String[] args) throws ParseException, IOException{
		 Options options = setupOptions(args);	
		 CommandLineParser parser = new DefaultParser();
		 try {
			 CommandLine cmd = parser.parse(options, args);
			 PermutationGroupFunctions.molecularFormula = cmd.getOptionValue("molecularFormula");
			 PermutationGroupFunctions.filedir = cmd.getOptionValue("filedir");
			 setFileWriter();
			 if (cmd.hasOption("verbose")) PermutationGroupFunctions.verbose = true;
		  } catch (ParseException e) {
			  // TODO Auto-generated catch block
			  HelpFormatter formatter = new HelpFormatter();
			  formatter.setOptionComparator(null);
			  String header = "\nGenerates adjacency matrices for a given molecular formula."
						+ " The input is a molecular formula string."
						+ "For example 'C2H2O'."
						+ "Besides this formula, the directory is needed to be specified for the output"
						+ "file. \n\n";
			  String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory";
			  formatter.printHelp( "java -jar AlgorithmicGroupTheory.jar", header, options, footer, true );
			  throw new ParseException("Problem parsing command line");
		   }
		}
		
		private Options setupOptions(String[] args){
			Options options = new Options();
			Option molecularFormula = Option.builder("f")
				     .required(true)
				     .hasArg()
				     .longOpt("molecularFormula")
				     .desc("Molecular formula as a string (required)")
				     .build();
			options.addOption(molecularFormula);
			Option verbose = Option.builder("v")
				     .required(false)
				     .longOpt("verbose")
				     .desc("Print messages about matrix generation")
				     .build();
			options.addOption(verbose);	
			Option filedir = Option.builder("d")
				     .required(false)
				     .hasArg()
				     .longOpt("filedir")
				     .desc("Creates and store the output txt file in the directory (required)")
				     .build();
			options.addOption(filedir);
			return options;
	 }
		
	public static String inchiGeneration(IAtomContainer molecule) throws CDKException {
		String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(molecule).getInchi();	
		return inchi;
	}
	
	 public static void main(String[] args) throws CloneNotSupportedException, CDKException, IOException {   
		 // TODO Auto-generated method stub	
		 
		 ArrayList<Integer> degrees= new ArrayList<Integer>();
		 
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(4);
		 degrees.add(5);
		 degrees.add(1);
		 
		 
		 ArrayList<Integer> partition= new ArrayList<Integer>();
		 
		 partition.add(10);
		 partition.add(1);
		 partition.add(1);
		 
		 ArrayList<Integer> part= new ArrayList<Integer>();
		 part.add(1);
		 part.add(1);
		 part.add(2);
		 part.add(1);
		 part.add(1);
		 part.add(4);
		 
		 ArrayList<Integer> part2= new ArrayList<Integer>();
		 part2.add(1);
		 part2.add(1);
		 part2.add(1);
		 part2.add(1);
		 part2.add(1);
		 part2.add(1);
		 part2.add(4);
		 
		 
		 int[] row= new int[5];
		 row[0]=0;
		 row[1]=3;
		 row[2]=0;
		 row[3]=0;
		 row[4]=1;
		 
		 //System.out.println(part);
		 PermutationGroup group1=getYoungGroup(part,10);
		 //System.out.println(part2);
		 PermutationGroup group2=getYoungGroup(part2,10);
		 /**for(Permutation perm:group1.transversal(group2)) {
			 System.out.println(perm.toCycleString());
		 }**/
		 
		 System.out.println(permTreeCycles(3,part,part2));
		
		 /**
		 for(Permutation perm1: group1.all()) {
			 System.out.println("bu"+" "+perm1.toCycleString());
			 for(Permutation perm2: group2.all()) {
				 System.out.println(perm2.multiply(perm1));
			 }
		 }**/
		 
		 //canMatGen(degrees,partition);
		 //file.close();
		 //PermutationGroup hep= PermutationGroup.makeSymN(3);
		 //System.out.println(hep);
		 /**for(Permutation perm:hep.all()) {
			 //System.out.println(perm.toCycleString());
		 }**/
		 
		 /**for(int[] ar: cycleComponents(1,simp2)) {
			 System.out.println(Arrays.toString(ar));
		 }**/
		 
		 //System.out.println(finalInt);
		 //desBlockwiseCheck(simp2,arr);
		 //PermutationGroup grp= PermutationGroup.makeSymN(3);
		 //System.out.println(grp);
		 //canonicalMatrixGenerator(simple,simp2);
		 //PermutationGroupFunctions.fileWriter = new BufferedWriter(new FileWriter("C:\\Users\\mehme\\Desktop\\test.txt"));
		 //canonicalMatrixGenerator(simple,simp2);
		 //System.out.println(can.toCycleString());
		 //matrixGenerator(simple);
		 //canonicalMatrixGenerator(simple,simp2);
		 
		 /**IAtomContainer ac = new AtomContainer();
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("O"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));

	     ArrayList<IAtomContainer> acs= new ArrayList<IAtomContainer>();
	     ArrayList<String> inc= new ArrayList<String>();
	     for(int i=0;i<output.size();i++) {
	    	 //System.out.println(i+" "+Arrays.deepToString(output.get(i)));
	    	 IAtomContainer acm=mat2AC(output.get(i));
	    	 depict(acm,"C:\\Users\\mehme\\Desktop\\output\\"+i+".png");
	    	 String in= inchiGeneration(acm);
	    	 if(!inc.contains((String)in)){
	    		 inc.add(in);
	    		 acs.add(acm);
	    	 }
	    	 //depict(mat2AC(output.get(i)),"C:\\Users\\mehme\\Desktop\\"+i+".png");
	     }
	     int c=0;
	     for(IAtomContainer a: acs) {
	    	 depict(a,"C:\\Users\\mehme\\Desktop\\"+c+".png");
	    	 c++;
	     }
	     System.out.println(acs.size());**/
		 //list.add(4);
		 //list.add(5);
		 
		 /**for(Permutation perm:getYoungsubgroupGenerators(list,6)) {
			 System.out.println("young"+" "+perm.toCycleString());
		 }**/
		 
		 
		 /**IAtomContainer ac = new AtomContainer();
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("N"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
		 //PermutationGroup group=getPermutationGroup(list2,5);
		 //System.out.println(group.getSize());
		 canonicalMatrix(list,list2);
		 System.out.println(output.size());
		 
		 int[] w= new int[10];
		 w[0]=0;
		 w[1]=2;
		 w[2]=1;
		 w[3]=0;
		 w[4]=3;
		 w[5]=0;
		 w[6]=0;
		 w[7]=1;
		 w[8]=1;
		 w[9]=1;
		 
		 int[] arr= new int[4];
		 arr[0]=1;
		 arr[1]=1;
		 arr[2]=1;
		 arr[3]=2;
		 
		 
		 int[] v= new int[10];
		 v[0]=2;
		 v[1]=1;
		 v[2]=0;
		 v[3]=0;
		 v[4]=3;
		 v[5]=1;
		 v[6]=1;
		 v[7]=1;
		 v[8]=0;
		 v[9]=0;
		//System.out.println(blockBasedEqual(v,w,list));
		 //System.out.println(partition(list,1));
		 //System.out.println(LValue(list,4));
		 //canonicalMatrix2(mat,2);
		 
		 /**IAtomContainer ac = new AtomContainer();
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("C"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addAtom(new org.openscience.cdk.Atom("N"));
	     ac.addAtom(new org.openscience.cdk.Atom("H"));
	     ac.addBond(0, 2, Order.SINGLE);
	     ac.addBond(1, 3, Order.SINGLE);
	     ac.addBond(0, 4, Order.SINGLE);
	     ac.addBond(1, 5, Order.SINGLE);
	     ac.addBond(6, 7, Order.SINGLE);
	     AtomContainerDiscretePartitionRefiner refiner = PartitionRefinement.forAtoms().create();
		 PermutationGroup autoGroup = refiner.getAutomorphismGroup(ac);
		 for(Permutation perm: autoGroup.all()) {
			 System.out.println(perm.toCycleString());
		 }
		 autoGroup**/
		 /**PermutationGroupFunctions gen = null;
		 String[] args1= {"-f","C2H2O", "-v","-d", "C:\\Users\\mehme\\Desktop\\output.txt"};
		 try {
			 gen = new PermutationGroupFunctions();
			 gen.parseArgs(args1);
			 PermutationGroupFunctions.generate(PermutationGroupFunctions.molecularFormula, PermutationGroupFunctions.filedir);
		 } catch (Exception e) {
		 // We don't do anything here. Apache CLI will print a usage text.
			 if (PermutationGroupFunctions.verbose) e.getCause(); 
		 }**/
	 }
}
