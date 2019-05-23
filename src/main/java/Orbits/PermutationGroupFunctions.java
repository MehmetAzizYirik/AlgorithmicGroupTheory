package Orbits;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.interfaces.IAtomContainer;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class PermutationGroupFunctions {
	public static PermutationGroup group=PermutationGroup.makeSymN(4);
	public static int order=(int)group.order();
	public static int size=4;
	
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
		depiction.withCarbonSymbols().withSize(200, 200).withZoom(4).depict(molecule).writeTo(path);
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
	
	public static PermutationGroup generateGroup(ArrayList<Permutation> list) {
		List<Permutation> generators = new ArrayList<Permutation>();
        generators.addAll(list);
        return new PermutationGroup(size, generators);
	}
	
	/**
	 * Counts the number of graphs by the acting of the permutation group
	 * 
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
	 * Automorphism group of an integer set. Th groups are needed for the degree partition.
	 * But all automorphism groups of sets are already symmetry groups of the sets. But automorphism
	 * group of a group does not have to be equal to symmetry group. 
	 */
	public static void automorphismGroup(ArrayList<Integer> set) {
		
	}
	
	/**
	 * 
	 * Gluing Lemma
	 * 
	 */
	
	 public static void gluingLemma(ArrayList<Integer> Y) {
		 ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		 PermutationGroup groupX = PermutationGroup.makeSymN(subsets.size());
		 PermutationGroup groupY = PermutationGroup.makeSymN(Y.size());
		 ArrayList<ArrayList<Permutation>> productGroup = groupDirectProduct(groupX, groupY);
		 ArrayList<Integer> multiplicities=groupActionMaps(productGroup, subsets);
	 }
	 
	 /**
	  * For gluing lemma we need first the induced symmetry groups.
	  * If G is a symmetry group on a set X, the G induced is the 
	  * induced group by G on X. MOLGEN pg. 167
	  */
	 
	 public static void inducedPermutationGroup( PermutationGroup G, ArrayList<Integer> set) {
		 
	 }
	 
	 //TODO: Direct product can be a 2*n matrix 
	 
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
				 list.add(replace(s,perm));
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
	  * are considered as a multigraph. The definition of fixed map is given
	  * on page 81. 
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
}
