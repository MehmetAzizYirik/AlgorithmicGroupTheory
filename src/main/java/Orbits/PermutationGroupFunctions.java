package Orbits;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import org.openscience.cdk.group.PermutationGroup.Backtracker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class PermutationGroupFunctions {
	public static PermutationGroup group=PermutationGroup.makeSymN(4);
	public static int order=(int)group.order();
	public static int size=4;
	class group_data{
		Permutation I;
		Permutation T0;
		Permutation T1;
		Permutation T2;
		Permutation T3;
		Permutation C1;
		Permutation C2;
		Permutation T4;
		int T5;
		int T6;
		int T7;
		int T8;
		Permutation pi;
		int DoneEarly;
		int X;
	}
	
	class orb_data{
		Set S;
		Set S2;
		setlist C= new setlist();
		int X;
		int OptX;
		Permutation beta;
		Permutation f;
		Permutation Invf;
		Permutation F;
	}
	
	class set_data{
		Set FullSet;
		Set EmptySet;
	}
	
	class setlist{
		universe U = new universe();
		int size ;
		Set blocks;
	}
	
	class iso_data{
		setlist TQ = new setlist();
		setlist TP0 = new setlist(); 
		setlist TP1 = new setlist();
		setlist TP2 = new setlist();
		setlist TP3 = new setlist();
		setlist TR = new setlist();
		Permutation best;
		int BestExist;
		Set Tx;
		Set Ty;
		setlist L= new setlist();
		Set U;
		setlist S= new setlist();
		Set T;
		int N;
		Permutation f;
		Permutation Invf;
	}
	
	class universe{
		int NODES;
		int WORDS;
		int n;
		group_data GRP= new group_data();
		orb_data ORB = new orb_data(); 
		set_data SET = new set_data(); 
		iso_data ISO = new iso_data();
	}
		
	/**
	 * These functions for the lexicographical generation of subsets.
	 * The functions are explained in CAGES books chapter 2. 
	 */
	
	public static int presence(Set<Integer> subset, int element, int size) {
		int presence=0;
		if(subset.contains(size-element)) {
			presence=presence+1;
		}
		return presence;
	}
	
	public static Set<Integer> buildSet(int size){
		Set<Integer> set = new HashSet<Integer>();
		for(int i=0;i<size;i++) {
			set.add(i);
		}
		return set;
	}
	
	public static int[] characteristicVector(Set<Integer> subset, int size) {
		int[] vector = new int[size];
		Set<Integer> set=buildSet(size);
		int size2=size-1;
		for(Integer element:set) {
			vector[size2-element]=presence(subset,element,size2);
		}
		return vector;
	}
	
	public static int power(int a, int b) {
		return (int)Math.pow(a, b);
	}
	
	public static void setSize(PermutationGroup group) {
		Permutate.size= group.getSize();
	}
	
	public static int getLexRank(Set<Integer> subset) {
		int rank=0;
		int [] vector=characteristicVector(subset,size);
		int size2=size-1;
		for(int i=size2;i>=0;i--) {
			rank+=vector[i]*power(2,size2-i);
		}
		return rank;
	}
	
	public static Set<Integer> getLexUnrank(int rank) {
		Set<Integer> set = new HashSet<Integer>();
		for(int i=size;i>=1;i--) {
			if(rank%2==1) {
				set.add(i-1);
				rank=rank/2;
			}
		}
		return set;
	}
	
	//TODO: This orbitSet can be also build in an order. So everytime, just check the first set,
	//which has the lowest rank in the set of orbits.
	
	public static Set<Integer> predecessorSet(List<Set<Integer>> orbitSet, Set<Integer> orbit) {
		Set<Integer> predecessor= new HashSet<Integer>();
		int rank1 = getLexRank(orbit); // Rank of the input orbit.
		for(Set<Integer> set: orbitSet) {
			if(getLexRank(set)<rank1) {
				predecessor=set;
			}
		}
		return predecessor;
	}
	
	public static Set<Integer> maximumSet(List<Set<Integer>> orbitSet){
		return orbitSet.get(orbitSet.size());
	}
	
	/**
	 * First, learn how to rebuild the backtracker like Gill did.
	 * So the following backtrack functions are from Gill.
	 */
	public static void backtrack(int l, Permutation g, Backtracker backtracker) {
		if (backtracker.isFinished()) { // He coded this backtrack in general. If needed finish is checked.
            return;
        }
        if (l == group.getSize()) {
            backtracker.applyTo(g);
        } else {
            for (int i = 0; i < group.getSize(); i++) {
                Permutation h = group.get(l,i);
                if (h != null) { // Because he defined them as null. Need to understand sims chain
                	//In the book it says mult but the mult is done in backtrack function g.multiply(h)
                	backtrack(l + 1, g.multiply(h), backtracker);
                }
            }
        }
    }
	public static ArrayList<Integer> act(ArrayList<Integer> set, Permutation p) {
		ArrayList<Integer> modifiedSet = new ArrayList<Integer>();
		for(Integer element:set) {
			modifiedSet.add(p.get(element));
		}
		return modifiedSet;
	}
	
	//TODO: Continue from this backtrack function to build orbitbacktrack function.
	//TODO: Understand what R and S sets are in the book.
	/**public static void orbitBacktrack(int l, Permutation g, Backtracker backtracker) {
        if (backtracker.isFinished()) { 
            return;
        }
        if (l == group.getSize()) {
            Set<Integer> A= maximumSet(S);
        }
        while(A.size()!=0) {
        	Set<Integer> C= predecessorSet(S,A);
        	if(getLexRank(act(A,g))<getLexRank(A)) {
        		S.delete(A);
        		if(!S.contains(A)) {
        			S.add(A);
        		}
        		if(S.size()==group.getSize()) {
        			finished=true;
        			return;
        		}
        	A=C;
        	} else {
                for (int i = 0; i < group.getSize(); i++) {
                    Permutation h = group.get(l,i);
                    if (h != null) { // Because he defined them as null. Need to understand sims chain
                    	//In the book it says mult but the mult is done in backtrack function g.multiply(h)
                        orbitBacktrack(l + 1, g.multiply(h), backtracker);
                    }
                }
            }
        }
    }**/
	
	/**
     * Apply the backtracker to all permutations in the larger group.
     *
     * @param backtracker a hook for acting on the permutations
     */
    public static void apply(Backtracker backtracker) {
        backtrack(0, new Permutation(group.getSize()), backtracker);
    }
	
    
    /**
     * Generate a transversal of a subgroup in this group.
     *
     * @param subgroup the subgroup to use for the transversal
     * @return a list of permutations
     */
    
    /**
     * Schreier-Sims Algorithm returns the group element representation of the orbits
     * Uo means the elements of the group sending 0 to the orbit element.
     * And all the Us in the group are the left transversal of their Gs.
     *
     **/
    public static List<Permutation> transversal(final PermutationGroup subgroup) {
        final long m = group.order() / subgroup.order();
        final List<Permutation> results = new ArrayList<Permutation>();
        Backtracker transversalBacktracker = new Backtracker() {
            private boolean finished = false;

            public void applyTo(Permutation p) {
            	//System.out.println("applyto"+" "+p);
                for (Permutation f : results) {
                    Permutation h = f.invert().multiply(p); // why we have invert function and what is the use() exactly
                    //System.out.println("apply to invert function"+" "+h);
                    //TODO: What is this test case ? Is subgroup test just for transversal ? Or for all the backtrackers ?
                    if (subgroup.test(h) == group.getSize()) {
                        return;
                    }
                }
                results.add(p);
                if (results.size() >= m) {
                    this.finished = true;
                }
            }

            public boolean isFinished() {
                return finished;
            }

        };
        apply(transversalBacktracker);
        return results;
    }
    public static List<Permutation> all() {
        final List<Permutation> permutations = new ArrayList<Permutation>();
        Backtracker counter = new Backtracker() {

            public void applyTo(Permutation p) {
                permutations.add(p);
            }

            public boolean isFinished() {
                return false;
            }
        };
        apply(counter);
        return permutations;
    }
    
    /**
     * Apply the backtracker to all permutations in the larger group.
     *
     * @param backtracker a hook for acting on the permutations
     */
    public static void apply2(Backtracker backtracker) {
        //orbitBacktrack(0, new Permutation(group.getSize()), backtracker);
    }
    
	public static List<Permutation> orbitRepresentatives() {
		final List<Permutation> representations = new ArrayList<Permutation>();
		Backtracker orbitRepresentation = new Backtracker() {

            public void applyTo(Permutation p) {
                representations.add(p);
            }

            public boolean isFinished() {
                return false;
            }
        };
        
        apply(orbitRepresentation);
        return representations;
	}
	
	/***********************************************************/
	public static IAtomContainer getAtomContainer(String elements, int[][] bonds) {
        IChemObjectBuilder builder =  SilentChemObjectBuilder.getInstance();
        IAtomContainer ac = builder.newInstance(IAtomContainer.class); 
        for (int i = 0; i < elements.length(); i++) {
            String eS = String.valueOf(elements.charAt(i));
            ac.addAtom(builder.newInstance(IAtom.class, eS));
        }
        for (int[] bond : bonds) {
            int o = bond[2];
            ac.addBond(bond[0], bond[1], IBond.Order.values()[o - 1]);
        }
        return ac;
    }
	
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
	/*
	 * To have a qunique version, I put ArrayList<Integer> rather than Set.
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
	 * Gneerating the k-disjoint subsets.
	 * @param k
	 * @param set
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
	public static ArrayList<ArrayList<Integer>> cloneList(ArrayList<ArrayList<Integer>> list){
		ArrayList<ArrayList<Integer>> result= new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> i:list) {
			result.add(i);
		}
		return result;
	}
	
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
	
	public static int checkStab(Permutation g, ArrayList<Integer> set){
		int no=0;
        if(set.equals(act(set,g))) {
			no++;
		}
		return no;
	}
	
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
	
	public static int combination(int m, int n) {
		return factorial(m) / (factorial(n) * factorial(m - n));
	}

	public static int factorial(int i){
		if (i==0){
			return 1;
		}
		return i * factorial(i - 1);
	}
	
	public static void recPartition(int k, int m, int B, int N) {
		int n=100;
		int[] c= new int[1000];
		int[] T= new int[size];
		int[] V2= new int[size];
		int[] V3= new int[size];
		int prod;
		if(m==0) {
			for(int i=1;i<=n;i++) {
				c[i]=0;
			}
			for(int i=1;i<=N;i++) {
				c[V2[i]]=c[V2[i]]+1;
			}
			prod=1;
			for(int i=1;i<=n;i++) {
				prod=prod*combination(T[i],c[i]);
			}
			V3[k]=V3[k]+prod;
		}else {
			for(int i=1;i<Math.min(B, m);i++) {
				V2[N+1]=i;
				recPartition(k,m-i,i,N+1);
			}
		}
		
	}
	/**public static void norbUse(int n,Permutation g){
		//TODO: You need to define the type as V1 but as global. Then also V2 and V3 as global. So all the functions generate this arrays.
		int OrderG,k,deg;
		int[] V1, V2;
		
		deg=G->U->n;
		G->U->GRP->DoneEarly=false;
		G->U->GRP->X=limit;
		for(k=0;k<=deg;k++) G->U->GRP->T5[k]=0;
		Run(stdout,G,NorbUse);
		OrderG=GroupOrder(G);
		for(k=0;k<=deg;k++) N[k]=G->U->GRP->T5[k]/OrderG;
		
		V1=type(n,g);
		for(int k=0;k<=n;k++) {
			V2[1]=k;
			recPartition(n,k,k,k,0);
		}
	}**/
	
	/**public static void Run( PermutationGroup G, void (* Use)(FILE *,perm) )
	{
	  G->U->GRP->DoneEarly=false;
	  RunBacktrack(F,0,G,G->U->GRP->I,Use);
	}**/
	
	
	public static int[] N;
	/**public static int[] norb(int n) {
		doneEarly=false;
		for(int k=0;k<n;k++) {
			V3[k]=0;
		}
		for(Permutation p:group.all()) {
			norbUse(n,p);
		}
		N = new int[n];
		for(int j=0;j<n;j++) {
			N[j]=V3[j]/size;
		}
		return N;
	}**/
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
	/**public static boolean minCheck(ArrayList<ArrayList<Integer>> list, ArrayList<Integer> orbit) {
		boolean check=true;
		for(ArrayList<Integer> set:list) {
			
		}
	}**/
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
	
	/**
	 * 
	 * 
	 **/
	public static ArrayList<ArrayList<Integer>> getOrbits2(PermutationGroup group, ArrayList<ArrayList<Integer>> set){
		ArrayList<ArrayList<Integer>> orbits = new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> s:set) {
			orbits.addAll(getOrbits(group,s));
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
	public static ArrayList<Integer> getBase(PermutationGroup group) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		for(int i=0;i<group.getSize();i++) {
			base.add(i);
		}
		return base;
	}
	
	/**public static int numberOfOrbits(PermutationGroup group) {
		int count=0;
		Set<Set<Integer>> subsets= subSet(getBase(group));
		for(Set<Integer> set:subsets) {
			//System.out.println(getOrbits(group,set).size()+" "+getStabilizers(group,set).size());
			count+=(getStabilizers(group,set).size()/group.all().size());
		}
		return count;
	}**/
	
	public static int numberOfOrbit(PermutationGroup group, ArrayList<ArrayList<Integer>> subset) {
		int count=0;
		for(Permutation perm:group.all()) {
        	for(ArrayList<Integer> set:subset) {
        		count=count+checkStab(perm,set);
        		if(checkStab(perm,set)!=0) {
        			System.out.println(perm.toCycleString()+" "+set); // Why ?
        		}
        	}
        }
		return count/group.all().size();
	}
	public static int numberOfOrbits2(int k) {
		int count=0;
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(k,getBase(group));
		for(Permutation perm:group.all()) {
			//System.out.println(getOrbits(group,set).size()+" "+getStabilizers(group,set).size());
			count+=getStabilizers2(perm,subsets);
		}
		return count/group.all().size();
	}
		
	/**
	 * Molecule depiction generator
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
	
	public static Multimap<HashMap<Integer, Integer>, Permutation> conjugacyClasses(PermutationGroup group) {
		Multimap<HashMap<Integer, Integer>, Permutation> multiMap = ArrayListMultimap.create();
		HashMap<Permutation, HashMap<Integer, Integer>> map=cycleType(group);
		for(Permutation entry:map.keySet()) {
			multiMap.put(map.get(entry),entry);
		}
		return multiMap;
	}
	
	/**
	 * 
	 * Cauchy Frobenious Theorem - Number of Orbits
	 *
	 **/
	
	public static PermutationGroup generateGroup(Permutation perm) {
		List<Permutation> generators = new ArrayList<Permutation>();
        generators.add(perm);
        return new PermutationGroup(size, generators);
	}
	
	public static PermutationGroup generateGroup(ArrayList<Permutation> list) {
		List<Permutation> generators = new ArrayList<Permutation>();
        generators.addAll(list);
        return new PermutationGroup(size, generators);
	}
	
	public static int numberOfGraphs(PermutationGroup group, int m) {
		int count=0;
		ArrayList<ArrayList<Integer>> subsets= ksubSet2(2,getBase(group));
		Multimap<HashMap<Integer, Integer>,Permutation> classes= conjugacyClasses(group);
		for(HashMap<Integer, Integer> key:classes.keySet()) {
			Permutation perm = (Permutation) classes.get(key).toArray()[0];
			System.out.println(perm.toCycleString());
			PermutationGroup group2=generateGroup(perm);
			//System.out.println(perm.toCycleString()+" "+numberOfOrbit(group2,subsets)+" "+group2.order());
			count=count+classes.get(key).size()*power(m,numberOfOrbit(group2,subsets));
		}
		return count/group.all().size();
	}
	
	/**
	 * Polya's enumeration method:
	 * 
	 * Group reduction method  for which I need cycle index polynomial first
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
	 **/
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
	
	public static int numberOfLengthOrbits(Permutation perm, Set<Set<Integer>> set,int i) {
		PermutationGroup group= generateGroup(perm);
		return getLengthBasedOrbits(group,set,i).size();
	}
	
	public static int[] buildZerosArray(int n) {
		int arr[] = new int[n];
		for(int i=0;i<arr.length;i++) {
		    arr[i] = 0;
		}
		return arr;
	}
	
	public static int[] normalizeArray(int[] arr) {
		for(int i=1;i<arr.length;i++) {
			arr[i]=arr[i]/i;
		}
		return arr;
	}
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
	 * Permutation Types Based on MOLGEN Book.
	 */
	
	public static ArrayList<Integer> cyclicFactors(Permutation perm) {
		int[] cyclic=type(5,perm);
		System.out.println("type"+" "+Arrays.toString(cyclic));
		ArrayList<Integer> lengths= new ArrayList<Integer>();
		for(int i=1;i<cyclic.length;i++) {
			if(cyclic[i]>0) {
				lengths.add(cyclic[i]);
			}
		}
		Collections.sort(lengths,DES_ORDER);
		return lengths;
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
	public static HashMap<Integer,Integer> nTuple(Permutation perm) {
		ArrayList<Integer> list= cyclicFactors(perm);
		Collections.sort(list, DES_ORDER);
		System.out.println("cyclic"+" "+list);
		int n=sum(list);
		System.out.println("n"+" "+n);
		HashMap<Integer,Integer> map= new HashMap<Integer,Integer>();
		for(int i=1;i<=n;i++) {
			map.put(i,count(list,i));
		}
		removeLastZeros(map,n);
		return map;
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
	
	public static final Comparator<Integer> DES_ORDER = new Comparator<Integer>() {
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
	
	/**
	 * Partition with ArrayList
	 * @param n
	 * @param d
	 * @param depth
	 * @return
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
	
	public static List<ArrayList<Integer>> buildArray2(int n,int d, int depth){
		List<ArrayList<Integer>> array= new ArrayList<ArrayList<Integer>>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(ArrayList<Integer> item: partition2(n-i,d,depth+1)) {
					item=addElement2(item,i);
			        if(item.size()==d) {
			        	if(sum(item)==n) { //TODO: BE CAREFUL!
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
			        	if(sum(item)==3) { //TODO: BE CAREFUL!
			        		array.add(item);
			        	}
			        }else {
			        	array.add(item);
			        }
			}
		}
		return array;
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
	 * 
	 * @param power power on the summation
	 * @param terms: How many variables you sum
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
	
	public static HashMap<ArrayList<Integer>,Integer> monomial(int var, int inPower){
		HashMap<ArrayList<Integer>,Integer> map= new HashMap<ArrayList<Integer>,Integer>();
		for(int i=0;i<var;i++) {
			map.put(multiplyList(simpleList(var,i),inPower), 1);
		}
		return map;
	}
	
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
	/**
	 * Group Reduction Function. The formula given in MOLGEN Book pg 77.
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
	 * Coset and Double Cosets on Groups - MOLGEN 1.36
	 * @param args
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
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
	//TODO: Theoretically check. Need to be tested.
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
	 * Group stabilizers are list of permutations. That is why we have this version.
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
	 * Group action for direct product groups
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
	 * Double Cosets: Orbits for the group build by the group action of a product group
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
	 */
	public static ArrayList<Permutation> doubleCoset(ArrayList<ArrayList<Permutation>> product, Permutation perm) {
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
	 * G/Gx Gx is the stabilizer
	 */
	public static ArrayList<Permutation> stabilizerOrbits(PermutationGroup group, PermutationGroup subGroup, ArrayList<Integer> set) {
		ArrayList<Permutation> stabilizers= getStabilizers(group,set);
		System.out.println("stabilizers"+" "+stabilizers);
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
	
	public static ArrayList<ArrayList<Permutation>> stabilizerOrbitsList(PermutationGroup group, PermutationGroup subGroup, ArrayList<ArrayList<Integer>> sets) {
		ArrayList<ArrayList<Permutation>> result= new ArrayList<ArrayList<Permutation>>();
		for(ArrayList<Integer> set: sets) {
			System.out.println("set"+" "+set);
			result.add(stabilizerOrbits(group,subGroup,set));
		}
		return result;
	}
	
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
	}
	
	/**
	 * MOLGEN - minimal orbit representation
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
	 * Counting the frequency of integers in the lists.
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
	
	/**
	 * Return the source index for the bond adding.
	 **/
	
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
	 */
	
	public static ArrayList<ArrayList<Integer>> minBondAdd(ArrayList<ArrayList<Integer>> list) {
		ArrayList<Integer> pair= new ArrayList<Integer>();
		pair.add(minSourceIndex(list));
		pair.add(minTargetIndex(list));
		list.add(pair);
		return list;
	}
	
	public static int[] OptX = new int[size];
	public static ArrayList<ArrayList<Integer>> C= new ArrayList<ArrayList<Integer>>();
	
	public static int FindLargestElement(ArrayList<Integer> set) {
		int result= set.get(0);
		for(int i=0;i<set.size();i++) {
			if(set.get(i)>result) result=set.get(i);
		}
		return result;
	}
	/**public static void Apply(Permutation g,ArrayList<Integer> A,ArrayList<Integer> B){
		int u,n,h,z;
		n=size;
		ArrayList<Integer> set = new ArrayList<Integer>();
		h=0;
		u=FindLargestElement(A);
		z=SetOrder(A);
		while(h<z){
			h=h+1;
			SetInsert(g->V[u],g->U->ORB->S2);
			if(h<z) u=FindPrevElement(u,A);
		}
		GetSet(B,g->U->ORB->S2);
	}**/
	
	public static void MinRep(PermutationGroup G, ArrayList<Integer> A) {
		int deg=size;
		int k=0; 
		for(int i=0;i<deg;i++) { 
			if(A.contains(i)) {
				OptX[k++]=i;
			}
		}
		C.add(A);
		MinRepBT(k,G,0);
		ArrayList<Integer> ans= new ArrayList<Integer>();
		for(int i=0;i<k;i++) {
			ans.add(OptX[i]);
		}
	}
	
	public static void SetPerm(int k,int[] A,Permutation g){
		int i,j,n;
		n=size;
		int[] V=g.getValues();
		//GetEmptySet(ORB->S);
		ArrayList<Integer> set= new ArrayList<Integer>();
		for(i=0;i<k;i++){
			set.add(A[i]);
			V[i]=A[i];
		}
		for(j=0;j<n;j++) { 
			if(!set.contains(j)){
				V[i]=j;
				i=i+1;
			}
		}
        Permutation p = new Permutation(V);
        g=p;
	}
	public static int[] X= new int[size];
	public static Permutation beta= new Permutation(size);
	/**public static void MinRepBT(int k,PermutationGroup G,int ell){
		int i,j,x,m,r,start,n;
		Permutation g = null;
		n=size;
		if (ell==0){
			start=0;
		}else{
			start=X[ell-1];
		}

		m=n;
		for(x=start;x<n;x++) { 
			if(C.get(ell).contains(x)){
				r=0;
				while((r<ell)&&(X[r]==OptX[r])) r=r+1;
				if( (r<ell) && X[r]>=OptX[r]) return;
				X[ell]=x; 
				SetPerm(ell+1,X,beta);
				G.changeBase(beta);

				for(j=0;j<=x;j++)
					if((g=G.get(ell,j))!=null) break;
					if (g.getValues()[x]<=m) {
						m=g.getValues()[x];
						X[ell]=m;
	 
						if ( k==ell+1) {
							for(i=r;(i<k) && (X[i]==OptX[i]);i++) {
								if( (i!=k) && X[i]<OptX[i]) {
									for(i=0;i<k;i++) OptX[i]=X[i];	
								}
							}
						}else{
							Apply(g,C.get(ell),C.get(ell+1));
							SetDelete(m,C->blocks[ell+1]);
							MinRepBT(k,G,ell+1);
						}
					}
			}
		}
	}**/
	
	/**
	 * MOLGEN: Homomorphism principle
	 * 
	 * Mapping m multigraph to (m-1) multigraph
	 */
	
	public static int bond(ArrayList<Integer> pair) {
		return 1;
	}
	
	public static int[] bondMultiplicity(int n) {
		ArrayList<ArrayList<Integer>> sets=ksubSet2(2,getBase(group));
		int[] mult= new int[sets.size()];
		for(int i=0;i<sets.size();i++) {
			mult[i]=0;
		}
		return mult;
	}
	
	public static int homomorphism(ArrayList<Integer> pair, int m) {
		if(bond(pair)<(m-1)) {
			return bond(pair);
		}else {
			return (m-2); // (m-2) is a graph or the integer ?
		}
	}
	/**
	 * In the bond list of inverse homomorphism, increasing the bonds.
	 * @param bond
	 * @param index
	 */
	
	public static void bondIncrease(int[] bond, int index) {
		if(index!=0) {
			for(int i=0;i<=(index-1);i++) {
					bond[bond.length-1-i]=bond[bond.length-1-i]+1;
			}	
		}
	}
	
	/**
	 * From molgen book, the inverse homomorphism for m-multiplicity. 
	 * @param bond
	 * @return
	 */
	
	public static int[] cloneArray(int[] arr) {
		int[] result= new int[arr.length];
		for(int i=0;i<arr.length;i++) {
			result[i]=arr[i];
		}
		return result;
	}
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
	 * 
	 * Generation of m-multigraphs. For a given set of numbers increasing them 
	 * in an order until they reach the upper limit, m.
	 * 
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
	/**
	 * Young Group construction
	 * @param set
	 * @param group
	 */
	
	/**
	 * Build the Young subgroup for the given disjoint set
	 * @param set
	 * @param group
	 * @return
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
	 * @param sets
	 * @param group
	 * @return
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
	 */
	
	public static ArrayList<ArrayList<Permutation>> youngSubgroupRepresentation(int k, PermutationGroup group) {
		ArrayList<ArrayList<ArrayList<Integer>>> kDisjointSets= disjoint(k,getBase(group));
		return youngSubgroups(kDisjointSets,group);
	}
	
	//TODO: Look for better representation of product group construction rather than the list.
	//TODO: This group representation also for the stab group and young subgroups.
	public static void main(String[] args) throws CloneNotSupportedException, CDKException, IOException {   
		ArrayList<ArrayList<Integer>> subSets= ksubSet2(2,getBase(group));
		//System.out.println(subSets);
		
		int[] arr= new int[3];
		arr[0]=0;
		arr[1]=0;
		arr[2]=0;
		
		int[] arr2= new int[3];
		arr2[0]=0;
		arr2[1]=0;
		arr2[2]=0;
		ArrayList<Integer> list= new ArrayList<Integer>();
		list.add(0);
		list.add(1);
		list.add(2);
		
		//disjoint(2,list);
		disjoint(2,list);
		for(ArrayList<ArrayList<Integer>> l:disjoints) {
			//System.out.println(l);
		}
		//System.out.println(Arrays.equals(arr,arr2));
		//increaseBond(arr,3);
		
		//System.out.println("--------------------------------------");
		//printListArray(inc);
		//printListArray(inverseHomomorphism(bond));
				
		List<Permutation> generators = new ArrayList<Permutation>();

        // p1 is (0, 1)
		int size=4;
        int[] p1 = new int[size];
        p1[0] = 1;
        p1[1] = 0;
        for (int i = 2; i < size; i++) {
            p1[i] = i;
        }

        // p2 is (1, 2, ...., n, 0)
        int[] p2 = new int[size];
        p2[0] = 1;
        for (int i = 1; i < size - 1; i++) {
            p2[i] = i + 1;
        }
        p2[size - 1] = 0;

        generators.add(new Permutation(p1));
        generators.add(new Permutation(p2));
		PermutationGroup group = PermutationGroup.makeSymN(size);
		//System.out.println(group);
		//System.out.println(group.getLeftTransversal(0));
		//System.out.println(group.getLeftTransversal(1));
		//System.out.println(group.getLeftTransversal(2));
		//System.out.println(group.getLeftTransversal(3));
		//System.out.println(group);
		Permutation perm = new Permutation(1,0,2,3);
		PermutationGroup h= generateGroup(perm);
		//System.out.println(group.transversal(h));
		/**ArrayList<ArrayList<Integer>> subSets= ksubSet2(2,getBase(group));
		for(ArrayList<Integer> set: subSets) {
			ArrayList<ArrayList<Integer>> orbit=getOrbits(group, set);
			System.out.println(orbit);
		}
		
		
		int size = 4;
	    // Sym(n) : make the total symmetry group
	    PermutationGroup group = PermutationGroup.makeSymN(size);
	    
	    for(int i=0;i<size;i++) {
			System.out.println(group.getLeftTransversal(i));
		}
	    
	    // Aut(G) : make the automorphism group for a graph
	    Permutation p1 = new Permutation(2, 1, 0, 3);
	    Permutation p2 = new Permutation(0, 3, 2, 1);
	    List<Permutation> generators = new ArrayList<Permutation>();
	    generators.add(p1);
	    generators.add(p2);
	    PermutationGroup subgroup = new PermutationGroup(size, generators);

	    // generate the traversal
	    List<Permutation> transversal = group.transversal(subgroup);

	    int subgroupOrder = (int) subgroup.order();
	    int groupOrder = (int) group.order();
	    int transversalSize = transversal.size();
	   
	    // check that |Aut(G)| / |Sym(N)| = |Transversal|
	    System.out.println(factorial(size)+" "+groupOrder);
	    System.out.println((groupOrder / subgroupOrder)+" "+transversalSize);
	    System.out.println(subgroupOrder);
		for(int i=0;i<4;i++) {
			System.out.println(group.getLeftTransversal(i));
		}**/
	    
	    //ArrayList<ArrayList<Integer>> subSets= ksubSet2(2,getBase(group));
		//constructTrans(group,group2);
		
	    //stabilizerOrbitsList(group,group2, subSets);
	}
}
