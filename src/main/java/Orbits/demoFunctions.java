/**
 * MIT License
 *
 * Copyright (c) 2018 Mehmet Aziz Yirik
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
 * This class is for the demo functions of the PeermutationGroupFunctions and other
 * classes from C.A.G.E.S. book. [1]
 * 
 * [1] Kreher, D.L. and Stinson, D.R., 1998. Combinatorial algorithms: generation,
 * enumeration, and search (Vol. 7). CRC press.
 * 
 * @author Mehmet Aziz Yirik
 */



package Orbits;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.group.PermutationGroup.Backtracker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

public class demoFunctions {
	
	/**
	 * 
	 * These functions were coded for the re-implementation of "Combinatorial
	 * Algorithms: Generation, Enumeration and Search" (C.A.G.E.S.) Book.
	 * 
	 */
	
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
	
	public static int[] multiply(int[] perm1, int[] perm2){
		int size=perm1.length;
		int[] newperm = new int[size];
		for(int i=0;i<size;i++) {
			newperm[i]=perm1[perm2[i]];
		}
		return newperm;
	}
	
	public static int[] identity(int size) {
		int[] identity= new int[size];
		for(int i=0;i<size;i++) {
			identity[i]=i;
		}
		return identity;
	}
	public static int[] powerOf(int[] perm,int n) {
		if(n==0) {
			return identity(perm.length);
		}else if(n==1) {
			return perm;
		}else {
			int[] power=identity(perm.length);
			while(n!=0) {
				power=multiply(power,perm);
				--n;
			}
			return power;
		}
	}
	
	public static int[] generator(int n) {
		int[] perm= new int[n];
		perm[n-1]=0;
		for(int i=0;i<n-1;i++) {
			perm[i]=i+1;
		}
		return perm;
	}
	
	public static int[] transposition(int n) {
		int[] perm= new int[n];
		perm[0]=1;
		perm[1]=0;
		for(int i=2;i<n;i++) {
			perm[i]=i;
		}
		return perm;
	}
	public static List<int[]> permutationGroup(int n){
		List<int[]> permutations = new ArrayList<int[]>();
		int[] transposition = transposition(n);
		int[] generator= generator(n);
		for(int i=0; i<n;i++) {
			int[] add= powerOf(generator,i);
			permutations.add(add);
			permutations.add(multiply(add,transposition));
		}
		return permutations;
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
	
	
	public static void backtrack(int l, Permutation g, Backtracker backtracker) {
		if (backtracker.isFinished()) { 
            return;
        }
        if (l == group.getSize()) {
            backtracker.applyTo(g);
        } else {
            for (int i = 0; i < group.getSize(); i++) {
                Permutation h = group.get(l,i);
                if (h != null) { 
                	backtrack(l + 1, g.multiply(h), backtracker);
                }
            }
        }
    }
	
	
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
                for (Permutation f : results) {
                    Permutation h = f.invert().multiply(p); 
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
	
	
}
