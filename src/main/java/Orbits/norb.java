package Orbits;

import java.util.Arrays;
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
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import Faulon.Atom;



public class norb {
	 	static Group grp;
		public static PermutationGroup group=PermutationGroup.makeSymN(4);
		//public static boolean DoneEarly;
		//public static int X;
		public static int size= group.getSize();
		//public static int[] T5= new int[size+1];
		//public static int[] T6= new int[size+1];
		//public static int[] T7= new int[size+1];
		//public static boolean[] T8 = new boolean[size+1];
		
		public static int choose(int n, int k) { 
			int C[][] = new int[n+1][k+1]; 
			int i, j;
			for (i = 0; i <= n; i++){ 
				for (j = 0; j <= Math.min(i, k); j++){ 
					if (j == 0 || j == i) {
						C[i][j] = 1; 
					}else {
						C[i][j] = C[i-1][j-1] + C[i-1][j];
					}
				}
			} 
			return C[n][k]; 
	    } 
		
		public static void type(Permutation g) {
			int i,j,ell;
			for(i=0;i<size;i++) {
				grp.T8[i]=true;
				grp.T7[i+1]=0;
			}
			for(i=0;i<size;i++) {
				if(grp.T8[i]) {
					ell=1;
					j=i;
					grp.T8[j]=false;
					while(grp.T8[g.get(j)]) {
						ell=ell+1;
						j=g.get(j);
						grp.T8[j]=false;
					}
					grp.T7[ell]=grp.T7[ell]+1;
				}
			}
		}
		
		public static void RecPartition(int k,int m,int B,int N){ 
			int[] c= new int[1000];
			int i,prod,n;
			n=size;
			if(m==0){
				for(i=1;i<=n;i++) {
					c[i]=0;
				}
				for(i=1;i<=N;i++) {
					c[grp.T6[i]]=c[grp.T6[i]]+1;
				}
				prod=1; 
				for(i=1;i<=n;i++) {
					prod=prod*choose(grp.T7[i],c[i]);
				}
				grp.T5[k]=grp.T5[k]+prod;
			}else{
				for(i=1;i<=Math.min(B,m);i=i+1){
					grp.T6[N+1]=i;
					System.out.println(k+" "+(m-i)+" "+i+" "+(N+1));
					RecPartition(k,m-i,i,N+1);
				}
			}
		}
		
		public static void norbBack() {
			grp.DoneEarly=false;
	        Backtracker norbUse = new Backtracker() {
	        	public void applyTo (Permutation g) {
	        		type(g);
	        		System.out.println("x"+" "+grp.X);
	        		for(int k=0;k<=grp.X;k++){
	        			grp.T6[1]=k;
	        			RecPartition(k,k,k,0);
	        		}
	        	}
	            public boolean isFinished() {
	                return grp.DoneEarly;
	            }
	        };
	        apply(norbUse);
	    }
		
		public static void apply(Backtracker backtracker) {
	        backtrack(0, grp.I, backtracker);
	    }
		
		public static void backtrack(int l, Permutation g, Backtracker backtracker) {
			if (grp.DoneEarly) { 
	            return;
	        }
	        if (l == size) {
	            backtracker.applyTo(g);
	        } else {
	            for (int i = 0; i < size; i++) {
	                Permutation h = group.get(l,i);
	                if (h != null) { 
	                	grp.pi.get(l).setTo(g.multiply(h));
	                	backtrack(l + 1, g.multiply(h), backtracker);
	                }
	            }
	        }
	    }
		
		public static void nOrb(PermutationGroup G, int[] N, int limit) {
			int k,deg,OrderG;
			deg=size;
			grp.DoneEarly=false;
			grp.X=limit;
			System.out.println("limit"+" "+limit);
			System.out.println("size"+" "+size);			
			for(k=0;k<=deg;k++) grp.T5[k]=0;
			norbBack();
			OrderG=Math.toIntExact(group.order()); 
			System.out.println("Order"+" "+OrderG);
			System.out.println(Arrays.toString(grp.T5));
			for(k=0;k<=deg;k++) N[k]=grp.T5[k]/OrderG;
			System.out.println(Arrays.toString(N));
		}
		
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
			norb.size= group.getSize();
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
		
		
		public static Set<Integer> act(Set<Integer> set, Permutation p) {
			Set<Integer> modifiedSet = new HashSet<Integer>();
			for(Integer element:set) {
				modifiedSet.add(p.get(element));
			}
			return modifiedSet;
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
		
		public static Set<Set<Integer>> ksubSet(int k,Set<Integer> set ) {
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
		
		public static void main(String[] args) {
			IAtomContainer acontainer = new org.openscience.cdk.AtomContainer();
			IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));	
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));
			acontainer.addAtom(builder.newInstance(IAtom.class, "C"));
			
			
			acontainer.getAtom(0).setImplicitHydrogenCount(1);
			acontainer.getAtom(1).setImplicitHydrogenCount(1);
			acontainer.getAtom(2).setImplicitHydrogenCount(2);
			acontainer.getAtom(3).setImplicitHydrogenCount(3);
			acontainer.getAtom(4).setImplicitHydrogenCount(1);
			acontainer.getAtom(5).setImplicitHydrogenCount(1);
			MoleculeSignature moleculeSignature = new MoleculeSignature(acontainer);
			String canonicalSignature = moleculeSignature.toCanonicalString();
			List<org.openscience.cdk.signature.Orbit> orbits = moleculeSignature.calculateOrbits();
			System.out.println(orbits);
			
			
			/**PermutationGroup groupp= PermutationGroup.makeSymN(4);
			System.out.println(group.getSize());
			norb.group=groupp;
			grp= new Group();
			int[] N= new int[size+1];
			
			nOrb(groupp,N,4);**/
		}
	}
