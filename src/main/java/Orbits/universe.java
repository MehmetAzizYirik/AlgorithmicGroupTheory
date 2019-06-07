package Orbits;

import java.util.HashSet;
import java.util.Set;

public class universe {
	public static int n=4;
	public static int WORDSIZE=32;
	public static int WORDS=ComputeWords(n);

	HashSet<Integer> SET = SetInit(U);
	GRP = GroupInit(U);
	ORB = OrbInit(U);
	public static int ComputeWords(int n) {
		return (n/WORDSIZE+1);
	}
	 
	public static Set<Integer> SetInit(){
		int i;
		int u;
		int j;
		setData set = new setData();
		for(u=0;u<n;u++){
			j=WORDSIZE-1-(u%WORDSIZE);
			i=u/WORDSIZE;
			set.FullSet->V[i]=SET->FullSet->V[i]|((UINT)1<<j);
		}
		for(i=0;i<WORDS;i++)  SET->EmptySet->V[i]=0;
		return (SET);
	}
}
