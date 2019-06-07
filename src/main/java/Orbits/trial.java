package Orbits;
import java.util.*;

public class trial
{
	//
	// The goal of the program is to verify Burnside's lemma.
	// Burnside's lemma is useful for counting things with symmetry.
	//
	// This program uses the lemma to compute the number of necklaces
	// that can be formed by n = 5 beads with c = 10 colors such that
	// all reflections and rotations of the same coloring is counted 
	// only once.
	//
	// Here is the key ingredient in the proof for the Burnside's lemma
	// 
	// y in Orbit of x means y = g x, define a map from X to 2^G as follow:
	//
	// Map an element y to the set of group element g such that y = g x
	// y = h x => y = g g^{-1} h x = g (g^{-1} h) x, 
	// Notice g^{-1} h x = g^{-1} g x = x, therefore for every group element 
	// y maps to, it is a coset element.
	// On the other hand, every coset element is valid to be mapped to by y, 
	// Therefore we can map injectively an element in the orbit to a coset
	// Therefore, the size of the orbit is the number of cosets = |G|/|Gx|
	//
	// Last but not least, if for each element, we add up the size its stablizer group, 
	// we get this:
	//
	//   sum for each x, sum |Gx|
	// = for each orbit, for each orbit element, sum |Gx|
	// = for each orbit, sum |G|/|Gx| * |Gx|
	// = for each orbit, sum |G|
	// = number of orbit * |G|
	// 
	// So this is the Burnside's lemma, you can calculate the number of orbits
	// by summing the size of all the stablizer group and then divide it by the
	// size of the group.
	//
	// Note that when we compute sum for each x, sum |Gx|, it is equivalent to 
	// flip to loop and ask, for each permutation, sum the number of configuration it fixes.
	//
	// This is where the Polya enumeration formula comes in. After we factor the 
	// permutation into disjoint cycles, we figure that in order for a configuration is to be fixed
	// It must have the same color for all nodes within a cycle, therefore, 
	// for each cycle, we can make c choices for the color, and therefore we can count the number
	// of configuration get fixed by a permutation without enumerating the configurations at all.
	//

	public static void main(String[] args)
	{
		ArrayList<ArrayList<Integer>> group_elements = new ArrayList<ArrayList<Integer>>();
		int n = 4;
		int c = 4; // Polya's enumeration work much better when c is large - enumerating the configuration space is expensive
		for (int offset = 0; offset < n; offset++)
		{
			ArrayList<Integer> forward = new ArrayList<Integer>();
			for (int i = 0; i < n; i++)
			{
				forward.add((i + offset) % n);
			}
			group_elements.add(forward);
			ArrayList<Integer> backward = new ArrayList<Integer>();
			for (int i = 0; i < n; i++)
			{
				backward.add(((n - i) + offset) % n);
			}
			group_elements.add(backward);
		}

		boolean burnside = true;
		boolean polya = true;

		System.out.print("Begin Group");
		System.out.print("\n");
		for (int i = 0; i < group_elements.size(); i++)
		{
			for (int j = 0; j < group_elements.get(i).size(); j++)
			{
				System.out.print(group_elements.get(i).get(j));
				System.out.print(" ");
			}
			System.out.print("\n");
		}
		System.out.print("End Group");
		System.out.print("\n");
		if (burnside)
		{
			System.out.print("Starting brute force verification of Burnside's lemma");
			System.out.print("\n");
			// This version of the code directly exercise the Burnside's lemma as is:
			// This is inefficient because we have to go through all configurations
			ArrayList<Integer> config = new ArrayList<Integer>();
			ArrayList<Integer> mirror = new ArrayList<Integer>();
			VectorHelper.resize(config, n);
			VectorHelper.resize(mirror, n);

			int all_stablizer_count = 0;
			for (int k = 0; k < group_elements.size(); k++)
			{
				for (int j = 0; j < group_elements.get(k).size(); j++)
				{
					System.out.print(group_elements.get(k).get(j));
					System.out.print(" ");
				}
				int fixed_count = 0;
				int seq = 0;

				// A simple odometer to loop through all the config
				while (true)
				{
					int cur = seq;
					boolean last = true;
					for (int i = 0; i < n; i++)
					{
						config.set(i, cur % c);
						cur = cur / c;
						last = last & (config.get(i) == c - 1);
					}

					// Here we have got a configuration - check if it is fixed by the current permutation

					boolean is_fixed = true;
					for (int l = 0; is_fixed && l < n; l++)
					{
						is_fixed = is_fixed && (config.get(group_elements.get(k).get(l)) == config.get(l));
					}

					if (is_fixed)
					{
						fixed_count++;
					}

					if (last)
					{
						break;
					}
					seq++;
				}
				System.out.print(" fixes ");
				System.out.print(fixed_count);
				System.out.print(" configurations");
				System.out.print("\n");
				all_stablizer_count += fixed_count;
			}
			System.out.print(all_stablizer_count);
			System.out.print(" ");
			System.out.print(group_elements.size());
			System.out.print("\n");
		}

		if (!polya)
		{
			System.out.print("Starting polya's enumeration formula");
			System.out.print("\n");
			int all_stablizer_count = 0;
			for (int i = 0; i < group_elements.size(); i++)
			{
				ArrayList<Boolean> used = new ArrayList<Boolean>();
				VectorHelper.resize(used, group_elements.get(i).size());
				for (int j = 0; j < used.size(); j++)
				{
					used.set(j, false);
				}
				for (int j = 0; j < group_elements.get(i).size(); j++)
				{
					System.out.print(group_elements.get(i).get(j));
					System.out.print(" ");
				}
				int fixed_count = 1;
				for (int j = 0; j < group_elements.get(i).size(); j++)
				{
					int k = j;
					int element_in_cycle = 0;
					System.out.print("(");
					while (!used.get(k))
					{
						element_in_cycle++;
						System.out.print(k);
						System.out.print(" ");
						used.set(k, true);
						k = group_elements.get(i).get(k);
					}
					System.out.print(") [");
					System.out.print(element_in_cycle);
					System.out.print("] ");
					if (element_in_cycle != 0)
					{
						fixed_count *= c;
					}
				}
				System.out.print(" fixes ");
				System.out.print(fixed_count);
				System.out.print(" configurations");
				System.out.print("\n");
				all_stablizer_count += fixed_count;
			}
			System.out.print(all_stablizer_count);
			System.out.print(" ");
			System.out.print(group_elements.size());
			System.out.print("\n");
		}
	}
}