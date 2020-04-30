package AlgorithmicGroupTheory;

import org.openscience.cdk.group.AtomContainerDiscretePartitionRefiner;
import org.openscience.cdk.group.PartitionRefinement;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class EnumerationMethods {
	public static IChemObjectBuilder builder =  SilentChemObjectBuilder.getInstance();
	
	/**
	 * General Functions
	 */
	
	/**
	 * Calculating automorphism group for a given molecular formula.
	 * @param formula String Molecular formula
	 * @return PermutationGroup automorphism group of the molecule
	 */
	 
	public static PermutationGroup automorphismGroup(String formula) {
		IAtomContainer ac= MolecularFormulaManipulator.getAtomContainer(formula, builder);
		AtomContainerDiscretePartitionRefiner refiner = PartitionRefinement.forAtoms().create();
	    PermutationGroup autoGroup = refiner.getAutomorphismGroup(ac);
	    return autoGroup;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
