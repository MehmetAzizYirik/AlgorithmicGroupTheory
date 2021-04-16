package AlgorithmicGroupTheory;

import static org.junit.Assert.*;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.stream.IntStream;

import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;

public class MorgenAverageTest {
	
	/**
	 * This is the junit test class for MORGEN. Randomly selected 40 molecular
	 * formulas are used. The number of generated structures are checked. The
	 * number of isomers are also tested with MOLGEN algorithm. MORGEN generates 
	 * same number of isomers like MOLGEN.
	 */
	
        private MORGEN morgen;

        @Before
        public void before() {
            morgen = new MORGEN();
            MORGEN.count.set(0);
        }

	@Test

	public void test_C3Cl2H4() {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 7, "C3Cl2H4")).average().getAsDouble();
		System.out.println("C3Cl2H4, average time - " + average);
	}

	private long execMorgen(long i, int expectedCount, String formula) {
		Instant start = Instant.now();
		morgen = new MORGEN();
		MORGEN.count.set(0);
		morgen.formula = formula;
		try {
			morgen.run();
			assertEquals(expectedCount, MORGEN.count.get());
		} catch (Exception e) {
			e.printStackTrace();
		}
		Instant finish = Instant.now();
		return Duration.between(start, finish).toMillis();
	}

	@Test

	public void test_C2NO2H5() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 84, "C2NO2H5")).average().getAsDouble();
		System.out.println("C2NO2H5, average time - " + average);
	}

	@Test

	public void test_C6H6() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 217, "C6H6")).average().getAsDouble();
		System.out.println("C6H6, average time - " + average);
	}

        @Test

	public void test_C3O3H4() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 152, "C3O3H4")).average().getAsDouble();
		System.out.println("C3O3H4, average time - " + average);
	}

	@Test

	public void test_Cl2C5H4() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 217, "Cl2C5H4")).average().getAsDouble();
		System.out.println("Cl2C5H4, average time - " + average);
	}

	@Test

	public void test_C5H9ClO() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 334, "C5H9ClO")).average().getAsDouble();
		System.out.println("C5H9ClO, average time - " + average);
	}

	@Test

	public void test_C7H6N2O() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 10_504_307, "C7H6N2O")).average().getAsDouble();
		System.out.println("C7H6N2O, average time - " + average);
	}

	@Test

	public void test_C6OF2H12() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 536, "C6OF2H12")).average().getAsDouble();
		System.out.println("C6OF2H12, average time - " + average);
	}

	@Test

	public void test_C7H10() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 575, "C7H10")).average().getAsDouble();
		System.out.println("C7H10, average time - " + average);
	}

	@Test

	public void test_C6O2H12() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 1313, "C6O2H12")).average().getAsDouble();
		System.out.println("C6O2H12, average time - " + average);
	}

	@Test

	public void test_C5H2N2O3() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 1_197_634, "C5H2N2O3")).average().getAsDouble();
		System.out.println("C5H2N2O3, average time - " + average);
	}

	@Test

	public void test_F2P3BrNO2H() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 1958, "F2P3BrNO2H")).average().getAsDouble();
		System.out.println("F2P3BrNO2H, average time - " + average);
	}

	@Test

	public void test_C6OH6() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 2237, "C6OH6")).average().getAsDouble();
		System.out.println("C6OH6, average time - " + average);
	}

	@Test

	public void test_C5H6BrN() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 2325, "C5H6BrN")).average().getAsDouble();
		System.out.println("C5H6BrN, average time - " + average);
	}

	@Test

	public void test_C6H7F2I() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 3523, "C6H7F2I")).average().getAsDouble();
		System.out.println("C6H7F2I, average time - " + average);
	}

	@Test

	public void test_C5F2O2H2() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 7094, "C5F2O2H2")).average().getAsDouble();
		System.out.println("C5F2O2H2, average time - " + average);
	}

	@Test

	public void test_C7OH10() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 7166, "C7OH10")).average().getAsDouble();
		System.out.println("C7OH10, average time - " + average);
	}

	@Test

	public void test_C4ClHF2O3() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 7346, "C4ClHF2O3")).average().getAsDouble();
		System.out.println("C4ClHF2O3, average time - " + average);
	}

	@Test

	public void test_C4O5H6() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 8070, "C4O5H6")).average().getAsDouble();
		System.out.println("C4O5H6, average time - " + average);
	}

	@Test

	public void test_C5ClHF2O2() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 12400, "C5ClHF2O2")).average().getAsDouble();
		System.out.println("C5ClHF2O2, average time - " + average);
	}

	@Test

	public void test_C5H10BrF2OP() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 15009, "C5H10BrF2OP")).average().getAsDouble();
		System.out.println("C5H10BrF2OP, average time - " + average);
	}

	@Test

	public void test_C9H12() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 19983, "C9H12")).average().getAsDouble();
		System.out.println("C9H12, average time - " + average);
	}

	@Test

	public void test_C6H10O2Br2() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 24201, "C6H10O2Br2")).average().getAsDouble();
		System.out.println("C6H10O2Br2, average time - " + average);
	}
	
	@Test
	
	public void test_C10H16() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 24938, "C10H16")).average().getAsDouble();
		System.out.println("C10H16, average time - " + average);
	}

	@Test

	public void test_C6H6ClOI() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 30728, "C6H6ClOI")).average().getAsDouble();
		System.out.println("C6H6ClOI, average time - " + average);
	}

	@Test

	public void test_C4H5O2Br2N() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 217, "C6H6")).average().getAsDouble();
		System.out.println("C6H6, average time - " + average);
		morgen.formula="C4H5O2Br2N";
		morgen.run();
		assertEquals(41067,(long) MORGEN.count.get());
	}

	@Test

	public void test_C4H10NOSP() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 217, "C6H6")).average().getAsDouble();
		System.out.println("C6H6, average time - " + average);
		morgen.formula="C4H10NOSP";
		morgen.run();
		assertEquals(52151,(long) MORGEN.count.get());
	}

	@Test

	public void test_C7O2H10() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 54641, "C7O2H10")).average().getAsDouble();
		System.out.println("C7O2H10, average time - " + average);
	}

	@Test

	public void test_P3O3NCl2() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 665, "P3O3NCl2")).average().getAsDouble();
		System.out.println("P3O3NCl2, average time - " + average);
	}

	@Test

	public void test_C5H5SI5() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 2619, "C5H5SI5")).average().getAsDouble();
		System.out.println("C5H5SI5, average time - " + average);
	}

	@Test

	public void test_C3O3NH5() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 2644, "C3O3NH5")).average().getAsDouble();
		System.out.println("C3O3NH5, average time - " + average);
	}

	@Test

	public void test_C5H9ClOS() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 3763, "C5H9ClOS")).average().getAsDouble();
		System.out.println("C5H9ClOS, average time - " + average);
	}

	@Test

	public void test_C3NO2SH7() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 3838, "C3NO2SH7")).average().getAsDouble();
		System.out.println("C3NO2SH7, average time - " + average);
	}

	@Test

	public void test_C4H8Cl3O2P() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 9313, "C4H8Cl3O2P")).average().getAsDouble();
		System.out.println("C4H8Cl3O2P, average time - " + average);
	}

	@Test

	public void test_C5H2F2SO() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 13446, "C5H2F2SO")).average().getAsDouble();
		System.out.println("C5H2F2SO, average time - " + average);
	}

	@Test

	public void test_C7H11ClS() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 15093, "C7H11ClS")).average().getAsDouble();
		System.out.println("C7H11ClS, average time - " + average);
	}

	@Test

	public void test_C4NO3H7() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 18469, "C4NO3H7")).average().getAsDouble();
		System.out.println("C4NO3H7, average time - " + average);
	}

	@Test

	public void test_C4H5O2F2P() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 41067, "C4H5O2F2P")).average().getAsDouble();
		System.out.println("C4H5O2F2P, average time - " + average);
	}

	@Test

	public void test_C3N3O2H7() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 45626, "C3N3O2H7")).average().getAsDouble();
		System.out.println("C3N3O2H7, average time - " + average);
	}

	@Test

	public void test_C5N3H9() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 46125, "C5N3H9")).average().getAsDouble();
		System.out.println("C5N3H9, average time - " + average);
	}

	@Test

	public void test_C3O6PH5() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 51323, "C3O6PH5")).average().getAsDouble();
		System.out.println("C3O6PH5, average time - " + average);
	}

	@Test

	public void test_C5H5POBr2() throws IOException, CDKException, CloneNotSupportedException {
		double average = IntStream.rangeClosed(1, 5).asLongStream().map(
				i -> execMorgen(i, 62886, "C5H5POBr2")).average().getAsDouble();
		System.out.println("C5H5POBr2, average time - " + average);
	}
}
