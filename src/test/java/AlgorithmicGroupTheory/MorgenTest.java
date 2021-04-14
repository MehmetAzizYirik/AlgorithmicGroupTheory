package AlgorithmicGroupTheory;

import static org.junit.Assert.*;

import java.io.IOException;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;

public class MorgenTest {

    /**
     * This is the junit test class for MORGEN. Randomly selected 40 molecular formulas are used.
     * The number of generated structures are checked. The number of isomers are also tested with
     * MOLGEN algorithm. MORGEN generates same number of isomers like MOLGEN.
     */
    private MORGEN morgen;

    @Before
    public void before() {
        morgen = new MORGEN();
    }

    @Test
    public void test_C3Cl2H4() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3Cl2H4";
        morgen.run();
        assertEquals(7l, MORGEN.count);
    }

    @Test
    public void test_C2NO2H5() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C2NO2H5";
        morgen.run();
        assertEquals(84, MORGEN.count);
    }

    @Test
    public void test_C6H6() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6H6";
        morgen.run();
        assertEquals(217, MORGEN.count);
    }

    @Test
    public void test_C3O3H4() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3O3H4";
        morgen.run();
        assertEquals(152, MORGEN.count);
    }

    @Test
    public void test_Cl2C5H4() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "Cl2C5H4";
        morgen.run();
        assertEquals(217, MORGEN.count);
    }

    @Test
    public void test_C5H9ClO() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H9ClO";
        morgen.run();
        assertEquals(334, MORGEN.count);
    }

    @Test
    public void test_C7H6N2O() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C7H6N2O";
        morgen.run();
        assertEquals(10_504_307, (long) MORGEN.count.get());
    }

    @Test
    public void test_C6OF2H12() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6OF2H12";
        morgen.run();
        assertEquals(536, MORGEN.count);
    }

    @Test
    public void test_C7H10() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C7H10";
        morgen.run();
        assertEquals(575, MORGEN.count);
    }

    @Test
    public void test_C6O2H12() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6O2H12";
        morgen.run();
        assertEquals(1313, MORGEN.count);
    }

    @Test
    public void test_C5H2N2O3() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H2N2O3";
        morgen.run();
        assertEquals(1_197_634, (long) MORGEN.count.get());
    }

    @Test
    public void test_F2P3BrNO2H() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "F2P3BrNO2H";
        morgen.run();
        assertEquals(1958, MORGEN.count);
    }

    @Test
    public void test_C6OH6() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6OH6";
        morgen.run();
        assertEquals(2237, MORGEN.count);
    }

    @Test
    public void test_C5H6BrN() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H6BrN";
        morgen.run();
        assertEquals(2325, MORGEN.count);
    }

    @Test
    public void test_C6H7F2I() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6H7F2I";
        morgen.run();
        assertEquals(3523, MORGEN.count);
    }

    @Test
    public void test_C5F2O2H2() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5F2O2H2";
        morgen.run();
        assertEquals(7094, MORGEN.count);
    }

    @Test
    public void test_C7OH10() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C7OH10";
        morgen.run();
        assertEquals(7166, MORGEN.count);
    }

    @Test
    public void test_C4ClHF2O3() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4ClHF2O3";
        morgen.run();
        assertEquals(7346, MORGEN.count);
    }

    @Test
    public void test_C4O5H6() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4O5H6";
        morgen.run();
        assertEquals(8070, MORGEN.count);
    }

    @Test
    public void test_C5ClHF2O2() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5ClHF2O2";
        morgen.run();
        assertEquals(12400, MORGEN.count);
    }

    @Test
    public void test_C5H10BrF2OP() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H10BrF2OP";
        morgen.run();
        assertEquals(15009, MORGEN.count);
    }

    @Test
    public void test_C9H12() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C9H12";
        morgen.run();
        assertEquals(19983, MORGEN.count);
    }

    @Test
    public void test_C6H10O2Br2() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6H10O2Br2";
        morgen.run();
        assertEquals(24201, MORGEN.count);
    }

    @Test
    public void test_C10H16() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C10H16";
        morgen.run();
        assertEquals(24938, MORGEN.count);
    }

    @Test
    public void test_C6H6ClOI() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C6H6ClOI";
        morgen.run();
        assertEquals(30728, MORGEN.count);
    }

    @Test
    public void test_C4H5O2Br2N() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4H5O2Br2N";
        morgen.run();
        assertEquals(41067, MORGEN.count);
    }

    @Test
    public void test_C4H10NOSP() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4H10NOSP";
        morgen.run();
        assertEquals(52151, MORGEN.count);
    }

    @Test
    public void test_C7O2H10() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C7O2H10";
        morgen.run();
        assertEquals(54641, MORGEN.count);
    }

    @Test
    public void test_P3O3NCl2() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "P3O3NCl2";
        morgen.run();
        assertEquals(665, MORGEN.count);
    }

    @Test
    public void test_C5H5SI5() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H5SI5";
        morgen.run();
        assertEquals(2619, MORGEN.count);
    }

    @Test
    public void test_C3O3NH5() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3O3NH5";
        morgen.run();
        assertEquals(2644, MORGEN.count);
    }

    @Test
    public void test_C5H9ClOS() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H9ClOS";
        morgen.run();
        assertEquals(3763, MORGEN.count);
    }

    @Test
    public void test_C3NO2SH7() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3NO2SH7";
        morgen.run();
        assertEquals(3838, MORGEN.count);
    }

    @Test
    public void test_C4H8Cl3O2P() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4H8Cl3O2P";
        morgen.run();
        assertEquals(9313, MORGEN.count);
    }

    @Test
    public void test_C5H2F2SO() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H2F2SO";
        morgen.run();
        assertEquals(13446, MORGEN.count);
    }

    @Test
    public void test_C7H11ClS() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C7H11ClS";
        morgen.run();
        assertEquals(15093, MORGEN.count);
    }

    @Test
    public void test_C4NO3H7() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4NO3H7";
        morgen.run();
        assertEquals(18469, MORGEN.count);
    }

    @Test
    public void test_C4H5O2F2P() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C4H5O2F2P";
        morgen.run();
        assertEquals(41067, MORGEN.count);
    }

    @Test
    public void test_C3N3O2H7() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3N3O2H7";
        morgen.run();
        assertEquals(45626, MORGEN.count);
    }

    @Test
    public void test_C5N3H9() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5N3H9";
        morgen.run();
        assertEquals(46125, MORGEN.count);
    }

    @Test
    public void test_C3O6PH5() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C3O6PH5";
        morgen.run();
        assertEquals(51323, MORGEN.count);
    }

    @Test
    public void test_C5H5POBr2() throws IOException, CDKException, CloneNotSupportedException {
        morgen.formula = "C5H5POBr2";
        morgen.run();
        assertEquals(62886, MORGEN.count);
    }
}
