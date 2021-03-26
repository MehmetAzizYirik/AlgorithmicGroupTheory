package AlgorithmicGroupTheory;

import java.util.Arrays;
import java.util.List;
import org.junit.Test;

public class MorgenTest {

    @Test
    public void main() {
        String[] argument = {"-f", "C5N3H9.", "-v", "-d", "result1"};
//        String[] argument2 = {"-f", "C2H5NO2", "-v", "-d", "result2"};
//        String[] argument3 = {"-f", "C4H7NO3", "-v", "-d", "result3"};
//        String[] argument4 = {"-f", "C5H9N3", "-v", "-d", "result4"};
//        List<String[]> arguments = Arrays.asList(argument1, /*argument2,*/ argument3, argument4);
        MORGEN.main(argument);
    }
}
