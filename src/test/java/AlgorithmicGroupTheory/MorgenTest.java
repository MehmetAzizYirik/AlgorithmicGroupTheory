package AlgorithmicGroupTheory;

import org.junit.Test;

public class MorgenTest {

    @Test
    public void main() throws Exception {
        String[] argument = {"-f", "C4OH8", "-v", "-d", "result1"};
        System.out.println("First call");
        MORGEN.main(argument);
        System.out.println("Second call");
        MORGEN.main(argument);
    }
}
