package AlgorithmicGroupTheory;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class FormulaGenerator {
    private String carbon;
    private String hydrogen;
    private String nitrogen;
    private String oxygen;
    private String phosphorus;

    private void parseArgs(String[] args) throws ParseException {
        Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            this.carbon = cmd.getOptionValue("carbon");
            this.hydrogen = cmd.getOptionValue("hydrogen");
            this.nitrogen = cmd.getOptionValue("nitrogen");
            this.oxygen = cmd.getOptionValue("oxygen");
            this.phosphorus = cmd.getOptionValue("phosphorus");
            if (Objects.isNull(carbon) && Objects.isNull(hydrogen)) {
                displayHelp(options);
            }
        } catch (ParseException e) {
            displayHelp(options);
            throw new ParseException("Problem parsing command line");
        }
    }

    private void displayHelp(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(null);
        String header =
                "\nGenerates formulas for a given components."
                        + " The input is component strings."
                        + "For example '-c 6 -h 10 -n 3 -o 2 -p 0'.\n\n";
        String footer =
                "\nPlease report issues at https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory";
        formatter.printHelp(
                "java -cp MORGEN.jar AlgorithmicGroupTheory.FormulaGenerator",
                header,
                options,
                footer,
                true);
    }

    private Options setupOptions(String[] args) {
        Options options = new Options();
        Option carbon =
                Option.builder("c")
                        .required(false)
                        .hasArg()
                        .longOpt("carbon")
                        .desc("Carbon")
                        .build();
        options.addOption(carbon);
        Option hydrogen =
                Option.builder("h")
                        .required(false)
                        .hasArg()
                        .longOpt("hydrogen")
                        .desc("Hydrogen")
                        .build();
        options.addOption(hydrogen);
        Option nitrogen =
                Option.builder("n")
                        .required(false)
                        .hasArg()
                        .longOpt("nitrogen")
                        .desc("Nitrogen")
                        .build();
        options.addOption(nitrogen);
        Option oxygen =
                Option.builder("o")
                        .required(false)
                        .hasArg()
                        .longOpt("oxygen")
                        .desc("Oxygen")
                        .build();
        options.addOption(oxygen);
        Option phosphorus =
                Option.builder("p")
                        .required(false)
                        .hasArg()
                        .longOpt("phosphorus")
                        .desc("Phosphorus")
                        .build();
        options.addOption(phosphorus);
        return options;
    }

    public static void main(String[] args) throws IOException, ParseException {
        FormulaGenerator gen = new FormulaGenerator();
        gen.parseArgs(args);
        gen.generateFormulas();
    }

    private void generateFormulas() throws IOException {
        Integer carbonValue = getIntegerValue(carbon);
        Integer hydrogenValue = getIntegerValue(hydrogen);
        Integer oxygenValue = getIntegerValue(oxygen);
        Integer nitrogenValue = getIntegerValue(nitrogen);
        Integer phosphorusValue = getIntegerValue(phosphorus);

        if (carbonValue == 0 && hydrogenValue == 0) {
            System.err.println("Carbon or Hydrogen should be defined.");
        } else if (carbonValue > 6) {
            System.err.println("Carbon maximum is 6.");
        } else if (hydrogenValue > 10) {
            System.err.println("Hydrogen maximum is 10.");
        } else {
            List<String> generatedFormulas = new ArrayList<>();
            for (int index = 1; index <= carbonValue; index += 1) {
                if (hydrogenValue > 0) {
                    for (int index2 = 1; index2 <= hydrogenValue; index2 += 1) {
                        generatedFormulas.add(
                                getItem("C", index)
                                        + getItem("H", index2)
                                        + getItem("O", oxygenValue)
                                        + getItem("N", nitrogenValue)
                                        + getItem("P", phosphorusValue));
                    }
                } else {
                    generatedFormulas.add(
                            getItem("C", index)
                                    + getItem("O", oxygenValue)
                                    + getItem("N", nitrogenValue)
                                    + getItem("P", phosphorusValue));
                }
            }
            Files.write(
                    Paths.get("generated-formulas.txt"),
                    String.join("\n", generatedFormulas).getBytes(StandardCharsets.UTF_8));
            System.out.println(generatedFormulas.size() + " formulas have been generated.");
        }
    }

    private String getItem(String letter, Integer index) {
        return Objects.isNull(index) || index == 0 ? "" : (index == 1 ? letter : letter + index);
    }

    private Integer getIntegerValue(String intValue) {
        return Objects.nonNull(intValue) && intValue.matches("\\d+")
                ? Integer.parseInt(intValue)
                : 0;
    }
}
