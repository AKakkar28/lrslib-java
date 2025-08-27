package com.vertexenumeration;

import java.util.*;

public final class OptionsParser {

    public static final class Parsed {
        public final LrsDat dat;
        public final String inputPath;
        private Parsed(LrsDat d, String p){ dat=d; inputPath=p; }
    }

    public static Parsed parse(String[] args){
        LrsDat.Builder b = new LrsDat.Builder();
        String input = null;

        for (int i=0; i<args.length; i++) {
            String a = args[i];
            switch (a) {
                case "-v": b.mode(LrsDat.Mode.VE); break;        // H->V (default)
                case "-h": b.mode(LrsDat.Mode.CH); break;        // V->H
                case "-redund": b.redund(true); break;
                case "-minrep": b.minrep(true); break;
                case "-printcobasis": b.printCobasis(true); break;
                case "-integer": b.integerInput(true); break;
                case "-seed": b.seed(Long.parseLong(args[++i])); break;
                case "-maxdepth": b.maxDepth(Integer.parseInt(args[++i])); break;
                case "-threads": b.threads(Integer.parseInt(args[++i])); break;
                case "-eliminate": {
                    for (String tok : nextCSV(args[++i])) b.addEliminate(Integer.parseInt(tok));
                    break;
                }
                case "-project": {
                    for (String tok : nextCSV(args[++i])) b.addProject(Integer.parseInt(tok));
                    break;
                }
                case "-linset":
                case "-linearity": {
                    for (String tok : nextCSV(args[++i])) b.addLinearity(Integer.parseInt(tok));
                    break;
                }
                default:
                    if (a.startsWith("-")) throw new IllegalArgumentException("Unknown option: " + a);
                    if (input != null) throw new IllegalArgumentException("Multiple inputs: " + a);
                    input = a;
            }
        }
        if (input == null) throw new IllegalArgumentException("Missing input file");
        return new Parsed(b.build(), input);
    }

    private static List<String> nextCSV(String s){
        String[] parts = s.split(",");
        List<String> out = new ArrayList<>(parts.length);
        for (String p : parts) {
            String t = p.trim();
            if (!t.isEmpty()) out.add(t);
        }
        return out;
    }
}
