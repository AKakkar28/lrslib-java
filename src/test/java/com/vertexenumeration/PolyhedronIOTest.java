package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

/**
 * Tests reading and writing of polyhedra in lrslib format.  These tests
 * ensure that the {@link Polyhedron#readFromFile(String)} and
 * {@link Polyhedron#write(PrintWriter)} methods correctly parse and
 * serialise the standard input format used by lrslib.
 */
public class PolyhedronIOTest {

    @Test
    public void testReadWriteSquareHRep() throws IOException {
        // A simple 2D square in H-representation (4 constraints, 2 variables + homogenising column)
        String input = String.join(System.lineSeparator(),
                "H-representation",
                "begin",
                "4 3 integer",
                " 1  1  0",
                " 1  0  1",
                " 1 -1  0",
                " 1  0 -1",
                "end",
                "");
        // Write the input to a temporary file
        Path temp = Files.createTempFile("poly", ".ine");
        Files.writeString(temp, input);
        // Read the polyhedron back in
        Polyhedron h = Polyhedron.readFromFile(temp.toString());
        assertEquals(Polyhedron.Type.H, h.getType());
        assertEquals(4, h.getRowCount());
        assertEquals(3, h.getColCount());
        assertTrue(h.isIntegerData());
        // Write out again and verify header content
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        h.write(pw);
        pw.flush();
        String out = sw.toString();
        assertTrue(out.contains("H-representation"));
        assertTrue(out.contains("4 3 integer"));
        // Clean up temp file
        Files.deleteIfExists(temp);
    }
}