package ninja.mspp.proformabuilder.test;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ninja.mspp.proformabuilder.FragmentIonType;
import ninja.mspp.proformabuilder.ProForma;
import ninja.mspp.proformabuilder.ProFormaException;
import ninja.mspp.proformabuilder.ProFormaParser;
import org.junit.Test;

/**
 * Compares masses, m/z and b/y ions with values computed by pyteomics 5.0.3a5
 * ({@code pyteomics_reference.tsv}: ProForma examples of the specification and 2,000 random
 * modified peptides, generated with the same unimod.xml).
 * <p>
 * pyteomics uses an older NIST element table than this library (current NIST values), so
 * differences up to a few 1e-6 Da are expected for molecules with many S or P atoms. Inputs for
 * which pyteomics deviates from the ProForma specification are skipped and checked against
 * specification-derived values in {@link ProFormaTest}.
 * </p>
 */
public class PyteomicsComparisonTest {

	/** Monoisotopic tolerance for compositions (pyteomics uses an older NIST table). */
	static final double TOL = 3e-6;
	/**
	 * Relative tolerance of average masses: this library uses the current NIST atomic weights
	 * (rounded to 6 decimals), pyteomics an older table, which differ by up to ~2.5e-7 Da per atom
	 * (about 2e-8 per Da for peptides and glycans).
	 */
	static final double AVG_REL_TOL = 1e-7;
	/** Tolerance when pyteomics used its tabulated Unimod masses (6 decimals) instead of compositions. */
	static final double TABULATED_TOL = 6e-6;

	static final Set<String> KNOWN_DIFFERENCES = new HashSet<>(Arrays.asList(
			"BZJX[+1]",
			"<13C>PEPTIDE", "<13C><15N>ATPEILTVNSIGQLK", "<D>ATPEILTVNSIGQLK",
			"<12C>ATPEILTVNS[Formula:C-1[13C1]]IGQLK",
			"<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK",
			"<[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG",
			"<[Gln->pyro-Glu]@N-term:Q><[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG",
			"SEQUEN[Glycan:{+203.079}1Hex2]CE"));

	@Test
	public void agreesWithPyteomics() throws IOException {
		int compared = 0;
		int skipped = 0;
		double maxDiff = 0;
		List<String> failures = new ArrayList<>();
		try ( BufferedReader r = new BufferedReader(new InputStreamReader(
				getClass().getResourceAsStream("pyteomics_reference.tsv"), StandardCharsets.UTF_8)) ) {
			String line;
			while ( (line = r.readLine()) != null ) {
				if ( line.startsWith("#") || line.trim().isEmpty() )
					continue;
				String[] f = java.util.Arrays.copyOf(line.split("\t", -1), 7);
				for ( int k = 0; k < f.length; k++ )
					if ( f[k] == null )
						f[k] = "";
				String input = f[0];
				if ( !f[1].equals("OK") || KNOWN_DIFFERENCES.contains(input) || isChimeric(input) ) {
					skipped++;
					continue;
				}
				ProForma p;
				double mono;
				try {
					p = ProForma.parse(input);
					mono = p.getMonoisotopicMass().doubleValue();
				} catch ( ProFormaException e ) {
					skipped++;
					continue;
				}
				compared++;
				boolean glycan = input.contains("Glycan:");
				boolean compositionKnown = !f[3].isEmpty();
				double tol = compositionKnown ? TOL : TABULATED_TOL;
				if ( !glycan || compositionKnown ) {
					maxDiff = Math.max(maxDiff, check(failures, input, "mono", f[2], mono, tol));
					if ( !f[4].isEmpty() )
						check(failures, input, "m/z", f[4], p.getMz().doubleValue(), tol);
				}
				if ( compositionKnown && !p.isAverageMassApproximate() && !input.contains(":z") )
					check(failures, input, "avg", f[3], p.getAverageMass().doubleValue(),
							AVG_REL_TOL * Double.parseDouble(f[3]));
				if ( !f[5].isEmpty() && p.getModel().getChains().size() == 1 && !glycan ) {
					series(failures, input, "b", f[5], p.getFragments(FragmentIonType.B, 1, false, true));
					series(failures, input, "y", f[6], p.getFragments(FragmentIonType.Y, 1, false, true));
				}
			}
		}
		System.out.printf("pyteomics comparison: %d compared, %d skipped, %d differences, max mono diff %.2e%n",
				compared, skipped, failures.size(), maxDiff);
		for ( int k = 0; k < Math.min(20, failures.size()); k++ )
			System.out.println(failures.get(k));
		assertTrue("too few comparisons: " + compared, compared > 1500);
		assertTrue(failures.size() + " differences, first: " + (failures.isEmpty() ? "" : failures.get(0)), failures.isEmpty());
	}

	private static boolean isChimeric(String s) {
		if ( s.indexOf('+') < 0 )
			return false;
		try {
			return ProFormaParser.parseChimeric(s).size() > 1;
		} catch ( ProFormaException e ) {
			return false;
		}
	}

	private static void series(List<String> failures, String input, String name, String ref, List<BigDecimal> actual) {
		String[] parts = ref.split(",");
		if ( parts.length != actual.size() ) {
			failures.add(input + ": " + name + " ion count " + parts.length + " vs " + actual.size());
			return;
		}
		for ( int k = 0; k < parts.length; k++ )
			check(failures, input, name + (k + 1), parts[k], actual.get(k).doubleValue(), TABULATED_TOL);
	}

	private static double check(List<String> failures, String input, String what, String expected, double actual, double tol) {
		double diff = Math.abs(Double.parseDouble(expected) - actual);
		if ( diff > tol )
			failures.add(String.format("%s: %s pyteomics=%s java=%.7f diff=%.2e", input, what, expected, actual, diff));
		return diff;
	}
}
