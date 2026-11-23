package ninja.mspp.proformabuilder.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import org.glycoinfo.ms.GlycanMassUtility.om.IMassElement;
import org.glycoinfo.ms.GlycanMassUtility.om.IonCloud;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;
import org.glycoinfo.ms.GlycanMassUtility.om.Peptide;
import ninja.mspp.proformabuilder.Compositions;
import ninja.mspp.proformabuilder.FragmentIonType;
import ninja.mspp.proformabuilder.GlycanComposition;
import ninja.mspp.proformabuilder.Model.PeptidoformIon;
import ninja.mspp.proformabuilder.ProForma;
import ninja.mspp.proformabuilder.ProFormaException;
import ninja.mspp.proformabuilder.ProFormaFormula;
import ninja.mspp.proformabuilder.ProFormaParser;
import ninja.mspp.proformabuilder.ProFormaWriter;
import org.glycoinfo.ms.GlycanMassUtility.utils.CompoundParserType;
import org.glycoinfo.ms.GlycanMassUtility.utils.MassUtils;
import org.junit.Test;

public class ProFormaTest {

	static final double D = 1e-6;
	static final double H2O = Compositions.H2O.getMonoisotopicMass().doubleValue();

	static double mono(String s) {
		return ProForma.parse(s).getMonoisotopicMass().doubleValue();
	}

	static double formula(String f) {
		return ProFormaFormula.parse(f).getMonoisotopicMass().doubleValue();
	}

	static void assertThrows(String input, Runnable r) {
		try {
			r.run();
			fail("IllegalArgumentException expected for " + input);
		} catch ( IllegalArgumentException e ) {
			// expected (ProFormaException is a subclass)
		}
	}

	static List<String> resourceLines(String name) throws IOException {
		List<String> out = new ArrayList<>();
		try ( BufferedReader r = new BufferedReader(new InputStreamReader(
				ProFormaTest.class.getResourceAsStream(name), StandardCharsets.UTF_8)) ) {
			String line;
			while ( (line = r.readLine()) != null )
				if ( !line.trim().isEmpty() && !line.startsWith("#") )
					out.add(line);
		}
		return out;
	}

	@Test
	public void referencePeptide() {
		ProForma p = ProForma.parse("PEPS[Phospho]TC[Carbamidomethyl]M[Oxidation]K/2");
		assertEquals("PEPSTCMK", p.getSequence());
		assertEquals("C38H65N10O18P1S2", Compositions.toProFormaFormula(p.getComposition()));
		assertEquals(1044.365734, p.getMonoisotopicMass().doubleValue(), D);
		assertEquals(1045.0838, p.getAverageMass().doubleValue(), 1e-4);
		assertEquals(523.190143, p.getMz().doubleValue(), D);
		assertFalse(p.isAverageMassApproximate());
		double[] b = {98.0600, 227.1026, 324.1554, 491.1538, 592.2014, 752.2321, 899.2675};
		double[] y = {147.1128, 294.1482, 454.1789, 555.2265, 722.2249, 819.2777, 948.3202};
		List<BigDecimal> fb = p.getFragments(FragmentIonType.B, 1);
		List<BigDecimal> fy = p.getFragments(FragmentIonType.Y, 1);
		for ( int k = 0; k < b.length; k++ ) {
			assertEquals(b[k], fb.get(k).doubleValue(), 1e-4);
			assertEquals(y[k], fy.get(k).doubleValue(), 1e-4);
		}
		assertEquals("PEPS[Phospho]TC[Carbamidomethyl]M[Oxidation]K/2", p.toString());
	}

	@Test
	public void equivalentNotations() {
		double e = 1044.365734;
		assertEquals(e, mono("PEPS[UNIMOD:21]TC[UNIMOD:4]M[UNIMOD:35]K"), D);
		assertEquals(e, mono("PEPS[U:Phospho]TC[U:Carbamidomethyl]M[U:Oxidation]K"), D);
		assertEquals(e, mono("PEPS[Formula:HPO3]TC[Formula:C2H3NO]M[Formula:O]K"), D);
		assertEquals(e, mono("<[Carbamidomethyl]@C>PEPS[Phospho]TCM[Oxidation]K"), D);
		assertEquals(e, mono("[Phospho]?PEPSTC[Carbamidomethyl]M[Oxidation]K"), D);
		assertEquals(e, mono("{Phospho}PEPSTC[Carbamidomethyl]M[Oxidation]K"), D);
		assertEquals(e, mono("pepS[phospho]tc[carbamidomethyl]m[oxidation]k"), D);
		// Unimod full names (e.g. Phosphorylation for Phospho) are not ProForma names
		assertThrows("Phosphorylation", new Runnable() { public void run() { mono("PEPS[Phosphorylation]TCMK"); } });
		assertThrows("phosphorylation", new Runnable() { public void run() { mono("PEPS[phosphorylation]TCMK"); } });
		assertEquals(e, mono("PEPS[+79.966331]TC[+57.021464]M[+15.994915]K"), 1e-5);
		assertTrue(ProForma.parse("PEPS[+79.966331]TCMK").isAverageMassApproximate());
		assertEquals(840.386513, mono("[Acetyl]-PEPTIDE-[Amidated]"), D);
	}

	@Test
	public void formulas() {
		assertEquals(633.368958, mono("PEPK[Formula:[13C6][15N2]C-6N-2]R"), D);
		assertEquals(8.014199, formula("[13C6][15N2]C-6N-2"), D);
		Molecule c13 = ProFormaFormula.parse("[13C1]");
		assertEquals(c13.getMonoisotopicMass(), c13.getAverageMass());
		assertEquals(formula("C2H4O2"), formula("C2 H4 O2"), 0);
		assertEquals(formula("[2H1]"), formula("D"), 0);
		assertEquals(63.929142, formula("Zn"), 1e-6); // element added from the NIST table
		assertThrows("Xy2", new Runnable() { public void run() { ProFormaFormula.parse("Xy2"); } });
	}

	@Test
	public void globalIsotopesAndFixedModifications() {
		double d13c = formula("[13C1]") - formula("C");
		double d15n = formula("[15N1]") - formula("N");
		assertEquals(799.359964 + 34 * d13c, mono("<13C>PEPTIDE"), 1e-5);
		assertEquals(799.359964 + 34 * d13c + 7 * d15n, mono("<13C><15N>PEPTIDE"), 1e-5);
		assertEquals(mono("ATPEILTVNS[Formula:C-1[13C1]]IGQLK"), mono("<12C>ATPEILTVNS[Formula:C-1[13C1]]IGQLK"), D);
		double tmt = mono("[TMT6plex]-K") - mono("K");
		assertEquals(mono("ATPEILTCNSIGCLK") + 2 * tmt, mono("<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK"), D);
		double ox = formula("O");
		assertEquals(mono("QATPEILTWCNSIGCLKG") + 2 * ox, mono("<[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG"), D);
		assertEquals(mono("MTPEILTCNSIGCLK") + 3 * ox, mono("<[Oxidation]@C,M>MTPEILTCNSIGCLK"), D);
	}

	static double avg(String s) {
		return ProForma.parse(s).getAverageMass().doubleValue();
	}

	static double avgFormula(String f) {
		return ProFormaFormula.parse(f).getAverageMass().doubleValue();
	}

	/** Average masses, checked against values computed independently with the existing classes of this library. */
	@Test
	public void averageMasses() {
		// unmodified peptide: same as Peptide (799.823893 in the README)
		assertEquals(Peptide.parse("PEPTIDE").getAverageMass().doubleValue(), avg("PEPTIDE"), 1e-9);
		assertEquals(799.823893, avg("PEPTIDE"), 1e-6);
		// modification by name: peptide + average mass of the delta composition
		assertEquals(Peptide.parse("PEPTMIDE").getAverageMass().doubleValue() + avgFormula("O"), avg("PEPTM[Oxidation]IDE"), 1e-9);
		assertEquals(avg("PEPTMIDE") + avgFormula("O"), avg("PEPTM[MOD:00425]IDE"), 1e-9);
		// reference peptide (pyteomics: 1045.0838)
		assertEquals(1045.0838, avg("PEPS[Phospho]TC[Carbamidomethyl]M[Oxidation]K"), 1e-4);
		// glycopeptide: same as the glycopeptide parser of this library
		IMassElement gpep = CompoundParserType.GPEP.parse("NEEYNK|Hex(5)HexNAc(4)Neu5Ac(1)");
		assertEquals(gpep.getAverageMass().doubleValue(), avg("NEEYN[Glycan:Hex5HexNAc4NeuAc1]K"), 1e-9);
		// specified isotopes contribute their exact mass also to the average mass
		double label = 6 * formula("[13C1]") + 2 * formula("[15N1]") - 6 * avgFormula("C") - 2 * avgFormula("N");
		assertEquals(avg("PEPKR") + label, avg("PEPK[Formula:[13C6][15N2]C-6N-2]R"), 1e-9);
		// global isotope <13C>: every carbon is 13C (PEPTIDE has 34 carbons)
		assertEquals(avg("PEPTIDE") + 34 * (formula("[13C1]") - avgFormula("C")), avg("<13C>PEPTIDE"), 1e-9);
		// charged formula: electrons removed as for the monoisotopic mass
		double electron = 0.0005486;
		assertEquals(avg("PEPTIDE") + avgFormula("Zn") - 2 * electron, avg("PEPT[Formula:Zn:z+2]IDE"), 1e-9);
		// mass deltas: the average mass is approximate (the delta is added as given) and flagged
		ProForma delta = ProForma.parse("PEPTM[+15.994915]IDE");
		assertTrue(delta.isAverageMassApproximate());
		assertEquals(avg("PEPTMIDE") + 15.994915, delta.getAverageMass().doubleValue(), 1e-9);
		assertFalse(ProForma.parse("PEPTM[Oxidation]IDE").isAverageMassApproximate());
	}

	@Test
	public void ambiguousResidues() {
		assertThrows("B", new Runnable() { public void run() { mono("B"); } });
		List<BigDecimal> b = ProForma.parse("B").getMonoisotopicMassCandidates();
		assertEquals(2, b.size());
		assertEquals(mono("N"), b.get(0).doubleValue(), D);
		assertEquals(mono("D"), b.get(1).doubleValue(), D);
		assertEquals(3, ProForma.parse("BB").getMonoisotopicMassCandidates().size());
		assertEquals(3, ProForma.parse("BZJX[+1]").getMonoisotopicMassCandidates().size()); // D+Q = N+E
		assertEquals(mono("I"), mono("J"), 0);
		assertEquals(1 + H2O, mono("X[+1]"), D);
	}

	@Test
	public void glycans() {
		double pep = mono("NEEYNK");
		assertEquals(795.339897, pep, D);
		assertEquals(1927.570468, mono("NEEYN[Glycan:HexS2Hex4]K"), D);
		assertEquals(892.317215, mono("NEEYN[Glycan:HexNAc 2 Hex 3]K") - pep, D);
		assertEquals(mono("N[Glycan:Hex2HexNAc1]"), mono("N[Glycan:Hex2HexNAc]"), D);
		assertEquals(242.009638, mono("X[Glycan:HexS]") - H2O, D);
		assertEquals(158.021523, mono("X[Glycan:en,aHex]") - H2O, D);
		assertEquals(79.956815, mono("X[Glycan:Sulfo]") - H2O, D);
		// Sug, Tri and Tet are parsed and written, but have no mass
		for ( final String name : GlycanComposition.UNSUPPORTED_NAMES ) {
			assertEquals("X[Glycan:" + name + "]", ProFormaWriter.write(ProFormaParser.parse("X[Glycan:" + name + "]")));
			assertThrows(name, new Runnable() { public void run() { mono("X[Glycan:" + name + "]"); } });
		}
		assertEquals(mono("SEQUEN[Glycan:HexNAc1Hex2]CE"), mono("SEQUEN[Glycan:{C8H13N1O5}1Hex2]CE"), D);
		assertEquals(333.105981, mono("X[Glycan:{C13H19N1O9}1]") - H2O, D);
		assertThrows("Sialic Acid", new Runnable() { public void run() { ProForma.parse("N[Glycan:Sialic Acid1]"); } });
		assertThrows("HexA", new Runnable() { public void run() { ProForma.parse("N[Glycan:HexA1]"); } });
		// same value as the glycopeptide parser of this library
		IMassElement gpep = CompoundParserType.GPEP.parse("NEEYNK|Hex(5)HexNAc(4)Neu5Ac(1)");
		assertEquals(gpep.getMonoisotopicMass().doubleValue(), mono("NEEYN[Glycan:Hex5HexNAc4NeuAc1]K"), D);
	}

	/**
	 * Every supported ProForma monosaccharide name resolves in the dictionaries to the composition listed for glypy.
	 */
	@Test
	public void glycanNamesAgreeWithGlypy() throws IOException {
		int n = 0;
		for ( String line : resourceLines("glycan_names_glypy.tsv") ) {
			String[] f = java.util.Arrays.copyOf(line.split("\t", -1), 6);
			if ( GlycanComposition.UNSUPPORTED_NAMES.contains(f[0]) )
				continue;
			assertEquals(f[0], ProFormaFormula.parse(f[3]), GlycanComposition.residue(f[0]));
			n++;
		}
		assertEquals(GlycanComposition.NAMES.size() - GlycanComposition.UNSUPPORTED_NAMES.size(), n);
	}

	@Test
	public void psiModAndResid() {
		double ox = mono("PEM[Oxidation]AT");
		assertEquals(ox, mono("PEM[monohydroxylated residue]AT"), D);
		assertEquals(ox, mono("PEM[MOD:00425]AT"), D);
		assertEquals(mono("EM[Oxidation]EVEES[Phospho]PEK"), mono("EM[M:L-methionine sulfoxide]EVEES[M:O-phospho-L-serine]PEK"), D);
		assertEquals(mono("ED[Amidated]ATR"), mono("ED[MOD:00093]ATR"), D);
		assertEquals(ProFormaFormula.parse("C6H7N2O3").getMonoisotopicMass().doubleValue() + H2O, mono("X[MOD:01928]"), D);
		assertThrows("MOD:00693", new Runnable() { public void run() { mono("PEN[MOD:00693]AT"); } });
		assertEquals(586.161465, mono("EM[R:L-methionine sulfone]EM[RESID:AA0581]"), D);
		assertEquals(mono("Q[Gln->pyro-Glu]PEPTIDE"), mono("Q[RESID:AA0031]PEPTIDE"), D);
		assertEquals(mono("E[Glu->pyro-Glu]PEPTIDE"), mono("E[RESID:AA0031]PEPTIDE"), D);
		assertThrows("[RESID:AA0031]-", new Runnable() { public void run() { mono("[RESID:AA0031]-QPEPTIDE"); } });
		assertThrows("X:formaldehyde", new Runnable() { public void run() { mono("EVTK[X:formaldehyde]LEK"); } });
	}

	@Test
	public void chargeStates() {
		double m = mono("PEPTIDE");
		double proton = 1.00727646677;
		double electron = 0.0005486;
		assertEquals((m + 2 * proton) / 2, ProForma.parse("PEPTIDE/2").getMz().doubleValue(), D);
		double na = formula("Na") - electron;
		assertEquals((m + 2 * na) / 2, ProForma.parse("PEPTIDE/[Na:z+1^2]").getMz().doubleValue(), D);
		assertEquals((m + na + proton) / 2, ProForma.parse("PEPTIDE/[Na:z+1,H:z+1]").getMz().doubleValue(), D);
		assertEquals((m + 2 * electron) / 2, ProForma.parse("PEPTIDE/-2").getMz().doubleValue(), D);
		assertEquals(Integer.valueOf(4), ProForma.parse("PEPT[Formula:Zn:z+2]IDE/2").getCharge());
		assertEquals(Integer.valueOf(4), ProForma.parse("PEPT[Formula:Zn:z+2]IDE/[Na:z+1^2]").getCharge());
		assertThrows("no charge", new Runnable() { public void run() { ProForma.parse("PEPTIDE").getMz(); } });
		// ProForma is an IMassElement: it works with IonCloud and MassUtils like the other parsers
		BigDecimal mz = MassUtils.computeMz(ProForma.parse("PEPTIDE"), IonCloud.parse("2H"));
		assertEquals((m + 2 * proton) / 2, mz.doubleValue(), D);
	}

	@Test
	public void negativeCharges() {
		double m = mono("PEPTIDE");
		double proton = 1.00727646677;
		double electron = 0.0005486;
		assertEquals(-1, Compositions.charge(-1).getCharges());
		assertEquals(electron, Compositions.charge(-1).getMonoisotopicMass().doubleValue(), 1e-9);
		assertEquals(-3, Compositions.charge(-3).getCharges());
		assertEquals(3 * electron, Compositions.charge(-3).getMonoisotopicMass().doubleValue(), 1e-9);
		// [M-H]- and [M-2H]2- written with charge carriers
		assertEquals(Integer.valueOf(-1), ProForma.parse("PEPTIDE/[H-1:z-1]").getCharge());
		assertEquals(m - proton, ProForma.parse("PEPTIDE/[H-1:z-1]").getMz().doubleValue(), D);
		assertEquals(Integer.valueOf(-2), ProForma.parse("PEPTIDE/[H-1:z-1^2]").getCharge());
		assertEquals(ProForma.parse("PEPTIDE").getMz(-2).doubleValue(),
				ProForma.parse("PEPTIDE/[H-1:z-1^2]").getMz().doubleValue(), D);
		assertEquals(Integer.valueOf(-2), ProForma.parse("PEPTIDE/[e-:z-1^2]").getCharge());
		assertEquals(ProForma.parse("PEPTIDE/-2").getMz().doubleValue(),
				ProForma.parse("PEPTIDE/[e-:z-1^2]").getMz().doubleValue(), D);
		// negatively charged formula
		assertEquals(Integer.valueOf(-1), ProForma.parse("PEPT[Formula:Cl:z-1]IDE").getCharge());
		assertEquals(m + formula("Cl") + electron, mono("PEPT[Formula:Cl:z-1]IDE"), D);
	}

	@Test
	public void mzWithGivenCharge() {
		double m = mono("PEPTIDE");
		double proton = 1.00727646677;
		// protons are added for a positive charge and removed for a negative charge
		assertEquals((m + 2 * proton) / 2, ProForma.parse("PEPTIDE").getMz(2).doubleValue(), D);
		assertEquals((m - 2 * proton) / 2, ProForma.parse("PEPTIDE").getMz(-2).doubleValue(), D);
		// same as the ion representation of GlycanMassUtility (whose m/z has the sign of the charge)
		assertEquals(Math.abs(MassUtils.computeMz(ProForma.parse("PEPTIDE"), IonCloud.parse("-2H")).doubleValue()),
				ProForma.parse("PEPTIDE").getMz(-2).doubleValue(), D);
		// unlike the charge state /-2 of the string, which adds electrons
		assertTrue(ProForma.parse("PEPTIDE/-2").getMz().doubleValue() > ProForma.parse("PEPTIDE").getMz(-2).doubleValue());
		assertThrows("0", new Runnable() { public void run() { ProForma.parse("PEPTIDE").getMz(0); } });
	}

	@Test
	public void fragments() {
		double[] b = {98.06004032, 227.1026334, 324.15539725, 425.20307572, 538.2871397, 653.31408272};
		double[] y = {148.06043424, 263.08737726, 376.17144124, 477.21911971, 574.27188356, 703.31447664};
		List<BigDecimal> fb = ProForma.parse("PEPTIDE").getFragments(FragmentIonType.B, 1);
		List<BigDecimal> fy = ProForma.parse("PEPTIDE").getFragments(FragmentIonType.Y, 1);
		for ( int k = 0; k < b.length; k++ ) {
			assertEquals(b[k], fb.get(k).doubleValue(), 1e-6);
			assertEquals(y[k], fy.get(k).doubleValue(), 1e-6);
		}
		assertEquals(fb, ProForma.parse("{Hex}PEPTIDE").getFragments(FragmentIonType.B, 1));
	}

	@Test
	public void parseOrNull() {
		IMassElement e = ProForma.parseOrNull(" ELVIS[Phospho|+79.966331]K ");
		assertEquals(mono("ELVIS[Phospho]K"), e.getMonoisotopicMass().doubleValue(), D);
		assertNull(ProForma.parseOrNull("   "));
	}

	@Test
	public void equalsAndHashCode() {
		ProForma p1 = ProForma.parse("EM[Oxidation]EVN[Glycan:HexNAc2Hex3]ES[Phospho]PEK/2");
		ProForma p2 = ProForma.parse("EM[Oxidation]EVN[Glycan:HexNAc2Hex3]ES[Phospho]PEK/2");
		assertEquals(p1, p2);
		assertEquals(p1.hashCode(), p2.hashCode());
		assertEquals(p1, ProForma.parse("em[Oxidation]evn[Glycan:HexNAc2Hex3]es[Phospho]pek/2"));
		assertFalse(p1.equals(ProForma.parse("EM[Oxidation]EVN[Glycan:HexNAc2Hex3]ES[Phospho]PEK/3")));
		assertFalse(p1.equals(ProForma.parse("EMEVN[Glycan:HexNAc2Hex3]ES[Phospho]PEK/2")));

		java.util.Set<ProForma> set = new java.util.HashSet<>();
		set.add(p1);
		assertTrue(set.contains(p2));
	}

	@Test
	public void specExamplesParseAndRoundTrip() throws IOException {
		int n = 0;
		for ( String s : resourceLines("spec_examples.txt") ) {
			List<PeptidoformIon> ions = ProFormaParser.parseChimeric(s);
			String written = ProFormaWriter.writeChimeric(ions);
			assertEquals("round trip of " + s, s.toUpperCase().replace(" ", ""), written.toUpperCase().replace(" ", ""));
			assertEquals("re-parse of " + s, ions, ProFormaParser.parseChimeric(written));
			n++;
		}
		assertTrue(n > 70);
		ProFormaParser.parse("[Oxidation|Limit:2]^4?PEPTIDE");
		ProFormaParser.parse("PTI(MERMERME)[+32|Position:E]PTIDE");
	}

	@Test
	public void invalidStrings() {
		String[] bad = {"", "PEPT1DE", "PEP[Phospho", "PEP]TIDE", "PEPTIDE/", "PEP(TIDE", "[Phospho]PEPTIDE",
				"PEPTIDE-", "<99C>PEPTIDE", "PEP[Formula:Xx]TIDE", "PEPTIDE//", "PEP TIDE"};
		for ( final String s : bad )
			assertThrows(s, new Runnable() { public void run() { ProFormaParser.parseChimeric(s); } });
		assertThrows("unknown", new Runnable() { public void run() { mono("PEP[NoSuchModification]TIDE"); } });
	}
}
