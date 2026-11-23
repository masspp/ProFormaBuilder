package ninja.mspp.proformabuilder;

import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Peptide fragment ion series. The neutral fragment is
 * {@code shift + H2O + sum(residues and modifications)} (same compositions as pyteomics' {@code std_ion_comp}).
 */
public enum FragmentIonType {
	A("a", "H-2O-2C-1", false),
	B("b", "H-2O-1", false),
	C("c", "H1O-1N1", false),
	X("x", "H-2O1C1", true),
	Y("y", null, true),
	Z("z", "H-3N-1", true),
	C_DOT("c-dot", "H2O-1N1", false),
	C_MINUS_1("c-1", "O-1N1", false),
	Z_DOT("z-dot", "H-2N-1", true),
	Z_PLUS_1("z+1", "H-2N-1", true),
	Z_PLUS_2("z+2", "H-1N-1", true);

	private final String label;
	private final String shiftFormula;
	private final boolean cTerminal;
	private Molecule shift;

	FragmentIonType(String label, String shiftFormula, boolean cTerminal) {
		this.label = label;
		this.shiftFormula = shiftFormula;
		this.cTerminal = cTerminal;
	}

	public String getLabel() {
		return this.label;
	}

	/** Composition shift relative to {@code H2O + residues}. */
	public synchronized Molecule getShift() {
		if ( this.shift == null )
			this.shift = this.shiftFormula == null ? Compositions.EMPTY : ProFormaFormula.parse(this.shiftFormula);
		return this.shift;
	}

	/** {@code true} for x, y, z series (counted from the C-terminus). */
	public boolean isCTerminal() {
		return this.cTerminal;
	}

	/** Looks up by label ({@code "b"}, {@code "z-dot"}). */
	public static FragmentIonType forLabel(String label) {
		for ( FragmentIonType t : values() )
			if ( t.label.equals(label) )
				return t;
		throw new IllegalArgumentException("Unknown ion type: " + label);
	}
}
