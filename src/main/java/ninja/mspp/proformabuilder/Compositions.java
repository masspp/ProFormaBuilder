package ninja.mspp.proformabuilder;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.dict.atom.IAtom;
import org.glycoinfo.ms.GlycanMassUtility.dict.atom.IsotopeType;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Helpers for building and transforming {@link Molecule} objects (elemental compositions).
 */
public final class Compositions {

	/** Empty composition. */
	public static final Molecule EMPTY = new Molecule();
	/** Water. */
	public static final Molecule H2O = Molecule.parse("H2O");

	private static final Map<String, Molecule> ATOMS = new ConcurrentHashMap<>();

	private Compositions() { }

	/** Returns the element with the symbol (case sensitive), or {@code null}. No error is logged. */
	public static AtomType element(String symbol) {
		for ( AtomType a : AtomType.values() )
			if ( a.getSymbol().equals(symbol) && a != AtomType.e )
				return a;
		return null;
	}

	/** Returns the isotope of an element, or {@code null}. No error is logged. */
	public static IsotopeType isotope(AtomType element, int massNumber) {
		for ( IsotopeType i : IsotopeType.values() )
			if ( i.getAtomSymbol().equals(element.getSymbol()) && i.getMassNumber() == massNumber )
				return i;
		return null;
	}

	/** A molecule consisting of one atom (an element, or an isotope if {@code massNumber > 0}). */
	public static Molecule atom(String symbol, int massNumber) {
		String key = massNumber == 0 ? symbol : symbol + "^" + massNumber;
		Molecule m = ATOMS.get(key);
		if ( m == null ) {
			m = Molecule.parse(key);
			ATOMS.put(key, m);
		}
		return m;
	}

	private static Molecule atom(IAtom a) {
		if ( a instanceof IsotopeType )
			return atom(a.getAtomSymbol(), ((IsotopeType) a).getMassNumber());
		return atom(a.getSymbol(), 0);
	}

	public static Molecule plus(Molecule a, Molecule b) {
		return a.addMolecules(b, 1);
	}

	public static Molecule times(Molecule a, int n) {
		return EMPTY.addMolecules(a, n);
	}

	/**
	 * Replaces every natural (unspecified) atom of an element by the given isotope.
	 * Explicitly specified isotopes are not changed (ProForma 2.1 section 11.3.1).
	 */
	public static Molecule withGlobalIsotope(Molecule m, AtomType element, IsotopeType isotope) {
		Molecule result = new Molecule();
		for ( Map.Entry<IAtom, Integer> e : m.getRawAtomTree().entrySet() ) {
			IAtom a = e.getKey();
			Molecule one = (a == element) ? atom(isotope) : atom(a);
			result = result.addMolecules(one, e.getValue());
		}
		if ( m.getCharges() != 0 )
			result = result.addMolecules(charge(m.getCharges()), 1);
		return result;
	}

	/** A molecule that only carries a charge (its mass is the electron mass correction). */
	public static Molecule charge(int z) {
		// Molecule.parse() reads a leading '-' as the sign of the whole formula, not as a charge,
		// so a negative charge is made by multiplying a positive one
		return times(Molecule.parse("+"), z);
	}

	/**
	 * Formula in ProForma notation: natural elements with explicit counts ({@code C2H3N1O1}),
	 * then isotopes ({@code [13C6]}), then a charge suffix ({@code :z+2}) if charged.
	 */
	public static String toProFormaFormula(Molecule m) {
		StringBuilder nat = new StringBuilder();
		StringBuilder iso = new StringBuilder();
		for ( Map.Entry<IAtom, Integer> e : m.getRawAtomTree().entrySet() ) {
			IAtom a = e.getKey();
			if ( a instanceof IsotopeType )
				iso.append('[').append(((IsotopeType) a).getMassNumber()).append(a.getAtomSymbol())
				   .append(e.getValue()).append(']');
			else
				nat.append(a.getSymbol()).append(e.getValue());
		}
		String s = nat.append(iso).toString();
		if ( m.getCharges() != 0 )
			s += ":z" + (m.getCharges() > 0 ? "+" : "") + m.getCharges();
		return s;
	}
}
