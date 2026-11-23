package ninja.mspp.proformabuilder;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Parser for molecular formulas in ProForma notation (ProForma 2.1 section 7.4).
 * <ul>
 *   <li>Elements are case sensitive.</li>
 *   <li>Counts may be omitted (1) and may be negative ({@code C-6}).</li>
 *   <li>Isotopes are written in brackets: {@code [13C6]}, {@code [15N-2]}.</li>
 *   <li>{@code D} and {@code T} are accepted as shorthand for {@code 2H} and {@code 3H}.</li>
 *   <li>Spaces are allowed around elements.</li>
 * </ul>
 * Unknown elements, isotopes or characters raise {@link IllegalArgumentException}; nothing is skipped.
 */
public final class ProFormaFormula {

	private final String s;
	private int i;

	private ProFormaFormula(String s) {
		this.s = s;
	}

	public static Molecule parse(String formula) {
		return new ProFormaFormula(formula).run();
	}

	private Molecule run() {
		Molecule m = new Molecule();
		skipSpaces();
		if ( atEnd() )
			throw error("Empty formula");
		while ( !atEnd() ) {
			if ( s.charAt(i) == '[' ) {
				i++;
				skipSpaces();
				int massNumber = readUnsigned();
				if ( massNumber < 0 )
					throw error("Mass number expected in isotope");
				String symbol = readSymbol();
				int count = readCount();
				skipSpaces();
				expect(']');
				m = m.addMolecules(atom(symbol, massNumber), count);
			} else {
				String symbol = readSymbol();
				int count = readCount();
				m = m.addMolecules(atom(symbol, 0), count);
			}
			skipSpaces();
		}
		return m;
	}

	private Molecule atom(String symbol, int massNumber) {
		AtomType element = Compositions.element(symbol);
		if ( element == null || element == AtomType.D ) {
			if ( symbol.equals("D") && massNumber == 0 )
				return Compositions.atom("H", 2);
			if ( symbol.equals("T") && massNumber == 0 )
				return Compositions.atom("H", 3);
			throw error("Unknown element '" + symbol + "'");
		}
		if ( massNumber != 0 && Compositions.isotope(element, massNumber) == null )
			throw error("Unknown isotope " + massNumber + symbol);
		return Compositions.atom(symbol, massNumber);
	}

	private boolean atEnd() {
		return i >= s.length();
	}

	private void skipSpaces() {
		while ( i < s.length() && s.charAt(i) == ' ' )
			i++;
	}

	private void expect(char c) {
		if ( atEnd() || s.charAt(i) != c )
			throw error("'" + c + "' expected");
		i++;
	}

	private String readSymbol() {
		if ( atEnd() || !Character.isUpperCase(s.charAt(i)) )
			throw error("Element symbol expected");
		int start = i++;
		if ( i < s.length() && Character.isLowerCase(s.charAt(i)) )
			i++;
		return s.substring(start, i);
	}

	private int readUnsigned() {
		int start = i;
		while ( i < s.length() && Character.isDigit(s.charAt(i)) )
			i++;
		return start == i ? -1 : Integer.parseInt(s.substring(start, i));
	}

	private int readCount() {
		int sign = 1;
		if ( i < s.length() && s.charAt(i) == '-' ) {
			sign = -1;
			i++;
		}
		int n = readUnsigned();
		if ( n < 0 ) {
			if ( sign < 0 )
				throw error("Number expected after '-'");
			return 1;
		}
		return sign * n;
	}

	private IllegalArgumentException error(String msg) {
		return new IllegalArgumentException(msg + " at position " + i + " in formula '" + s + "'");
	}
}
