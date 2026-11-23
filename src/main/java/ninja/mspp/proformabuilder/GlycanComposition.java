package ninja.mspp.proformabuilder;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.glycoinfo.ms.GlycanMassUtility.dict.residue.ResidueType;
import org.glycoinfo.ms.GlycanMassUtility.dict.residue.monosaccharide.MonosaccharideDictionary;
import org.glycoinfo.ms.GlycanMassUtility.dict.residue.substituent.SubstituentDictionary;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Glycan composition of a {@code Glycan:} tag (ProForma 2.1 section 10.2).
 * <ul>
 *   <li>Monosaccharide names followed by a count; a count of 1 may be omitted; spaces are allowed
 *       around names and numbers ({@code HexNAc 2 Hex 3}).</li>
 *   <li>Names are those of the section 10.2 table and of the {@code monosaccharide_name} list of
 *       the formal grammar. Other names (e.g. {@code Sialic Acid}, {@code Fucose}) are errors.</li>
 *   <li>Components not in the list are written as a formula or a mass in braces:
 *       {@code {C13H19N1O9}1}, {@code {+203.079}1}; charged formulas {@code {C8H13N1O5Na1:z+1}1} (Level 3).</li>
 * </ul>
 * Names are matched case sensitively by longest match ({@code HexNAcS2} is HexNAcS x 2).
 * The residue composition of a name is taken from {@link MonosaccharideDictionary} or
 * {@link SubstituentDictionary}, where the ProForma names are registered as names or synonyms:
 * residue = free monosaccharide (or substituent) - H2O.
 */
public final class GlycanComposition {

	/** Names in the ProForma 2.1 section 10.2 table (26 names). */
	public static final List<String> TABLE_NAMES = Collections.unmodifiableList(Arrays.asList(
			"Hex", "HexNAc", "HexS", "HexP", "HexNAcS", "HexN", "HexNS", "dHex", "aHex", "en,aHex",
			"Neu", "NeuAc", "NeuGc", "Sug", "Tri", "Tet", "Pen", "Hep", "Oct", "Non", "Dec",
			"Fuc", "Kdo", "Kdn", "Sulfo", "Phospho"));

	/** Names that appear only in the formal grammar ({@code monosaccharide_name}). */
	public static final List<String> GRAMMAR_ONLY_NAMES = Collections.unmodifiableList(Arrays.asList(
			"Neu5Ac", "Neu5Gc", "Sulfate", "Phosphate", "S"));

	/** All accepted names. */
	public static final Set<String> NAMES;
	static {
		Set<String> s = new LinkedHashSet<>(TABLE_NAMES);
		s.addAll(GRAMMAR_ONLY_NAMES);
		NAMES = Collections.unmodifiableSet(s);
	}

	/**
	 * Names which are parsed and written but have no mass, since GlycanMassUtility does not define them:
	 * Sugose does not exist as a monosaccharide, and trioses and tetroses rarely appear in glycans.
	 */
	public static final Set<String> UNSUPPORTED_NAMES = Collections.unmodifiableSet(
			new LinkedHashSet<>(Arrays.asList("Sug", "Tri", "Tet")));

	private static final Map<String, Molecule> RESIDUE_CACHE = new ConcurrentHashMap<>();

	/** One component of the composition. */
	public abstract static class Component extends Model.Value {
		private final int count;

		Component(int count) {
			this.count = count;
		}

		public int getCount() {
			return this.count;
		}
	}

	/** Named monosaccharide, e.g. {@code HexNAc2}. */
	public static final class Named extends Component {
		private final String name;

		Named(String name, int count) {
			super(count);
			this.name = name;
		}

		public String getName() {
			return this.name;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.name, getCount()};
		}
	}

	/** Formula in braces, e.g. {@code {C13H19N1O9}1}, optionally charged. */
	public static final class FormulaComponent extends Component {
		private final Molecule molecule;
		private final String text;

		FormulaComponent(Molecule molecule, String text, int count) {
			super(count);
			this.molecule = molecule;
			this.text = text;
		}

		public Molecule getMolecule() {
			return this.molecule;
		}

		public String getText() {
			return this.text;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.text, getCount()};
		}
	}

	/** Mass in braces, e.g. {@code {+203.079}1}. */
	public static final class MassComponent extends Component {
		private final BigDecimal mass;
		private final String text;

		MassComponent(BigDecimal mass, String text, int count) {
			super(count);
			this.mass = mass;
			this.text = text;
		}

		public BigDecimal getMass() {
			return this.mass;
		}

		public String getText() {
			return this.text;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.text, getCount()};
		}
	}

	private final List<Component> components;
	private final String text;

	private GlycanComposition(List<Component> components, String text) {
		this.components = Model.copy(components);
		this.text = text;
	}

	public List<Component> getComponents() {
		return this.components;
	}

	/** Source text as written. */
	public String getText() {
		return this.text;
	}

	/** Parses the part after {@code Glycan:}. */
	public static GlycanComposition parse(String s) {
		List<Component> out = new ArrayList<>();
		int i = 0;
		int n = s.length();
		while ( true ) {
			while ( i < n && s.charAt(i) == ' ' )
				i++;
			if ( i >= n )
				break;
			if ( s.charAt(i) == '{' ) {
				int close = s.indexOf('}', i);
				if ( close < 0 )
					throw new ProFormaException("Unclosed '{' in glycan composition '" + s + "'");
				String inner = s.substring(i + 1, close).trim();
				i = close + 1;
				int[] cnt = readCount(s, i);
				i = cnt[1];
				out.add(braceComponent(inner, cnt[0]));
				continue;
			}
			String name = longestName(s, i);
			if ( name == null )
				throw new ProFormaException("Unknown monosaccharide at '" + s.substring(i)
						+ "' in glycan composition '" + s + "' (allowed: ProForma 2.1 section 10.2 names)");
			i += name.length();
			int[] cnt = readCount(s, i);
			i = cnt[1];
			out.add(new Named(name, cnt[0]));
		}
		if ( out.isEmpty() )
			throw new ProFormaException("Empty glycan composition");
		return new GlycanComposition(out, s);
	}

	private static Component braceComponent(String inner, int count) {
		if ( inner.startsWith("+") || inner.startsWith("-") ) {
			try {
				return new MassComponent(new BigDecimal(inner), "{" + inner + "}", count);
			} catch ( NumberFormatException e ) {
				throw new ProFormaException("Invalid mass in glycan composition: {" + inner + "}");
			}
		}
		try {
			return new FormulaComponent(chargedFormula(inner), "{" + inner + "}", count);
		} catch ( ProFormaException e ) {
			throw e;
		} catch ( IllegalArgumentException e ) {
			throw new ProFormaException("Invalid formula in glycan composition: {" + inner + "}: " + e.getMessage());
		}
	}

	/** Formula with an optional charge suffix {@code :z+1} (ProForma 2.1 section 11.1). */
	static Molecule chargedFormula(String text) {
		int z = text.lastIndexOf(":z");
		if ( z < 0 )
			return ProFormaFormula.parse(text);
		int charge = parseCharge(text.substring(z + 2), text);
		return Compositions.plus(ProFormaFormula.parse(text.substring(0, z)), Compositions.charge(charge));
	}

	static int parseCharge(String s, String context) {
		try {
			return Integer.parseInt(s.startsWith("+") ? s.substring(1) : s);
		} catch ( NumberFormatException e ) {
			throw new ProFormaException("Invalid charge ':z" + s + "' in '" + context + "'");
		}
	}

	private static String longestName(String s, int i) {
		String best = null;
		for ( String name : NAMES )
			if ( s.startsWith(name, i) && (best == null || name.length() > best.length()) )
				best = name;
		return best;
	}

	/** Reads an optional count (default 1) after optional spaces; returns {count, newIndex}. */
	private static int[] readCount(String s, int i) {
		int j = i;
		while ( j < s.length() && s.charAt(j) == ' ' )
			j++;
		int start = j;
		while ( j < s.length() && Character.isDigit(s.charAt(j)) )
			j++;
		if ( start == j )
			return new int[] {1, i};
		return new int[] {Integer.parseInt(s.substring(start, j)), j};
	}

	/**
	 * Residue composition (condensed form, as in the ProForma 10.2 table) of one monosaccharide name,
	 * looked up in the monosaccharide and substituent dictionaries.
	 */
	public static Molecule residue(String name) {
		Molecule m = RESIDUE_CACHE.get(name);
		if ( m != null )
			return m;
		if ( UNSUPPORTED_NAMES.contains(name) )
			throw new ProFormaException("Monosaccharide '" + name + "' is not supported");
		ResidueType type = null;
		if ( MonosaccharideDictionary.findMonosaccharide(name) )
			type = MonosaccharideDictionary.getMonosaccharide(name);
		else if ( SubstituentDictionary.findSubstituent(name) )
			type = SubstituentDictionary.getSubstituent(name);
		if ( type == null || type.getMolecule() == null )
			throw new ProFormaException("Monosaccharide '" + name + "' is not defined in the dictionaries");
		m = type.getMolecule().addMolecules(Compositions.H2O, -1);
		RESIDUE_CACHE.put(name, m);
		return m;
	}

	/** Residue mass of the whole composition (condensed form, no extra water). */
	public MassValue getMass() {
		MassValue total = MassValue.ZERO;
		for ( Component c : this.components ) {
			MassValue m;
			if ( c instanceof Named )
				m = MassValue.of(residue(((Named) c).getName()));
			else if ( c instanceof FormulaComponent )
				m = MassValue.of(((FormulaComponent) c).getMolecule());
			else
				m = MassValue.ofDelta(((MassComponent) c).getMass());
			total = total.plus(m.times(c.getCount()));
		}
		return total;
	}

	@Override
	public boolean equals(Object o) {
		return o instanceof GlycanComposition && ((GlycanComposition) o).components.equals(this.components)
				&& ((GlycanComposition) o).text.equals(this.text);
	}

	@Override
	public int hashCode() {
		return this.components.hashCode() * 31 + this.text.hashCode();
	}

	@Override
	public String toString() {
		return this.text;
	}
}
