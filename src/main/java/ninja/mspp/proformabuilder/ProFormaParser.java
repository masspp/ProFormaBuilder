package ninja.mspp.proformabuilder;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.dict.atom.IsotopeType;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;
import ninja.mspp.proformabuilder.Descriptor.Cv;
import ninja.mspp.proformabuilder.Model.ChargeCarrier;
import ninja.mspp.proformabuilder.Model.ChargeState;
import ninja.mspp.proformabuilder.Model.FixedModification;
import ninja.mspp.proformabuilder.Model.GlobalIsotope;
import ninja.mspp.proformabuilder.Model.Group;
import ninja.mspp.proformabuilder.Model.Modification;
import ninja.mspp.proformabuilder.Model.Peptidoform;
import ninja.mspp.proformabuilder.Model.PeptidoformIon;
import ninja.mspp.proformabuilder.Model.Range;
import ninja.mspp.proformabuilder.Model.Residue;
import ninja.mspp.proformabuilder.Model.SequenceAmbiguity;
import ninja.mspp.proformabuilder.Model.Unlocalized;

/**
 * Hand-written recursive-descent parser for ProForma 2.1.
 * Invalid input raises {@link ProFormaException} with the position; nothing is skipped silently.
 */
public final class ProFormaParser {

	private final String s;
	private int i;

	private ProFormaParser(String s) {
		this.s = s;
	}

	/** Parses one peptidoform ion. Chimeric strings ({@code A/2+B/3}) are rejected; use {@link #parseChimeric}. */
	public static PeptidoformIon parse(String s) {
		List<PeptidoformIon> ions = parseChimeric(s);
		if ( ions.size() != 1 )
			throw new ProFormaException("Chimeric ProForma string ('+' between peptidoform ions); use parseChimeric");
		return ions.get(0);
	}

	/** Parses one or more peptidoform ions separated by {@code +} (ProForma 2.1 section 11.4). */
	public static List<PeptidoformIon> parseChimeric(String s) {
		if ( s == null || s.isEmpty() )
			throw new ProFormaException("Empty ProForma string");
		if ( s.indexOf('\n') >= 0 || s.indexOf('\r') >= 0 )
			throw new ProFormaException("Line breaks are not allowed in ProForma");
		ProFormaParser p = new ProFormaParser(s);
		List<PeptidoformIon> ions = new ArrayList<>();
		ions.add(p.ion());
		while ( !p.atEnd() ) {
			p.expect('+');
			ions.add(p.ion());
		}
		return ions;
	}

	// ------------------------------------------------------------------ ion level

	private PeptidoformIon ion() {
		List<GlobalIsotope> isotopes = new ArrayList<>();
		List<FixedModification> fixed = new ArrayList<>();
		while ( peekIs('<') ) {
			if ( peekIs(1, '[') )
				fixed.add(fixedModification());
			else
				isotopes.add(globalIsotope());
		}
		String ionName = null;
		if ( s.startsWith("(>>", i) && !s.startsWith("(>>>", i) ) {
			i += 3;
			ionName = readUntil(')');
			expect(')');
		}
		List<Peptidoform> chains = new ArrayList<>();
		chains.add(peptidoform());
		while ( s.startsWith("//", i) ) {
			i += 2;
			chains.add(peptidoform());
		}
		ChargeState charge = null;
		if ( peekIs('/') )
			charge = chargeState();
		return new PeptidoformIon(ionName, isotopes, fixed, chains, charge);
	}

	private GlobalIsotope globalIsotope() {
		int start = i;
		expect('<');
		String text = readUntil('>');
		expect('>');
		int j = 0;
		while ( j < text.length() && Character.isDigit(text.charAt(j)) )
			j++;
		String symbol = text.substring(j);
		int massNumber = j == 0 ? 0 : Integer.parseInt(text.substring(0, j));
		if ( massNumber == 0 && symbol.equals("D") ) {
			symbol = "H";
			massNumber = 2;
		} else if ( massNumber == 0 && symbol.equals("T") ) {
			symbol = "H";
			massNumber = 3;
		}
		AtomType element = massNumber == 0 ? null : Compositions.element(symbol);
		IsotopeType isotope = element == null ? null : Compositions.isotope(element, massNumber);
		if ( isotope == null )
			throw new ProFormaException("Invalid global isotope <" + text + ">", start);
		return new GlobalIsotope(element, isotope, text);
	}

	private FixedModification fixedModification() {
		expect('<');
		Modification mod = bracketModification();
		expect('@');
		String targetsText = readUntil('>');
		expect('>');
		List<String> targets = new ArrayList<>();
		for ( String t : targetsText.split(",") ) {
			String target = t.trim();
			String norm = normaliseTarget(target);
			if ( norm == null )
				throw new ProFormaException("Invalid fixed modification target '" + target + "'", i);
			targets.add(norm);
		}
		return new FixedModification(mod, targets);
	}

	private static String normaliseTarget(String t) {
		String u = t.toUpperCase(Locale.ROOT);
		if ( u.length() == 1 && AminoAcids.isValid(u.charAt(0)) )
			return u;
		if ( u.equals("N-TERM") )
			return "N-term";
		if ( u.equals("C-TERM") )
			return "C-term";
		if ( (u.startsWith("N-TERM:") || u.startsWith("C-TERM:")) && u.length() == 8 && AminoAcids.isValid(u.charAt(7)) )
			return (u.charAt(0) == 'N' ? "N-term:" : "C-term:") + u.charAt(7);
		return null;
	}

	private ChargeState chargeState() {
		expect('/');
		if ( peekIs('[') ) {
			i++;
			List<ChargeCarrier> carriers = new ArrayList<>();
			while ( true ) {
				int start = i;
				int depth = 0;
				while ( !atEnd() ) {
					char c = s.charAt(i);
					if ( c == '[' )
						depth++;
					else if ( c == ']' ) {
						if ( depth == 0 )
							break;
						depth--;
					} else if ( c == ',' && depth == 0 )
						break;
					i++;
				}
				carriers.add(chargeCarrier(s.substring(start, i).trim(), start));
				if ( peekIs(',') ) {
					i++;
					continue;
				}
				expect(']');
				break;
			}
			return new ChargeState(0, carriers);
		}
		int start = i;
		if ( peekIs('+') || peekIs('-') )
			i++;
		while ( !atEnd() && Character.isDigit(s.charAt(i)) )
			i++;
		String num = s.substring(start, i);
		if ( num.isEmpty() || num.equals("+") || num.equals("-") )
			throw new ProFormaException("Charge expected after '/'", start);
		return new ChargeState(Integer.parseInt(num.startsWith("+") ? num.substring(1) : num),
				new ArrayList<ChargeCarrier>());
	}

	private static ChargeCarrier chargeCarrier(String text, int pos) {
		int count = 1;
		String body = text;
		int caret = text.lastIndexOf('^');
		if ( caret >= 0 ) {
			try {
				count = Integer.parseInt(text.substring(caret + 1).trim());
			} catch ( NumberFormatException e ) {
				throw new ProFormaException("Invalid charge carrier count: " + text, pos);
			}
			body = text.substring(0, caret);
		}
		int z = body.indexOf(":z");
		if ( z < 0 )
			throw new ProFormaException("Charge carrier needs ':z<charge>': " + text, pos);
		int charge = GlycanComposition.parseCharge(body.substring(z + 2).trim(), text);
		String formula = body.substring(0, z).trim();
		Molecule m;
		if ( formula.equals("e-") || formula.equals("e") ) {
			m = Compositions.charge(charge); // electron carrier: mass handled by the charge
		} else {
			try {
				m = Compositions.plus(ProFormaFormula.parse(formula), Compositions.charge(charge));
			} catch ( IllegalArgumentException e ) {
				throw new ProFormaException("Invalid charge carrier '" + text + "': " + e.getMessage(), pos);
			}
		}
		return new ChargeCarrier(m, count, text);
	}

	// ------------------------------------------------------------------ peptidoform level

	private Peptidoform peptidoform() {
		String name = null;
		if ( s.startsWith("(>", i) && !s.startsWith("(>>", i) ) {
			i += 2;
			name = readUntil(')');
			expect(')');
		}
		List<Modification> labile = new ArrayList<>();
		List<Unlocalized> unlocalized = new ArrayList<>();
		List<Modification> nTerm = new ArrayList<>();

		// prefix: {labile} ... [unlocalized]^n? ... [N-term]-
		while ( true ) {
			if ( peekIs('{') ) {
				labile.add(braceModification());
				continue;
			}
			if ( peekIs('[') ) {
				List<Unlocalized> pending = new ArrayList<>();
				int start = i;
				while ( peekIs('[') ) {
					Modification m = bracketModification();
					int count = 1;
					if ( peekIs('^') ) {
						i++;
						count = readPositiveInt();
					}
					pending.add(new Unlocalized(m, count));
				}
				if ( peekIs('?') ) {
					i++;
					unlocalized.addAll(pending);
					continue;
				}
				if ( peekIs('-') ) {
					i++;
					for ( Unlocalized u : pending ) {
						if ( u.getCount() != 1 )
							throw new ProFormaException("'^' is not allowed on N-terminal modifications", start);
						nTerm.add(u.getModification());
					}
					break;
				}
				throw new ProFormaException("'?' or '-' expected after modifications before the sequence", i);
			}
			break;
		}

		List<Residue> residues = new ArrayList<>();
		List<Range> ranges = new ArrayList<>();
		List<SequenceAmbiguity> ambiguities = new ArrayList<>();
		int rangeStart = -1;
		boolean rangeIsAmbiguity = false;
		while ( !atEnd() ) {
			char c = s.charAt(i);
			if ( c == '(' ) {
				if ( rangeStart >= 0 )
					throw new ProFormaException("Nested or overlapping ranges are not allowed", i);
				i++;
				rangeIsAmbiguity = peekIs('?');
				if ( rangeIsAmbiguity )
					i++;
				rangeStart = residues.size();
				continue;
			}
			if ( c == ')' ) {
				if ( rangeStart < 0 )
					throw new ProFormaException("Unmatched ')'", i);
				i++;
				if ( residues.size() == rangeStart )
					throw new ProFormaException("Empty range", i);
				if ( rangeIsAmbiguity ) {
					ambiguities.add(new SequenceAmbiguity(rangeStart, residues.size()));
				} else {
					List<Modification> mods = new ArrayList<>();
					while ( peekIs('[') )
						mods.add(bracketModification());
					ranges.add(new Range(rangeStart, residues.size(), mods));
				}
				rangeStart = -1;
				continue;
			}
			if ( Character.isLetter(c) ) {
				char aa = Character.toUpperCase(c);
				if ( !AminoAcids.isValid(aa) )
					throw new ProFormaException("Invalid amino acid '" + c + "'", i);
				i++;
				List<Modification> mods = new ArrayList<>();
				while ( peekIs('[') )
					mods.add(bracketModification());
				residues.add(new Residue(aa, mods));
				continue;
			}
			break;
		}
		if ( rangeStart >= 0 )
			throw new ProFormaException("Unclosed '('", i);
		if ( residues.isEmpty() )
			throw new ProFormaException("Empty sequence", i);

		List<Modification> cTerm = new ArrayList<>();
		if ( peekIs('-') && peekIs(1, '[') ) {
			i++;
			while ( peekIs('[') )
				cTerm.add(bracketModification());
		}
		if ( !atEnd() && !peekIs('/') && !peekIs('+') && !s.startsWith("//", i) )
			throw new ProFormaException("Unexpected character '" + s.charAt(i) + "'", i);
		return new Peptidoform(name, labile, unlocalized, nTerm, residues, cTerm, ranges, ambiguities);
	}

	// ------------------------------------------------------------------ modifications

	private Modification bracketModification() {
		int start = i;
		expect('[');
		return modification(readBalanced('[', ']', start), start);
	}

	private Modification braceModification() {
		int start = i;
		expect('{');
		return modification(readBalanced('{', '}', start), start);
	}

	/** Reads up to the matching closing bracket (nested brackets allowed, e.g. {@code [Formula:[13C6]H]}). */
	private String readBalanced(char open, char close, int start) {
		int depth = 0;
		int begin = i;
		while ( !atEnd() ) {
			char c = s.charAt(i);
			if ( c == open ) {
				depth++;
			} else if ( c == close ) {
				if ( depth == 0 ) {
					String content = s.substring(begin, i);
					i++;
					return content;
				}
				depth--;
			}
			i++;
		}
		throw new ProFormaException("Unclosed '" + open + "'", start);
	}

	private static Modification modification(String content, int pos) {
		List<String> parts = splitTopLevel(content, '|');
		List<Descriptor> descriptors = new ArrayList<>();
		Group group = null;
		for ( String part : parts ) {
			String text = part;
			// group suffix (#g1, #g1(0.9), #XL1, #BRANCH) on any non-INFO part
			if ( !startsWithIgnoreCase(text, "INFO:") ) {
				int hash = text.lastIndexOf('#');
				if ( hash >= 0 ) {
					if ( group != null )
						throw new ProFormaException("More than one group label in a tag: [" + content + "]", pos);
					group = group(text.substring(hash + 1), pos);
					text = text.substring(0, hash);
					if ( text.isEmpty() )
						continue;
				}
			}
			if ( text.isEmpty() )
				throw new ProFormaException("Empty modification descriptor in [" + content + "]", pos);
			descriptors.add(descriptor(text, pos));
		}
		if ( descriptors.isEmpty() && group == null )
			throw new ProFormaException("Empty modification", pos);
		return new Modification(descriptors, group);
	}

	private static Group group(String text, int pos) {
		Double score = null;
		String scoreText = null;
		String label = text;
		int paren = text.indexOf('(');
		if ( paren >= 0 ) {
			if ( !text.endsWith(")") )
				throw new ProFormaException("Invalid localisation score in #" + text, pos);
			label = text.substring(0, paren);
			scoreText = text.substring(paren + 1, text.length() - 1);
			try {
				score = Double.valueOf(scoreText);
			} catch ( NumberFormatException e ) {
				throw new ProFormaException("Invalid localisation score in #" + text, pos);
			}
		}
		if ( label.isEmpty() )
			throw new ProFormaException("Invalid group label #" + text, pos);
		for ( char ch : label.toCharArray() )
			if ( !Character.isLetterOrDigit(ch) && ch != '_' )
				throw new ProFormaException("Invalid group label #" + text, pos);
		return new Group(label, score, scoreText);
	}

	private static Descriptor descriptor(String text, int pos) {
		int colon = text.indexOf(':');
		if ( colon > 0 ) {
			String prefix = text.substring(0, colon);
			String rest = text.substring(colon + 1);
			String upper = prefix.toUpperCase(Locale.ROOT);
			if ( upper.equals("INFO") )
				return new Descriptor.Info(rest, text);
			if ( upper.equals("FORMULA") ) {
				try {
					return new Descriptor.FormulaTag(GlycanComposition.chargedFormula(rest), text);
				} catch ( ProFormaException e ) {
					throw e;
				} catch ( IllegalArgumentException e ) {
					throw new ProFormaException("Invalid formula '" + text + "': " + e.getMessage(), pos);
				}
			}
			if ( upper.equals("GLYCAN") )
				return new Descriptor.GlycanTag(GlycanComposition.parse(rest), text);
			if ( upper.equals("POSITION") || upper.equals("LIMIT") )
				return new Descriptor.Placement(prefix, rest, text);
			Cv accession = Cv.byAccessionPrefix(prefix);
			if ( accession != null && !rest.isEmpty() && !isSignedNumber(rest) )
				return new Descriptor.Accession(accession, rest, text);
			Cv cv = Cv.byShortPrefix(prefix);
			if ( cv != null ) {
				if ( isSignedNumber(rest) )
					return new Descriptor.MassDelta(cv, new BigDecimal(rest), text);
				if ( cv == Cv.OBSERVED )
					throw new ProFormaException("Obs: needs a mass delta: " + text, pos);
				return new Descriptor.Named(cv, rest, text);
			}
		}
		if ( text.equalsIgnoreCase("CoMKP") || text.equalsIgnoreCase("CoMUP")
				|| text.equalsIgnoreCase("ColocaliseModificationsOfKnownPosition")
				|| text.equalsIgnoreCase("ColocaliseModificationsOfUnknownPosition") )
			return new Descriptor.Placement(text, null, text);
		if ( isSignedNumber(text) )
			return new Descriptor.MassDelta(null, new BigDecimal(text), text);
		return new Descriptor.Named(null, text, text);
	}

	private static boolean isSignedNumber(String t) {
		if ( t.length() < 2 || (t.charAt(0) != '+' && t.charAt(0) != '-') )
			return false;
		try {
			new BigDecimal(t);
			return true;
		} catch ( NumberFormatException e ) {
			return false;
		}
	}

	private static List<String> splitTopLevel(String s, char sep) {
		List<String> out = new ArrayList<>();
		int depth = 0;
		int start = 0;
		for ( int k = 0; k < s.length(); k++ ) {
			char c = s.charAt(k);
			if ( c == '[' || c == '{' || c == '(' )
				depth++;
			else if ( c == ']' || c == '}' || c == ')' )
				depth--;
			else if ( c == sep && depth == 0 ) {
				out.add(s.substring(start, k));
				start = k + 1;
			}
		}
		out.add(s.substring(start));
		return out;
	}

	private static boolean startsWithIgnoreCase(String s, String prefix) {
		return s.regionMatches(true, 0, prefix, 0, prefix.length());
	}

	// ------------------------------------------------------------------ low level

	private boolean atEnd() {
		return i >= s.length();
	}

	private boolean peekIs(char c) {
		return i < s.length() && s.charAt(i) == c;
	}

	private boolean peekIs(int offset, char c) {
		return i + offset < s.length() && s.charAt(i + offset) == c;
	}

	private void expect(char c) {
		if ( !peekIs(c) )
			throw new ProFormaException("'" + c + "' expected" + (atEnd() ? " but string ended"
					: " but found '" + s.charAt(i) + "'"), i);
		i++;
	}

	private String readUntil(char c) {
		int j = s.indexOf(c, i);
		if ( j < 0 )
			throw new ProFormaException("'" + c + "' expected", i);
		String r = s.substring(i, j);
		i = j;
		return r;
	}

	private int readPositiveInt() {
		int start = i;
		while ( !atEnd() && Character.isDigit(s.charAt(i)) )
			i++;
		if ( start == i )
			throw new ProFormaException("Number expected", start);
		return Integer.parseInt(s.substring(start, i));
	}
}
