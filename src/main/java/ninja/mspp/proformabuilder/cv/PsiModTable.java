package ninja.mspp.proformabuilder.cv;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;
import ninja.mspp.proformabuilder.ProFormaFormula;

/**
 * PSI-MOD terms: difference formula ({@code DiffFormula}) and mass ({@code DiffMono}), origin
 * residues, RESID names and accessions.
 * <p>
 * Non-standard amino acids have {@code DiffMono "none"}; for them the whole residue
 * ({@code Formula}, {@code MassMono}) is used, as they are meant to be written on {@code X}.
 * Some terms (e.g. MOD:00693 glycosylated residue) have no mass.
 * </p>
 * <p>
 * Read from {@code PSI-MOD.obo} (the unmodified file bundled with the library, or another file
 * given to {@link ModificationResolver}).
 * </p>
 */
public final class PsiModTable {

	/** One term. */
	public static final class Term {
		private final String id;
		private final String name;
		private final Molecule composition;
		private final BigDecimal monoMass;
		private final boolean residueSubstitute;
		private final boolean obsolete;
		private final List<String> origins;
		private final List<String> residNames;
		private final List<String> residIds;

		Term(String id, String name, Molecule composition, BigDecimal monoMass, boolean residueSubstitute,
				boolean obsolete, List<String> origins, List<String> residNames, List<String> residIds) {
			this.id = id;
			this.name = name;
			this.composition = composition;
			this.monoMass = monoMass;
			this.residueSubstitute = residueSubstitute;
			this.obsolete = obsolete;
			this.origins = Collections.unmodifiableList(origins);
			this.residNames = Collections.unmodifiableList(residNames);
			this.residIds = Collections.unmodifiableList(residIds);
		}

		public String getId() {
			return this.id;
		}

		public String getName() {
			return this.name;
		}

		/** Difference composition (whole residue for non-standard amino acids), or {@code null}. */
		public Molecule getComposition() {
			return this.composition;
		}

		/** Difference mass, or {@code null}. */
		public BigDecimal getMonoMass() {
			return this.monoMass;
		}

		/** {@code true} if the values are those of a whole residue (non-standard amino acid). */
		public boolean isResidueSubstitute() {
			return this.residueSubstitute;
		}

		public boolean isObsolete() {
			return this.obsolete;
		}

		/** Origin residues ({@code M}, {@code X}, ...). */
		public List<String> getOrigins() {
			return this.origins;
		}

		/** RESID names given as synonyms ({@code RESID-name}). */
		public List<String> getResidNames() {
			return this.residNames;
		}

		/** RESID accessions cited in the definition ({@code AA0581}). */
		public List<String> getResidIds() {
			return this.residIds;
		}

		public boolean hasMass() {
			return this.composition != null || this.monoMass != null;
		}
	}

	private static final Pattern QUOTED = Pattern.compile("\"([^\"]*)\"");
	private static final Pattern RESID_REF = Pattern.compile("RESID:(AA\\d{4})");

	private final Map<String, Term> byId = new LinkedHashMap<>();
	private final Map<String, Term> byName = new HashMap<>();
	private final Map<String, List<Term>> byLowerName = new HashMap<>();
	private final Map<String, List<Term>> byResidName = new HashMap<>();
	private final Map<String, List<Term>> byResidId = new HashMap<>();

	private PsiModTable() { }

	/** Reads {@code PSI-MOD.obo} (https://github.com/HUPO-PSI/psi-mod-CV). */
	public static PsiModTable fromObo(InputStream in) throws IOException {
		PsiModTable t = new PsiModTable();
		BufferedReader r = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
		Map<String, String> xref = new HashMap<>();
		String id = null;
		String name = null;
		String def = null;
		boolean obsolete = false;
		boolean inTerm = false;
		List<String> residNames = new ArrayList<>();
		while ( true ) {
			String line = r.readLine();
			boolean boundary = line == null || line.startsWith("[");
			if ( boundary && inTerm && id != null )
				t.add(build(id, name, def, obsolete, xref, residNames));
			if ( line == null )
				break;
			if ( boundary ) {
				inTerm = line.equals("[Term]");
				id = null;
				name = null;
				def = null;
				obsolete = false;
				xref = new HashMap<>();
				residNames = new ArrayList<>();
				continue;
			}
			if ( !inTerm )
				continue;
			if ( line.startsWith("id: ") ) {
				id = line.substring(4).trim();
			} else if ( line.startsWith("name: ") ) {
				name = line.substring(6).trim();
			} else if ( line.startsWith("def: ") ) {
				def = line;
			} else if ( line.startsWith("is_obsolete: true") ) {
				obsolete = true;
			} else if ( line.startsWith("synonym: ") && line.contains(" RESID-name ") ) {
				Matcher m = QUOTED.matcher(line);
				if ( m.find() )
					residNames.add(m.group(1));
			} else if ( line.startsWith("xref: ") ) {
				int colon = line.indexOf(':', 6);
				if ( colon > 0 ) {
					Matcher m = QUOTED.matcher(line.substring(colon));
					if ( m.find() )
						xref.put(line.substring(6, colon).trim(), m.group(1));
				}
			}
		}
		return t;
	}

	private static Term build(String id, String name, String def, boolean obsolete, Map<String, String> xref,
			List<String> residNames) {
		boolean substitute = "none".equals(xref.get("DiffMono"));
		String formula = substitute ? xref.get("Formula") : xref.get("DiffFormula");
		String mass = substitute ? xref.get("MassMono") : xref.get("DiffMono");
		Molecule comp = null;
		if ( formula != null && !formula.equals("none") ) {
			try {
				comp = parseOboFormula(formula);
			} catch ( IllegalArgumentException e ) {
				comp = null; // unknown element: fall back to the mass
			}
		}
		BigDecimal mono = null;
		if ( mass != null ) {
			try {
				mono = new BigDecimal(mass);
			} catch ( NumberFormatException e ) {
				mono = null; // "none"
			}
		}
		List<String> residIds = new ArrayList<>();
		if ( def != null ) {
			Matcher m = RESID_REF.matcher(def);
			while ( m.find() )
				residIds.add(m.group(1));
		}
		List<String> origins = new ArrayList<>();
		String origin = xref.get("Origin");
		if ( origin != null )
			for ( String o : origin.split(",") )
				origins.add(o.trim());
		return new Term(id, name, comp, mono, substitute, obsolete, origins, residNames, residIds);
	}

	/** PSI-MOD formula notation: {@code "C 2 H 3 N 1 O 1"}, isotopes {@code "(13)C 6"}, negative counts {@code "H -2"}. */
	static Molecule parseOboFormula(String f) {
		String[] t = f.trim().split("\\s+");
		if ( t.length % 2 != 0 )
			throw new IllegalArgumentException("Invalid PSI-MOD formula: " + f);
		Molecule m = new Molecule();
		for ( int k = 0; k < t.length; k += 2 ) {
			String el = t[k];
			int massNumber = 0;
			if ( el.startsWith("(") ) {
				int close = el.indexOf(')');
				massNumber = Integer.parseInt(el.substring(1, close));
				el = el.substring(close + 1);
			}
			Molecule atom = ProFormaFormula.parse(massNumber == 0 ? el : "[" + massNumber + el + "]");
			m = m.addMolecules(atom, Integer.parseInt(t[k + 1]));
		}
		return m;
	}

	private void add(Term t) {
		this.byId.put(t.id, t);
		if ( t.name != null ) {
			if ( !t.obsolete && !this.byName.containsKey(t.name) )
				this.byName.put(t.name, t);
			put(this.byLowerName, t.name.toLowerCase(Locale.ROOT), t);
		}
		if ( !t.obsolete ) {
			for ( String n : t.residNames )
				put(this.byResidName, n, t);
			for ( String rid : t.residIds )
				put(this.byResidId, rid, t);
		}
	}

	private static void put(Map<String, List<Term>> map, String key, Term t) {
		List<Term> l = map.get(key);
		if ( l == null ) {
			l = new ArrayList<>();
			map.put(key, l);
		}
		l.add(t);
	}

	/** Term by accession number ({@code "00425"} for MOD:00425), or {@code null}. */
	public Term byAccession(String digits) {
		try {
			return this.byId.get(String.format("MOD:%05d", Integer.parseInt(digits)));
		} catch ( NumberFormatException e ) {
			return null;
		}
	}

	/** Term by name (exact, then unique case-insensitive among non-obsolete terms), or {@code null}. */
	public Term byName(String name) {
		Term t = this.byName.get(name);
		if ( t != null )
			return t;
		List<Term> l = this.byLowerName.get(name.toLowerCase(Locale.ROOT));
		if ( l == null )
			return null;
		Term found = null;
		for ( Term x : l ) {
			if ( x.obsolete )
				continue;
			if ( found != null )
				return null;
			found = x;
		}
		return found;
	}

	/** Non-obsolete terms carrying the RESID name as a synonym. */
	public List<Term> byResidName(String name) {
		List<Term> l = this.byResidName.get(name);
		return l == null ? Collections.<Term>emptyList() : l;
	}

	/** Non-obsolete terms citing the RESID accession ({@code AA0581}) in their definition. */
	public List<Term> byResidId(String id) {
		List<Term> l = this.byResidId.get(id);
		return l == null ? Collections.<Term>emptyList() : l;
	}

	public int size() {
		return this.byId.size();
	}

	public Map<String, Term> getTerms() {
		return Collections.unmodifiableMap(this.byId);
	}
}
