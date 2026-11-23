package ninja.mspp.proformabuilder.cv;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ninja.mspp.proformabuilder.Descriptor;
import ninja.mspp.proformabuilder.Descriptor.Cv;
import ninja.mspp.proformabuilder.MassValue;
import ninja.mspp.proformabuilder.ProFormaException;

/**
 * Turns a modification descriptor into a mass contribution.
 * <ul>
 *   <li>Unimod names and accessions ({@code Phospho}, {@code U:Phospho}, {@code UNIMOD:21}): delta composition.</li>
 *   <li>PSI-MOD names and accessions ({@code M:...}, {@code MOD:00046}): difference composition.</li>
 *   <li>Names without prefix: Unimod first, then PSI-MOD (ProForma 2.1 section 6.2.1).</li>
 *   <li>RESID names and accessions ({@code R:...}, {@code RESID:AA0581}): through PSI-MOD (RESID names are
 *       PSI-MOD synonyms, RESID accessions are cited in PSI-MOD definitions), narrowed down by the residue
 *       carrying the modification; if the mass is still not unique an exception is raised.</li>
 *   <li>Mass deltas: mass only (average mass approximate). Formulas: composition (charged formulas carry a charge).
 *       Glycans: residue compositions from the monosaccharide dictionaries.</li>
 *   <li>XL-MOD, GNOme and BRNO names/accessions: not supported, an exception is raised.</li>
 * </ul>
 * The default resolver reads the unmodified {@code unimod.xml} and {@code PSI-MOD.obo} bundled with the
 * library (their source, version and retrieval date are listed in {@code CV_SOURCES.md} next to them);
 * the system properties {@value #UNIMOD_PROPERTY} and {@value #PSIMOD_PROPERTY} may name other files
 * to use instead (e.g. newer releases). The files are read once, when a name is first resolved.
 */
public final class ModificationResolver {

	/** System property naming an unimod.xml file to use instead of the bundled one. */
	public static final String UNIMOD_PROPERTY = "proformabuilder.unimod.xml";
	/** System property naming a PSI-MOD.obo file to use instead of the bundled one. */
	public static final String PSIMOD_PROPERTY = "proformabuilder.psimod.obo";

	private static volatile ModificationResolver defaultInstance;

	private final UnimodTable unimod;
	private final PsiModTable psimod;

	public ModificationResolver(UnimodTable unimod, PsiModTable psimod) {
		this.unimod = unimod;
		this.psimod = psimod;
	}

	/** Default resolver: bundled CV files, or the files named by the system properties. */
	public static ModificationResolver getDefault() {
		ModificationResolver r = defaultInstance;
		if ( r == null ) {
			synchronized ( ModificationResolver.class ) {
				if ( defaultInstance == null )
					defaultInstance = loadDefault();
				r = defaultInstance;
			}
		}
		return r;
	}

	/** Replaces the default resolver. */
	public static void setDefault(ModificationResolver resolver) {
		defaultInstance = resolver;
	}

	/** Creates a resolver from an unimod.xml and a PSI-MOD.obo file. */
	public static ModificationResolver fromFiles(Path unimodXml, Path psiModObo) throws IOException {
		try ( InputStream u = Files.newInputStream(unimodXml); InputStream p = Files.newInputStream(psiModObo) ) {
			return new ModificationResolver(UnimodTable.fromXml(u), PsiModTable.fromObo(p));
		}
	}

	private static ModificationResolver loadDefault() {
		try {
			UnimodTable u;
			String up = System.getProperty(UNIMOD_PROPERTY);
			if ( up != null ) {
				try ( InputStream in = Files.newInputStream(Paths.get(up)) ) {
					u = UnimodTable.fromXml(in);
				}
			} else {
				try ( InputStream in = resource("unimod.xml") ) {
					u = UnimodTable.fromXml(in);
				}
			}
			PsiModTable p;
			String pp = System.getProperty(PSIMOD_PROPERTY);
			if ( pp != null ) {
				try ( InputStream in = Files.newInputStream(Paths.get(pp)) ) {
					p = PsiModTable.fromObo(in);
				}
			} else {
				try ( InputStream in = resource("PSI-MOD.obo") ) {
					p = PsiModTable.fromObo(in);
				}
			}
			return new ModificationResolver(u, p);
		} catch ( IOException e ) {
			throw new UncheckedIOException(e);
		}
	}

	private static InputStream resource(String name) throws IOException {
		InputStream in = ModificationResolver.class.getResourceAsStream(name);
		if ( in == null )
			throw new IOException("Bundled CV file not found: " + name);
		return in;
	}

	public UnimodTable getUnimod() {
		return this.unimod;
	}

	public PsiModTable getPsiMod() {
		return this.psimod;
	}

	/** {@code true} if the descriptor contributes mass ({@code INFO} and placement rules do not). */
	public static boolean hasMass(Descriptor d) {
		return !(d instanceof Descriptor.Info) && !(d instanceof Descriptor.Placement);
	}

	/**
	 * Mass contribution of one descriptor.
	 * @param residue the residue carrying the modification (upper case), or {@code null}
	 *                (termini, labile, unknown position); used to narrow down RESID entries
	 */
	public MassValue resolve(Descriptor d, Character residue) {
		if ( d instanceof Descriptor.MassDelta )
			return MassValue.ofDelta(((Descriptor.MassDelta) d).getValue());
		if ( d instanceof Descriptor.FormulaTag )
			return MassValue.of(((Descriptor.FormulaTag) d).getMolecule());
		if ( d instanceof Descriptor.GlycanTag ) {
			try {
				return ((Descriptor.GlycanTag) d).getGlycan().getMass();
			} catch ( ProFormaException e ) {
				throw e;
			}
		}
		if ( d instanceof Descriptor.Named )
			return named((Descriptor.Named) d, residue);
		if ( d instanceof Descriptor.Accession )
			return accession((Descriptor.Accession) d, residue);
		return MassValue.ZERO;
	}

	private MassValue named(Descriptor.Named n, Character residue) {
		Cv cv = n.getCv();
		if ( cv == null || cv == Cv.UNIMOD ) {
			UnimodTable.Entry e = this.unimod.byName(n.getName());
			if ( e != null )
				return MassValue.of(e.getComposition());
			if ( cv == Cv.UNIMOD )
				throw new ProFormaException("Unknown Unimod name '" + n.getName() + "'");
		}
		if ( cv == null || cv == Cv.PSI_MOD ) {
			PsiModTable.Term t = this.psimod.byName(n.getName());
			if ( t != null )
				return psiModMass(t, n.getText());
			throw new ProFormaException("Unknown modification name '" + n.getName() + "'"
					+ (cv == null ? " (not found in Unimod or PSI-MOD)" : " in PSI-MOD"));
		}
		if ( cv == Cv.RESID ) {
			String name = n.getName().trim();
			List<PsiModTable.Term> candidates = new ArrayList<>(this.psimod.byResidName(name));
			PsiModTable.Term byName = this.psimod.byName(name);
			if ( byName != null && !candidates.contains(byName) )
				candidates.add(byName);
			return resid(candidates, residue, n.getText());
		}
		throw new ProFormaException("Modifications from " + cv + " are not supported: " + n.getText());
	}

	private MassValue accession(Descriptor.Accession a, Character residue) {
		Cv cv = a.getCv();
		if ( cv == Cv.UNIMOD ) {
			UnimodTable.Entry e;
			try {
				e = this.unimod.byId(Integer.parseInt(a.getId().trim()));
			} catch ( NumberFormatException ex ) {
				throw new ProFormaException("Invalid Unimod accession: " + a.getText());
			}
			if ( e == null )
				throw new ProFormaException("Unknown Unimod accession: " + a.getText());
			return MassValue.of(e.getComposition());
		}
		if ( cv == Cv.PSI_MOD ) {
			PsiModTable.Term t = this.psimod.byAccession(a.getId().trim());
			if ( t == null )
				throw new ProFormaException("Unknown PSI-MOD accession: " + a.getText());
			return psiModMass(t, a.getText());
		}
		if ( cv == Cv.RESID )
			return resid(this.psimod.byResidId(a.getId().trim()), residue, a.getText());
		throw new ProFormaException("Accessions from " + cv + " are not supported: " + a.getText());
	}

	/** RESID through PSI-MOD: prefer terms whose origin is the residue, then origin X; the mass must be unique. */
	private static MassValue resid(List<PsiModTable.Term> candidates, Character residue, String text) {
		List<PsiModTable.Term> withMass = new ArrayList<>();
		for ( PsiModTable.Term t : candidates )
			if ( !t.isObsolete() && t.hasMass() )
				withMass.add(t);
		if ( withMass.isEmpty() )
			throw new ProFormaException((candidates.isEmpty() ? "Unknown RESID entry (no PSI-MOD term refers to it): "
					: "RESID entry has no defined mass in PSI-MOD: ") + text);
		List<PsiModTable.Term> chosen = withMass;
		if ( residue != null ) {
			List<PsiModTable.Term> exact = new ArrayList<>();
			List<PsiModTable.Term> any = new ArrayList<>();
			for ( PsiModTable.Term t : withMass ) {
				if ( t.getOrigins().contains(residue.toString()) )
					exact.add(t);
				if ( t.getOrigins().contains("X") )
					any.add(t);
			}
			chosen = !exact.isEmpty() ? exact : !any.isEmpty() ? any : withMass;
		}
		Map<String, PsiModTable.Term> distinct = new LinkedHashMap<>();
		for ( PsiModTable.Term t : chosen ) {
			String key = t.getComposition() != null ? t.getComposition().toString() : t.getMonoMass().toPlainString();
			if ( !distinct.containsKey(key) )
				distinct.put(key, t);
		}
		if ( distinct.size() != 1 ) {
			List<String> ids = new ArrayList<>();
			for ( PsiModTable.Term t : distinct.values() )
				ids.add(t.getId() + " (origin " + t.getOrigins() + ")");
			throw new ProFormaException("RESID entry " + text + " has more than one mass"
					+ (residue == null ? "" : " on residue " + residue) + ": " + ids);
		}
		return psiModMass(distinct.values().iterator().next(), text);
	}

	private static MassValue psiModMass(PsiModTable.Term t, String text) {
		if ( t.isObsolete() )
			throw new ProFormaException("PSI-MOD term " + t.getId() + " is obsolete: " + text);
		if ( t.getComposition() != null )
			return MassValue.of(t.getComposition());
		if ( t.getMonoMass() != null )
			return MassValue.ofDelta(t.getMonoMass());
		throw new ProFormaException("PSI-MOD term " + t.getId() + " (" + t.getName() + ") has no defined mass: " + text);
	}
}
