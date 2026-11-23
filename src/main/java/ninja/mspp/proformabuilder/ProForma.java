package ninja.mspp.proformabuilder;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.dict.molecule.IonType;
import org.glycoinfo.ms.GlycanMassUtility.om.IMassElement;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;
import ninja.mspp.proformabuilder.Model.ChargeCarrier;
import ninja.mspp.proformabuilder.Model.ChargeState;
import ninja.mspp.proformabuilder.Model.FixedModification;
import ninja.mspp.proformabuilder.Model.GlobalIsotope;
import ninja.mspp.proformabuilder.Model.Modification;
import ninja.mspp.proformabuilder.Model.Peptidoform;
import ninja.mspp.proformabuilder.Model.PeptidoformIon;
import ninja.mspp.proformabuilder.Model.Range;
import ninja.mspp.proformabuilder.Model.Residue;
import ninja.mspp.proformabuilder.Model.Unlocalized;
import ninja.mspp.proformabuilder.cv.ModificationResolver;

/**
 * A parsed ProForma 2.1 peptidoform ion with mass calculation.
 * <pre>
 * ProForma p = ProForma.parse("PEPS[Phospho]TC[Carbamidomethyl]M[Oxidation]K/2");
 * p.getMonoisotopicMass();   // 1044.365734
 * p.getAverageMass();        // 1045.0838
 * p.getMz();                 // 523.190143
 * p.getFragments(FragmentIonType.B, 1);
 * </pre>
 * <p>
 * ProForma implements {@link IMassElement}: the masses are those of the neutral molecule
 * (local charges of charged formulas included), so it can be combined with {@code IonCloud} and
 * {@code MassUtils} like the other mass elements.
 * </p>
 * <ul>
 *   <li>Modification masses come from compositions (Unimod, PSI-MOD, formula, glycan); a tag with
 *       alternatives {@code [A|B]} uses the first descriptor that has a mass.</li>
 *   <li>Negative charges {@code /-2} in the string add electrons (ProForma 2.1 section 11.5), while a negative
 *       charge given to {@link #getMz(int)} removes protons.</li>
 *   <li>B and Z have two possible compositions: {@link #getMonoisotopicMass()} raises an exception and
 *       {@link #getMonoisotopicMassCandidates()} lists all possible masses.</li>
 *   <li>Fragments exclude labile modifications by default.</li>
 * </ul>
 */
public final class ProForma implements IMassElement {

	private static final int MZ_SCALE = 10;

	private final PeptidoformIon ion;
	private final ModificationResolver resolver;

	private ProForma(PeptidoformIon ion, ModificationResolver resolver) {
		this.ion = ion;
		this.resolver = resolver;
	}

	/** Parses one peptidoform ion with the default resolver (bundled Unimod and PSI-MOD). */
	public static ProForma parse(String s) {
		return new ProForma(ProFormaParser.parse(s), null);
	}

	/**
	 * As {@link #parse(String)}, but returns {@code null} for blank input
	 * (the convention of the parsers in GlycanMassUtility).
	 */
	public static ProForma parseOrNull(String s) {
		if ( s == null || s.trim().isEmpty() )
			return null;
		return parse(s.trim());
	}

	/** Parses one peptidoform ion with the given resolver. */
	public static ProForma parse(String s, ModificationResolver resolver) {
		return new ProForma(ProFormaParser.parse(s), resolver);
	}

	/** Parses peptidoform ions separated by {@code +} (chimeric spectra, ProForma 2.1 section 11.4). */
	public static List<ProForma> parseChimeric(String s) {
		List<ProForma> out = new ArrayList<>();
		for ( PeptidoformIon i : ProFormaParser.parseChimeric(s) )
			out.add(new ProForma(i, null));
		return out;
	}

	/** Wraps an existing model. */
	public static ProForma of(PeptidoformIon ion) {
		return new ProForma(ion, null);
	}

	private ModificationResolver resolver() {
		return this.resolver != null ? this.resolver : ModificationResolver.getDefault();
	}

	/** The parsed structure. */
	public PeptidoformIon getModel() {
		return this.ion;
	}

	/** Plain sequence; chains of a cross-linked ion are joined with {@code //}. */
	public String getSequence() {
		StringBuilder sb = new StringBuilder();
		for ( Peptidoform p : this.ion.getChains() ) {
			if ( sb.length() > 0 )
				sb.append("//");
			sb.append(p.getSequence());
		}
		return sb.toString();
	}

	/** Charge state as written, or {@code null}. */
	public ChargeState getChargeState() {
		return this.ion.getChargeState();
	}

	// ------------------------------------------------------------------ IMassElement

	/** The ProForma string. */
	@Override
	public String getName() {
		return toString();
	}

	/** Monoisotopic mass (same as {@link #getMonoisotopicMass()}). */
	@Override
	public BigDecimal getMass() {
		return getMonoisotopicMass();
	}

	/** Monoisotopic neutral mass. Raises {@link ProFormaException} for B and Z. */
	@Override
	public BigDecimal getMonoisotopicMass() {
		return getMassValue().getMonoisotopicMass();
	}

	/** Average neutral mass; approximate if {@link #isAverageMassApproximate()}. */
	@Override
	public BigDecimal getAverageMass() {
		return getMassValue().getAverageMass();
	}

	// ------------------------------------------------------------------ mass

	/** Total mass contribution (all chains, global isotopes applied). */
	public MassValue getMassValue() {
		return totalMass(null);
	}

	/** Elemental composition of the whole ion. Raises if a part is known only as a mass. */
	public Molecule getComposition() {
		MassValue m = getMassValue();
		if ( !m.hasExactComposition() )
			throw new ProFormaException("The composition is not exact: some modifications are given only as a mass");
		return m.getMolecule();
	}

	/** {@code true} if a modification is known only as a mass, so the average mass is approximate. */
	public boolean isAverageMassApproximate() {
		return getMassValue().isAverageApproximate();
	}

	/** All possible monoisotopic masses when B (D/N) or Z (E/Q) is present; one value otherwise. */
	public List<BigDecimal> getMonoisotopicMassCandidates() {
		List<int[]> ambiguous = new ArrayList<>();
		int offset = 0;
		for ( Peptidoform p : this.ion.getChains() ) {
			List<Residue> rs = p.getResidues();
			for ( int r = 0; r < rs.size(); r++ )
				if ( AminoAcids.isMassAmbiguous(rs.get(r).getAminoAcid()) )
					ambiguous.add(new int[] {offset + r});
			offset += rs.size();
		}
		int n = ambiguous.size();
		if ( n > 20 )
			throw new ProFormaException("Too many ambiguous residues (" + n + ")");
		TreeSet<BigDecimal> masses = new TreeSet<>();
		for ( int bits = 0; bits < (1 << n); bits++ ) {
			int[] choice = new int[offset];
			for ( int k = 0; k < n; k++ )
				choice[ambiguous.get(k)[0]] = (bits >> k) & 1;
			masses.add(totalMass(choice).getMonoisotopicMass().stripTrailingZeros());
		}
		return new ArrayList<>(masses);
	}

	private MassValue totalMass(int[] choice) {
		MassValue total = MassValue.ZERO;
		int offset = 0;
		for ( Peptidoform p : this.ion.getChains() ) {
			total = total.plus(chainMass(p, choice, offset));
			offset += p.getResidues().size();
		}
		return applyGlobalIsotopes(total);
	}

	private MassValue applyGlobalIsotopes(MassValue m) {
		for ( GlobalIsotope g : this.ion.getGlobalIsotopes() )
			m = m.withGlobalIsotope(g.getElement(), g.getIsotope());
		return m;
	}

	private MassValue chainMass(Peptidoform p, int[] choice, int offset) {
		MassValue m = MassValue.of(Compositions.H2O);
		List<Residue> rs = p.getResidues();
		for ( int r = 0; r < rs.size(); r++ ) {
			Residue res = rs.get(r);
			m = m.plus(MassValue.of(residueComposition(res.getAminoAcid(), choice, offset + r)));
			for ( Modification mod : res.getModifications() )
				m = m.plus(modificationMass(mod, res.getAminoAcid()));
			m = m.plus(fixedMass(p, r));
		}
		for ( Modification mod : p.getNTermModifications() )
			m = m.plus(modificationMass(mod, null));
		for ( Modification mod : p.getCTermModifications() )
			m = m.plus(modificationMass(mod, null));
		for ( Modification mod : p.getLabileModifications() )
			m = m.plus(modificationMass(mod, null));
		for ( Unlocalized u : p.getUnlocalizedModifications() )
			m = m.plus(modificationMass(u.getModification(), null).times(u.getCount()));
		for ( Range range : p.getRanges() )
			for ( Modification mod : range.getModifications() )
				m = m.plus(modificationMass(mod, null));
		return m;
	}

	private static Molecule residueComposition(char aa, int[] choice, int index) {
		if ( AminoAcids.isMassAmbiguous(aa) ) {
			if ( choice == null )
				throw new ProFormaException("Residue " + aa + " has two possible masses; use getMonoisotopicMassCandidates()");
			return AminoAcids.candidates(aa).get(choice[index]);
		}
		return AminoAcids.residue(aa);
	}

	/** Fixed modifications applying to residue {@code r} of the chain (including terminal targets). */
	private MassValue fixedMass(Peptidoform p, int r) {
		MassValue m = MassValue.ZERO;
		char aa = p.getResidues().get(r).getAminoAcid();
		boolean first = r == 0;
		boolean last = r == p.getResidues().size() - 1;
		for ( FixedModification f : this.ion.getFixedModifications() ) {
			for ( String t : f.getTargets() ) {
				boolean applies;
				if ( t.length() == 1 )
					applies = t.charAt(0) == aa;
				else if ( t.equals("N-term") )
					applies = first;
				else if ( t.equals("C-term") )
					applies = last;
				else if ( t.startsWith("N-term:") )
					applies = first && t.charAt(7) == aa;
				else
					applies = t.startsWith("C-term:") && last && t.charAt(7) == aa;
				if ( applies )
					m = m.plus(modificationMass(f.getModification(), aa));
			}
		}
		return m;
	}

	/** Mass of one tag: the first descriptor that has a mass; group references have none. */
	MassValue modificationMass(Modification mod, Character residue) {
		for ( Descriptor d : mod.getDescriptors() )
			if ( ModificationResolver.hasMass(d) )
				return resolver().resolve(d, residue);
		return MassValue.ZERO;
	}

	// ------------------------------------------------------------------ charge and m/z

	/** Total charge: charge state plus local charges of charged formulas; {@code null} if neither is given. */
	public Integer getCharge() {
		int local = getMassValue().getCharge();
		if ( this.ion.getChargeState() == null )
			return local != 0 ? Integer.valueOf(local) : null;
		return this.ion.getChargeState().getTotalCharge() + local;
	}

	/** m/z using the charge state of the string (and local charges). */
	public BigDecimal getMz() {
		MassValue m = getMassValue();
		ChargeState cs = this.ion.getChargeState();
		int z = m.getCharge();
		BigDecimal carrierMass = BigDecimal.ZERO;
		if ( cs != null ) {
			if ( cs.hasCarriers() ) {
				for ( ChargeCarrier c : cs.getCarriers() ) {
					MassValue carrier = applyGlobalIsotopes(MassValue.of(c.getMolecule()));
					carrierMass = carrierMass.add(carrier.getMonoisotopicMass().multiply(new BigDecimal(c.getCount())));
				}
			} else if ( cs.getCharge() > 0 ) {
				carrierMass = IonType.H.getMonoisotopicMass().multiply(new BigDecimal(cs.getCharge()));
			} else {
				// negative charge: electrons added (ProForma 2.1 section 11.5)
				carrierMass = AtomType.e.getMonoisotopicMass().multiply(new BigDecimal(-cs.getCharge()));
			}
			z += cs.getTotalCharge();
		}
		if ( z == 0 )
			throw new ProFormaException("No charge: give a charge state (/2) or use getMz(int)");
		return m.getMonoisotopicMass().add(carrierMass).divide(new BigDecimal(Math.abs(z)), MZ_SCALE, RoundingMode.HALF_UP);
	}

	/**
	 * m/z of the protonated or deprotonated ion, ignoring the charge state of the string.
	 * <p>
	 * A positive {@code z} adds {@code z} protons, i.e. [M+zH]<sup>z+</sup>, and a negative {@code z}
	 * removes {@code -z} protons, i.e. [M-zH]<sup>z-</sup>, as in negative ion mode mass spectrometry.
	 * Note that a negative charge state in the string (e.g. {@code /-2}) adds electrons instead,
	 * following ProForma 2.1 section 11.5, so {@link #getMz()} of {@code PEPTIDE/-2} differs from
	 * {@code getMz(-2)} of {@code PEPTIDE}.
	 * </p>
	 * @param z charge, not 0
	 * @return m/z
	 * @throws IllegalArgumentException if {@code z} is 0
	 */
	public BigDecimal getMz(int z) {
		if ( z == 0 )
			throw new IllegalArgumentException("Charge must not be 0");
		MassValue m = getMassValue();
		BigDecimal mass = m.getMonoisotopicMass().add(IonType.H.getMonoisotopicMass().multiply(new BigDecimal(z)));
		return mass.divide(new BigDecimal(Math.abs(z + m.getCharge())), MZ_SCALE, RoundingMode.HALF_UP);
	}

	// ------------------------------------------------------------------ fragments

	/** Fragment m/z values (labile modifications excluded, unlocalized ones included). */
	public List<BigDecimal> getFragments(FragmentIonType type, int charge) {
		return getFragments(type, charge, false, true);
	}

	/**
	 * Fragment ions of one series, in order of increasing length (b1, b2, ... / y1, y2, ...).
	 * @param charge             charge of the fragments; 0 returns neutral masses
	 * @param includeLabile      include labile modifications {@code {...}} in every fragment
	 * @param includeUnlocalized include modifications of unknown position {@code [...]?} in every fragment
	 */
	public List<BigDecimal> getFragments(FragmentIonType type, int charge,
			boolean includeLabile, boolean includeUnlocalized) {
		if ( this.ion.getChains().size() != 1 )
			throw new ProFormaException("Fragments are supported for a single chain only");
		Peptidoform p = this.ion.getChains().get(0);
		List<Residue> rs = p.getResidues();
		int n = rs.size();
		MassValue m = MassValue.of(Compositions.plus(type.getShift(), Compositions.H2O));
		for ( Modification mod : type.isCTerminal() ? p.getCTermModifications() : p.getNTermModifications() )
			m = m.plus(modificationMass(mod, null));
		if ( includeLabile )
			for ( Modification mod : p.getLabileModifications() )
				m = m.plus(modificationMass(mod, null));
		if ( includeUnlocalized )
			for ( Unlocalized u : p.getUnlocalizedModifications() )
				m = m.plus(modificationMass(u.getModification(), null).times(u.getCount()));

		BigDecimal proton = IonType.H.getMonoisotopicMass();
		List<Range> pending = new ArrayList<>(p.getRanges());
		List<BigDecimal> out = new ArrayList<>();
		for ( int k = 0; k < n - 1; k++ ) {
			int r = type.isCTerminal() ? n - 1 - k : k;
			Residue res = rs.get(r);
			m = m.plus(MassValue.of(residueComposition(res.getAminoAcid(), null, r)));
			for ( Modification mod : res.getModifications() )
				m = m.plus(modificationMass(mod, res.getAminoAcid()));
			m = m.plus(fixedMass(p, r));
			// a range modification is added once the fragment contains any residue of the range
			for ( Iterator<Range> it = pending.iterator(); it.hasNext(); ) {
				Range range = it.next();
				if ( r >= range.getStart() && r < range.getEnd() ) {
					for ( Modification mod : range.getModifications() )
						m = m.plus(modificationMass(mod, null));
					it.remove();
				}
			}
			BigDecimal neutral = applyGlobalIsotopes(m).getMonoisotopicMass();
			out.add(charge == 0 ? neutral : neutral.add(proton.multiply(new BigDecimal(charge)))
					.divide(new BigDecimal(Math.abs(charge)), MZ_SCALE, RoundingMode.HALF_UP));
		}
		return out;
	}

	// ------------------------------------------------------------------ output

	/** ProForma string (see {@link ProFormaWriter}). */
	@Override
	public String toString() {
		return ProFormaWriter.write(this.ion);
	}

	/**
	 * Two ProForma objects are equal if their parsed structures are equal and they use the same resolver
	 * (the default resolver when none is given).
	 * Strings that differ only in notation, e.g. in the case of residues, give equal objects.
	 */
	@Override
	public boolean equals(Object o) {
		if ( this == o )
			return true;
		if ( !(o instanceof ProForma) )
			return false;
		ProForma other = (ProForma) o;
		return this.ion.equals(other.ion) && this.resolver == other.resolver;
	}

	@Override
	public int hashCode() {
		return this.ion.hashCode();
	}
}
