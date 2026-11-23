package ninja.mspp.proformabuilder;

import java.math.BigDecimal;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.dict.atom.IsotopeType;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * A mass contribution: an elemental composition ({@link Molecule}, possibly charged) plus mass
 * known only as a number (ProForma mass deltas such as {@code [+79.966]}).
 * <p>
 * For mass-only parts the average mass is unknown; the given value is used for both masses and
 * {@link #isAverageApproximate()} becomes {@code true}.
 * </p>
 */
public final class MassValue {

	public static final MassValue ZERO = new MassValue(Compositions.EMPTY, BigDecimal.ZERO, false);

	private final Molecule molecule;
	private final BigDecimal extraMass;
	private final boolean averageApproximate;

	public MassValue(Molecule molecule, BigDecimal extraMass, boolean averageApproximate) {
		this.molecule = molecule;
		this.extraMass = extraMass;
		this.averageApproximate = averageApproximate;
	}

	public static MassValue of(Molecule m) {
		return new MassValue(m, BigDecimal.ZERO, false);
	}

	public static MassValue ofDelta(BigDecimal delta) {
		return new MassValue(Compositions.EMPTY, delta, true);
	}

	public MassValue plus(MassValue o) {
		return new MassValue(Compositions.plus(this.molecule, o.molecule), this.extraMass.add(o.extraMass),
				this.averageApproximate || o.averageApproximate);
	}

	public MassValue times(int n) {
		return new MassValue(Compositions.times(this.molecule, n), this.extraMass.multiply(new BigDecimal(n)),
				this.averageApproximate && n != 0);
	}

	public MassValue withGlobalIsotope(AtomType element, IsotopeType isotope) {
		return new MassValue(Compositions.withGlobalIsotope(this.molecule, element, isotope), this.extraMass,
				this.averageApproximate);
	}

	/** Elemental part (including its charge). */
	public Molecule getMolecule() {
		return this.molecule;
	}

	/** Mass-only part. */
	public BigDecimal getExtraMass() {
		return this.extraMass;
	}

	public boolean isAverageApproximate() {
		return this.averageApproximate;
	}

	/** Net charge of charged formulas (ProForma Level 3). */
	public int getCharge() {
		return this.molecule.getCharges();
	}

	/** Monoisotopic mass (electron masses of charges are already accounted for by {@link Molecule}). */
	public BigDecimal getMonoisotopicMass() {
		return this.molecule.getMonoisotopicMass().add(this.extraMass);
	}

	public BigDecimal getAverageMass() {
		return this.molecule.getAverageMass().add(this.extraMass);
	}

	/** {@code true} if the value is an elemental composition only. */
	public boolean hasExactComposition() {
		return !this.averageApproximate && this.extraMass.signum() == 0;
	}
}
