package ninja.mspp.proformabuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.glycoinfo.ms.GlycanMassUtility.dict.atom.AtomType;
import org.glycoinfo.ms.GlycanMassUtility.dict.atom.IsotopeType;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Parsed ProForma structure. The layers follow ProForma 2.1 section 5: a peptidoform ion consists
 * of one or more peptidoforms (chains joined by {@code //}), a charge state and global modifications.
 * All classes are immutable; equality compares all fields.
 */
public final class Model {

	private Model() { }

	static <T> List<T> copy(List<T> l) {
		return Collections.unmodifiableList(new ArrayList<T>(l));
	}

	/** Base class giving value-based equality from {@link #fields()}. */
	abstract static class Value {
		abstract Object[] fields();

		@Override
		public boolean equals(Object o) {
			return o != null && o.getClass() == getClass() && Arrays.equals(fields(), ((Value) o).fields());
		}

		@Override
		public int hashCode() {
			return Arrays.hashCode(fields());
		}

		@Override
		public String toString() {
			return getClass().getSimpleName() + Arrays.toString(fields());
		}
	}

	/** Group reference: {@code #g1}, {@code #g1(0.90)}, {@code #XL1}, {@code #BRANCH}. */
	public static final class Group extends Value {
		private final String label;
		private final Double score;
		private final String scoreText;

		Group(String label, Double score, String scoreText) {
			this.label = label;
			this.score = score;
			this.scoreText = scoreText;
		}

		/** Group name without '#'. */
		public String getLabel() {
			return this.label;
		}

		/** Localisation score, or {@code null}. */
		public Double getScore() {
			return this.score;
		}

		/** Score as written (e.g. {@code 0.90}), or {@code null}. */
		public String getScoreText() {
			return this.scoreText;
		}

		public boolean isCrossLink() {
			return this.label.startsWith("XL");
		}

		public boolean isBranch() {
			return this.label.equals("BRANCH");
		}

		@Override
		Object[] fields() {
			return new Object[] {this.label, this.score, this.scoreText};
		}
	}

	/** A modification tag: descriptors separated by {@code |} plus an optional group. */
	public static final class Modification extends Value {
		private final List<Descriptor> descriptors;
		private final Group group;

		Modification(List<Descriptor> descriptors, Group group) {
			this.descriptors = copy(descriptors);
			this.group = group;
		}

		public List<Descriptor> getDescriptors() {
			return this.descriptors;
		}

		/** Group, or {@code null}. */
		public Group getGroup() {
			return this.group;
		}

		/** {@code true} if this tag only refers to a group ({@code [#g1]}) and carries no mass itself. */
		public boolean isGroupReferenceOnly() {
			return this.descriptors.isEmpty();
		}

		@Override
		Object[] fields() {
			return new Object[] {this.descriptors, this.group};
		}
	}

	/** One residue (upper-case letter) with its modifications. */
	public static final class Residue extends Value {
		private final char aa;
		private final List<Modification> modifications;

		Residue(char aa, List<Modification> modifications) {
			this.aa = aa;
			this.modifications = copy(modifications);
		}

		public char getAminoAcid() {
			return this.aa;
		}

		public List<Modification> getModifications() {
			return this.modifications;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.aa, this.modifications};
		}
	}

	/** A modification of unknown position: {@code [Phospho]^2?}. */
	public static final class Unlocalized extends Value {
		private final Modification modification;
		private final int count;

		Unlocalized(Modification modification, int count) {
			this.modification = modification;
			this.count = count;
		}

		public Modification getModification() {
			return this.modification;
		}

		public int getCount() {
			return this.count;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.modification, this.count};
		}
	}

	/** A range of residues with modifications: {@code (ESFRMS)[+19.0523]}; end is exclusive. */
	public static final class Range extends Value {
		private final int start;
		private final int end;
		private final List<Modification> modifications;

		Range(int start, int end, List<Modification> modifications) {
			this.start = start;
			this.end = end;
			this.modifications = copy(modifications);
		}

		public int getStart() {
			return this.start;
		}

		public int getEnd() {
			return this.end;
		}

		public List<Modification> getModifications() {
			return this.modifications;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.start, this.end, this.modifications};
		}
	}

	/** Sequence ambiguity {@code (?DQ)}; end is exclusive. */
	public static final class SequenceAmbiguity extends Value {
		private final int start;
		private final int end;

		SequenceAmbiguity(int start, int end) {
			this.start = start;
			this.end = end;
		}

		public int getStart() {
			return this.start;
		}

		public int getEnd() {
			return this.end;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.start, this.end};
		}
	}

	/** One peptidoform (a single chain). */
	public static final class Peptidoform extends Value {
		private final String name;
		private final List<Modification> labile;
		private final List<Unlocalized> unlocalized;
		private final List<Modification> nTerm;
		private final List<Residue> residues;
		private final List<Modification> cTerm;
		private final List<Range> ranges;
		private final List<SequenceAmbiguity> ambiguities;

		Peptidoform(String name, List<Modification> labile, List<Unlocalized> unlocalized, List<Modification> nTerm,
				List<Residue> residues, List<Modification> cTerm, List<Range> ranges, List<SequenceAmbiguity> ambiguities) {
			this.name = name;
			this.labile = copy(labile);
			this.unlocalized = copy(unlocalized);
			this.nTerm = copy(nTerm);
			this.residues = copy(residues);
			this.cTerm = copy(cTerm);
			this.ranges = copy(ranges);
			this.ambiguities = copy(ambiguities);
		}

		/** Peptidoform name {@code (>name)}, or {@code null}. */
		public String getName() {
			return this.name;
		}

		public List<Modification> getLabileModifications() {
			return this.labile;
		}

		public List<Unlocalized> getUnlocalizedModifications() {
			return this.unlocalized;
		}

		public List<Modification> getNTermModifications() {
			return this.nTerm;
		}

		public List<Residue> getResidues() {
			return this.residues;
		}

		public List<Modification> getCTermModifications() {
			return this.cTerm;
		}

		public List<Range> getRanges() {
			return this.ranges;
		}

		public List<SequenceAmbiguity> getSequenceAmbiguities() {
			return this.ambiguities;
		}

		/** Plain sequence (upper case). */
		public String getSequence() {
			StringBuilder sb = new StringBuilder(this.residues.size());
			for ( Residue r : this.residues )
				sb.append(r.getAminoAcid());
			return sb.toString();
		}

		@Override
		Object[] fields() {
			return new Object[] {this.name, this.labile, this.unlocalized, this.nTerm, this.residues, this.cTerm,
					this.ranges, this.ambiguities};
		}
	}

	/** Global isotope replacement {@code <13C>}, {@code <15N>}, {@code <D>}. */
	public static final class GlobalIsotope extends Value {
		private final AtomType element;
		private final IsotopeType isotope;
		private final String text;

		GlobalIsotope(AtomType element, IsotopeType isotope, String text) {
			this.element = element;
			this.isotope = isotope;
			this.text = text;
		}

		public AtomType getElement() {
			return this.element;
		}

		public IsotopeType getIsotope() {
			return this.isotope;
		}

		/** Text between the angle brackets. */
		public String getText() {
			return this.text;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.element, this.isotope, this.text};
		}
	}

	/** Fixed modification {@code <[Carbamidomethyl]@C,N-term:Q>}. */
	public static final class FixedModification extends Value {
		private final Modification modification;
		private final List<String> targets;

		FixedModification(Modification modification, List<String> targets) {
			this.modification = modification;
			this.targets = copy(targets);
		}

		public Modification getModification() {
			return this.modification;
		}

		/** Residue letters, {@code N-term}, {@code C-term}, {@code N-term:Q}, {@code C-term:G}. */
		public List<String> getTargets() {
			return this.targets;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.modification, this.targets};
		}
	}

	/** Charge carrier {@code Na:z+1^2}: composition (with its charge) and number of carriers. */
	public static final class ChargeCarrier extends Value {
		private final Molecule molecule;
		private final int count;
		private final String text;

		ChargeCarrier(Molecule molecule, int count, String text) {
			this.molecule = molecule;
			this.count = count;
			this.text = text;
		}

		/** Carrier including its charge (electron masses accounted for). */
		public Molecule getMolecule() {
			return this.molecule;
		}

		public int getCount() {
			return this.count;
		}

		public String getText() {
			return this.text;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.molecule, this.count, this.text};
		}
	}

	/**
	 * Charge state: a plain charge ({@code /2}: protons for positive, electrons for negative charges,
	 * ProForma 2.1 section 11.5) or charge carriers ({@code /[Na:z+1,H:z+1]}).
	 */
	public static final class ChargeState extends Value {
		private final int charge;
		private final List<ChargeCarrier> carriers;

		ChargeState(int charge, List<ChargeCarrier> carriers) {
			this.charge = charge;
			this.carriers = copy(carriers);
		}

		/** Plain charge ({@code /2}); 0 when carriers are given. */
		public int getCharge() {
			return this.charge;
		}

		public List<ChargeCarrier> getCarriers() {
			return this.carriers;
		}

		public boolean hasCarriers() {
			return !this.carriers.isEmpty();
		}

		/** Total charge of the carriers, or the plain charge. */
		public int getTotalCharge() {
			if ( this.carriers.isEmpty() )
				return this.charge;
			int z = 0;
			for ( ChargeCarrier c : this.carriers )
				z += c.getMolecule().getCharges() * c.getCount();
			return z;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.charge, this.carriers};
		}
	}

	/** A peptidoform ion: chains joined by {@code //}, global modifications and charge. */
	public static final class PeptidoformIon extends Value {
		private final String name;
		private final List<GlobalIsotope> globalIsotopes;
		private final List<FixedModification> fixedModifications;
		private final List<Peptidoform> chains;
		private final ChargeState charge;

		PeptidoformIon(String name, List<GlobalIsotope> globalIsotopes, List<FixedModification> fixedModifications,
				List<Peptidoform> chains, ChargeState charge) {
			this.name = name;
			this.globalIsotopes = copy(globalIsotopes);
			this.fixedModifications = copy(fixedModifications);
			this.chains = copy(chains);
			this.charge = charge;
		}

		/** Ion name {@code (>>name)}, or {@code null}. */
		public String getName() {
			return this.name;
		}

		public List<GlobalIsotope> getGlobalIsotopes() {
			return this.globalIsotopes;
		}

		public List<FixedModification> getFixedModifications() {
			return this.fixedModifications;
		}

		public List<Peptidoform> getChains() {
			return this.chains;
		}

		/** Charge state, or {@code null} if not given. */
		public ChargeState getChargeState() {
			return this.charge;
		}

		@Override
		Object[] fields() {
			return new Object[] {this.name, this.globalIsotopes, this.fixedModifications, this.chains, this.charge};
		}
	}
}
