package ninja.mspp.proformabuilder;

import java.math.BigDecimal;

import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * One element of a modification tag. A tag {@code [A|B|C]} holds descriptors separated by
 * {@code |} (ProForma 2.1 sections 6.2 and 7.5). Each descriptor keeps its source text, so that
 * the writer reproduces the input; two descriptors are equal when their type and text are equal.
 */
public abstract class Descriptor {

	/** Controlled vocabularies (ProForma 2.1 section 6.2). */
	public enum Cv {
		UNIMOD("U", "UNIMOD"), PSI_MOD("M", "MOD"), RESID("R", "RESID"),
		XL_MOD("X", "XLMOD"), GNO("G", "GNO"), BRNO("B", "BRNO"), OBSERVED("Obs", null);

		private final String shortPrefix;
		private final String accessionPrefix;

		Cv(String shortPrefix, String accessionPrefix) {
			this.shortPrefix = shortPrefix;
			this.accessionPrefix = accessionPrefix;
		}

		public String getShortPrefix() {
			return this.shortPrefix;
		}

		public String getAccessionPrefix() {
			return this.accessionPrefix;
		}

		static Cv byShortPrefix(String p) {
			for ( Cv cv : values() )
				if ( cv.shortPrefix.equalsIgnoreCase(p) )
					return cv;
			return null;
		}

		static Cv byAccessionPrefix(String p) {
			for ( Cv cv : values() )
				if ( cv.accessionPrefix != null && cv.accessionPrefix.equalsIgnoreCase(p) )
					return cv;
			return null;
		}
	}

	private final String text;

	protected Descriptor(String text) {
		this.text = text;
	}

	/** Source text of this descriptor, as written. */
	public String getText() {
		return this.text;
	}

	@Override
	public boolean equals(Object o) {
		return o != null && o.getClass() == getClass() && ((Descriptor) o).text.equals(this.text);
	}

	@Override
	public int hashCode() {
		return getClass().hashCode() * 31 + this.text.hashCode();
	}

	@Override
	public String toString() {
		return this.text;
	}

	/** Modification name, optionally with a CV prefix: {@code Phospho}, {@code U:Phospho}. */
	public static final class Named extends Descriptor {
		private final Cv cv;
		private final String name;

		Named(Cv cv, String name, String text) {
			super(text);
			this.cv = cv;
			this.name = name;
		}

		/** CV prefix, or {@code null} if none was written. */
		public Cv getCv() {
			return this.cv;
		}

		public String getName() {
			return this.name;
		}
	}

	/** Accession: {@code UNIMOD:21}, {@code MOD:00046}, {@code RESID:AA0581}, {@code GNO:G59626AS}. */
	public static final class Accession extends Descriptor {
		private final Cv cv;
		private final String id;

		Accession(Cv cv, String id, String text) {
			super(text);
			this.cv = cv;
			this.id = id;
		}

		public Cv getCv() {
			return this.cv;
		}

		public String getId() {
			return this.id;
		}
	}

	/** Mass delta, optionally prefixed: {@code +15.995}, {@code U:+15.995}, {@code Obs:+79.978}. */
	public static final class MassDelta extends Descriptor {
		private final Cv cv;
		private final BigDecimal value;

		MassDelta(Cv cv, BigDecimal value, String text) {
			super(text);
			this.cv = cv;
			this.value = value;
		}

		public Cv getCv() {
			return this.cv;
		}

		public BigDecimal getValue() {
			return this.value;
		}
	}

	/** Molecular formula, optionally charged (Level 3): {@code Formula:HPO3}, {@code Formula:Zn1:z+2}. */
	public static final class FormulaTag extends Descriptor {
		private final Molecule molecule;

		FormulaTag(Molecule molecule, String text) {
			super(text);
			this.molecule = molecule;
		}

		/** Composition including the charge (electron masses accounted for). */
		public Molecule getMolecule() {
			return this.molecule;
		}
	}

	/** Glycan composition: {@code Glycan:Hex5HexNAc4NeuAc1}. */
	public static final class GlycanTag extends Descriptor {
		private final GlycanComposition glycan;

		GlycanTag(GlycanComposition glycan, String text) {
			super(text);
			this.glycan = glycan;
		}

		public GlycanComposition getGlycan() {
			return this.glycan;
		}
	}

	/** Free text: {@code INFO:any text}. Has no mass. */
	public static final class Info extends Descriptor {
		private final String value;

		Info(String value, String text) {
			super(text);
			this.value = value;
		}

		public String getValue() {
			return this.value;
		}
	}

	/** Placement control (Level 3): {@code Position:M}, {@code Limit:2}, {@code CoMKP}, {@code CoMUP}. No mass. */
	public static final class Placement extends Descriptor {
		private final String kind;
		private final String value;

		Placement(String kind, String value, String text) {
			super(text);
			this.kind = kind;
			this.value = value;
		}

		public String getKind() {
			return this.kind;
		}

		/** Value after the colon, or {@code null} for CoMKP/CoMUP. */
		public String getValue() {
			return this.value;
		}
	}
}
