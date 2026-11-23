package ninja.mspp.proformabuilder;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.glycoinfo.ms.GlycanMassUtility.dict.residue.aminoacid.AminoAcidDictionary;
import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;

/**
 * Amino acid residues for ProForma (free amino acid of {@link AminoAcidDictionary} minus H2O).
 * <p>Ambiguous residues of ProForma 2.1 section 7.1:</p>
 * <ul>
 *   <li>B = D or N, Z = E or Q: two possible compositions ({@link #candidates(char)}).</li>
 *   <li>J = I or L: same composition.</li>
 *   <li>X = unknown residue: empty composition; its mass is given by an attached modification.</li>
 * </ul>
 */
public final class AminoAcids {

	private AminoAcids() { }

	/** {@code true} if the upper-case letter is a ProForma amino acid symbol. */
	public static boolean isValid(char aa) {
		return AminoAcidDictionary.findAminoAcid(String.valueOf(aa)) || aa == 'B' || aa == 'Z' || aa == 'J' || aa == 'X';
	}

	/** {@code true} for residues with two possible compositions (B, Z). */
	public static boolean isMassAmbiguous(char aa) {
		return aa == 'B' || aa == 'Z';
	}

	/** Residue composition. Throws for B and Z. */
	public static Molecule residue(char aa) {
		switch ( aa ) {
		case 'J':
			return standard('L');
		case 'X':
			return Compositions.EMPTY;
		case 'B':
		case 'Z':
			throw new IllegalStateException("Residue " + aa + " has two possible compositions");
		default:
			return standard(aa);
		}
	}

	/** All possible compositions (two for B and Z, one otherwise). */
	public static List<Molecule> candidates(char aa) {
		if ( aa == 'B' )
			return Arrays.asList(standard('D'), standard('N'));
		if ( aa == 'Z' )
			return Arrays.asList(standard('E'), standard('Q'));
		return Collections.singletonList(residue(aa));
	}

	private static Molecule standard(char aa) {
		String s = String.valueOf(aa);
		if ( !AminoAcidDictionary.findAminoAcid(s) )
			throw new IllegalArgumentException("Unknown amino acid: " + aa);
		return AminoAcidDictionary.getAminoAcid(s).getMolecule().addMolecules(Compositions.H2O, -1);
	}
}
