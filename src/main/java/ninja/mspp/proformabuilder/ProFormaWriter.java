package ninja.mspp.proformabuilder;

import java.util.ArrayList;
import java.util.List;

import ninja.mspp.proformabuilder.Model.ChargeCarrier;
import ninja.mspp.proformabuilder.Model.ChargeState;
import ninja.mspp.proformabuilder.Model.FixedModification;
import ninja.mspp.proformabuilder.Model.GlobalIsotope;
import ninja.mspp.proformabuilder.Model.Modification;
import ninja.mspp.proformabuilder.Model.Peptidoform;
import ninja.mspp.proformabuilder.Model.PeptidoformIon;
import ninja.mspp.proformabuilder.Model.Range;
import ninja.mspp.proformabuilder.Model.Residue;
import ninja.mspp.proformabuilder.Model.SequenceAmbiguity;
import ninja.mspp.proformabuilder.Model.Unlocalized;

/**
 * Writes the model back to a ProForma string. Descriptors are written as they were read and
 * residues in upper case, so that {@code parse(write(x))} reproduces {@code x}.
 */
public final class ProFormaWriter {

	private ProFormaWriter() { }

	public static String write(PeptidoformIon ion) {
		StringBuilder sb = new StringBuilder();
		for ( GlobalIsotope g : ion.getGlobalIsotopes() )
			sb.append('<').append(g.getText()).append('>');
		for ( FixedModification f : ion.getFixedModifications() ) {
			sb.append('<').append(tag(f.getModification(), '[', ']')).append('@');
			sb.append(join(f.getTargets(), ",")).append('>');
		}
		if ( ion.getName() != null )
			sb.append("(>>").append(ion.getName()).append(')');
		List<String> chains = new ArrayList<>();
		for ( Peptidoform p : ion.getChains() )
			chains.add(peptidoform(p));
		sb.append(join(chains, "//"));
		if ( ion.getChargeState() != null )
			sb.append('/').append(charge(ion.getChargeState()));
		return sb.toString();
	}

	public static String writeChimeric(List<PeptidoformIon> ions) {
		List<String> parts = new ArrayList<>();
		for ( PeptidoformIon i : ions )
			parts.add(write(i));
		return join(parts, "+");
	}

	private static String charge(ChargeState cs) {
		if ( !cs.hasCarriers() )
			return Integer.toString(cs.getCharge());
		List<String> parts = new ArrayList<>();
		for ( ChargeCarrier c : cs.getCarriers() )
			parts.add(c.getText());
		return "[" + join(parts, ",") + "]";
	}

	private static String peptidoform(Peptidoform p) {
		StringBuilder sb = new StringBuilder();
		if ( p.getName() != null )
			sb.append("(>").append(p.getName()).append(')');
		for ( Modification m : p.getLabileModifications() )
			sb.append(tag(m, '{', '}'));
		if ( !p.getUnlocalizedModifications().isEmpty() ) {
			for ( Unlocalized u : p.getUnlocalizedModifications() ) {
				sb.append(tag(u.getModification(), '[', ']'));
				if ( u.getCount() != 1 )
					sb.append('^').append(u.getCount());
			}
			sb.append('?');
		}
		if ( !p.getNTermModifications().isEmpty() ) {
			for ( Modification m : p.getNTermModifications() )
				sb.append(tag(m, '[', ']'));
			sb.append('-');
		}
		List<Residue> rs = p.getResidues();
		for ( int r = 0; r < rs.size(); r++ ) {
			for ( Range range : p.getRanges() )
				if ( range.getStart() == r )
					sb.append('(');
			for ( SequenceAmbiguity a : p.getSequenceAmbiguities() )
				if ( a.getStart() == r )
					sb.append("(?");
			sb.append(rs.get(r).getAminoAcid());
			for ( Modification m : rs.get(r).getModifications() )
				sb.append(tag(m, '[', ']'));
			for ( SequenceAmbiguity a : p.getSequenceAmbiguities() )
				if ( a.getEnd() == r + 1 )
					sb.append(')');
			for ( Range range : p.getRanges() )
				if ( range.getEnd() == r + 1 ) {
					sb.append(')');
					for ( Modification m : range.getModifications() )
						sb.append(tag(m, '[', ']'));
				}
		}
		if ( !p.getCTermModifications().isEmpty() ) {
			sb.append('-');
			for ( Modification m : p.getCTermModifications() )
				sb.append(tag(m, '[', ']'));
		}
		return sb.toString();
	}

	/** One tag, e.g. {@code [Phospho#g1(0.9)]}; the group is written after the first descriptor. */
	static String tag(Modification m, char open, char close) {
		StringBuilder sb = new StringBuilder().append(open);
		List<Descriptor> ds = m.getDescriptors();
		if ( ds.isEmpty() ) {
			sb.append(group(m));
		} else {
			for ( int k = 0; k < ds.size(); k++ ) {
				if ( k > 0 )
					sb.append('|');
				sb.append(ds.get(k).getText());
				if ( k == 0 )
					sb.append(group(m));
			}
		}
		return sb.append(close).toString();
	}

	private static String group(Modification m) {
		if ( m.getGroup() == null )
			return "";
		String s = "#" + m.getGroup().getLabel();
		if ( m.getGroup().getScoreText() != null )
			s += "(" + m.getGroup().getScoreText() + ")";
		return s;
	}

	private static String join(List<String> parts, String sep) {
		StringBuilder sb = new StringBuilder();
		for ( int k = 0; k < parts.size(); k++ ) {
			if ( k > 0 )
				sb.append(sep);
			sb.append(parts.get(k));
		}
		return sb.toString();
	}
}
