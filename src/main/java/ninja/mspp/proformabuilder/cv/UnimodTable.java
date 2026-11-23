package ninja.mspp.proformabuilder.cv;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.glycoinfo.ms.GlycanMassUtility.om.Molecule;
import ninja.mspp.proformabuilder.ProFormaFormula;

/**
 * Unimod modifications: record id, title, full name and delta composition.
 * <p>
 * Read from {@code unimod.xml} (the unmodified file bundled with the library, or another file
 * given to {@link ModificationResolver}).
 * Names are looked up with the capitalisation of the CV first; a case-insensitive match is used
 * only when it is unique (ProForma 2.1 section 5).
 * </p>
 */
public final class UnimodTable {

	/** One Unimod entry. */
	public static final class Entry {
		private final int id;
		private final String title;
		private final String fullName;
		private final Molecule composition;

		Entry(int id, String title, String fullName, Molecule composition) {
			this.id = id;
			this.title = title;
			this.fullName = fullName;
			this.composition = composition;
		}

		public int getId() {
			return this.id;
		}

		public String getTitle() {
			return this.title;
		}

		public String getFullName() {
			return this.fullName;
		}

		/** Delta composition. */
		public Molecule getComposition() {
			return this.composition;
		}
	}

	private final Map<Integer, Entry> byId = new LinkedHashMap<>();
	private final Map<String, Entry> byTitle = new HashMap<>();
	private final Map<String, Entry> byFullName = new HashMap<>();
	private final Map<String, List<Entry>> byLowerName = new HashMap<>();
	private UnimodTable() { }

	/** Reads {@code unimod.xml} (http://www.unimod.org/xml/unimod.xml). */
	public static UnimodTable fromXml(InputStream in) throws IOException {
		UnimodTable t = new UnimodTable();
		try {
			XMLInputFactory f = XMLInputFactory.newFactory();
			f.setProperty(XMLInputFactory.SUPPORT_DTD, false);
			f.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, false);
			XMLStreamReader r = f.createXMLStreamReader(in);
			String title = null;
			String fullName = null;
			int id = -1;
			boolean inMod = false;
			boolean inDelta = false;
			Molecule comp = null;
			while ( r.hasNext() ) {
				int ev = r.next();
				if ( ev == XMLStreamConstants.START_ELEMENT ) {
					String local = r.getLocalName();
					if ( local.equals("mod") ) {
						inMod = true;
						title = r.getAttributeValue(null, "title");
						fullName = r.getAttributeValue(null, "full_name");
						String rid = r.getAttributeValue(null, "record_id");
						id = rid == null ? -1 : Integer.parseInt(rid);
						comp = new Molecule();
					} else if ( inMod && local.equals("delta") ) {
						inDelta = true;
					} else if ( inDelta && local.equals("element") ) {
						String symbol = r.getAttributeValue(null, "symbol");
						String num = r.getAttributeValue(null, "number");
						int count = num == null ? 1 : Integer.parseInt(num);
						comp = comp.addMolecules(unimodAtom(symbol), count);
					}
				} else if ( ev == XMLStreamConstants.END_ELEMENT ) {
					String local = r.getLocalName();
					if ( local.equals("delta") ) {
						inDelta = false;
					} else if ( local.equals("mod") && inMod ) {
						inMod = false;
						t.add(new Entry(id, title, fullName, comp));
					}
				}
			}
			r.close();
		} catch ( XMLStreamException e ) {
			throw new IOException("Invalid unimod.xml: " + e.getMessage(), e);
		}
		return t;
	}

	/** Unimod element symbols: {@code C}, {@code 13C}, {@code 2H}, {@code 18O}. */
	private static Molecule unimodAtom(String symbol) {
		int i = 0;
		while ( i < symbol.length() && Character.isDigit(symbol.charAt(i)) )
			i++;
		int massNumber = i == 0 ? 0 : Integer.parseInt(symbol.substring(0, i));
		String el = symbol.substring(i);
		return ProFormaFormula.parse(massNumber == 0 ? el : "[" + massNumber + el + "]");
	}

	private void add(Entry e) {
		this.byId.put(e.id, e);
		if ( !this.byTitle.containsKey(e.title) )
			this.byTitle.put(e.title, e);
		if ( e.fullName != null && !this.byFullName.containsKey(e.fullName) )
			this.byFullName.put(e.fullName, e);
		addLower(e.title, e);
		if ( e.fullName != null && !e.fullName.equalsIgnoreCase(e.title) )
			addLower(e.fullName, e);
	}

	private void addLower(String name, Entry e) {
		String key = name.toLowerCase(Locale.ROOT);
		List<Entry> l = this.byLowerName.get(key);
		if ( l == null ) {
			l = new ArrayList<>();
			this.byLowerName.put(key, l);
		}
		l.add(e);
	}

	/** Looks up by title, then full name, then a unique case-insensitive match; {@code null} if not found. */
	public Entry byName(String name) {
		Entry e = this.byTitle.get(name);
		if ( e == null )
			e = this.byFullName.get(name);
		if ( e == null ) {
			List<Entry> l = this.byLowerName.get(name.toLowerCase(Locale.ROOT));
			if ( l != null ) {
				Entry first = l.get(0);
				boolean unique = true;
				for ( Entry x : l )
					unique &= x.id == first.id;
				if ( unique )
					e = first;
			}
		}
		return e;
	}

	public Entry byId(int id) {
		return this.byId.get(id);
	}

	public int size() {
		return this.byId.size();
	}

	public Map<Integer, Entry> getEntries() {
		return Collections.unmodifiableMap(this.byId);
	}
}
