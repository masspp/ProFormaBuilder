/**
 * ProForma 2.1 parser and mass calculation.
 * <p>
 * Entry point: {@link ninja.mspp.proformabuilder.ProForma#parse(String)}.
 * </p>
 * <h2>Supported</h2>
 * <ul>
 *   <li>Level 1: residues (including U and O), Unimod and PSI-MOD names and accessions, mass deltas,
 *       N/C-terminal modifications, labile modifications {@code {...}}, several modifications per residue,
 *       {@code INFO} tags.</li>
 *   <li>Level 2: ambiguous residues B, Z (all possible masses), J and X; prefixed names and mass deltas;
 *       formulas with isotopes and negative counts; alternatives {@code [A|B]} (the first descriptor with a
 *       mass is used); modifications of unknown position {@code ?}, {@code ^n}; groups {@code #g1} with scores;
 *       ranges; sequence ambiguity {@code (?DQ)}.</li>
 *   <li>Level 2 + glycans: {@code Glycan:} compositions with the ProForma 2.1 monosaccharide names
 *       (section 10.2 table and formal grammar; other names are errors) and braces components.</li>
 *   <li>Level 2 + top-down: names {@code (>...)}, RESID names and accessions (resolved through PSI-MOD and
 *       narrowed down by the residue; an exception is raised if the mass is not unique).</li>
 *   <li>Level 3: charge {@code /2} (negative charges add electrons), charge carriers {@code /[Na:z+1^2]},
 *       charged formulas, fixed modifications (including {@code N-term:Q}), global isotopes ({@code <13C>},
 *       {@code <D>}), chimeric strings {@code +}, ion notation ({@code -[b-type-ion]}).</li>
 *   <li>Fragment ions a/b/c/x/y/z (and c-dot, z-dot, ...) for a single chain; labile modifications are
 *       excluded by default.</li>
 *   <li>Writing back ({@link ninja.mspp.proformabuilder.ProFormaWriter}): modification tags are
 *       written as they were read and residues in upper case (ProForma is case insensitive except for CV names and
 *       formulas), so parsing the written string gives the same structure.</li>
 * </ul>
 * <h2>Not supported</h2>
 * <ul>
 *   <li>XL-MOD and GNOme names/accessions (parsed, but resolving their mass raises an exception);
 *       fragment ions of cross-linked chains.</li>
 *   <li>The monosaccharide names Sug, Tri and Tet, which are not defined in GlycanMassUtility
 *       (parsed, but resolving their mass raises an exception).</li>
 *   <li>Placement controls ({@code Position:}, {@code Limit:}, {@code CoMKP}, {@code CoMUP}) are parsed only.</li>
 *   <li>Per-site assignment of {@code #g} groups to fragment ions; normalisation of glycan names on output;
 *       an output mode restricted to a compliance level.</li>
 * </ul>
 * <h2>Data</h2>
 * <p>
 * The unmodified {@code unimod.xml} and {@code PSI-MOD.obo} are bundled as resources of the package
 * {@code ninja.mspp.proformabuilder.cv}; {@code CV_SOURCES.md} next to them records their source URL, version, retrieval date,
 * SHA-256 and license (the Design Science License text of Unimod is bundled as well). Newer CV files can be
 * used instead through the system properties of
 * {@link ninja.mspp.proformabuilder.cv.ModificationResolver}.
 * </p>
 */
package ninja.mspp.proformabuilder;
