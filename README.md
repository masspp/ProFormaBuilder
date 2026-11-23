# ProFormaBuilder

A Java library for parsing [ProForma 2.1](https://www.psidev.info/proforma) peptidoform strings
and calculating their masses, m/z values and fragment ions.

Masses are calculated with [GlycanMassUtility](https://gitlab.com/glycoinfo/ms/glycanmassutility).

## Features

* Parse ProForma 2.1 strings and write them back
* Calculate monoisotopic and average masses, m/z values and fragment ions (a/b/c/x/y/z and variants)
* Resolve modifications by Unimod, PSI-MOD and RESID names or accessions, mass deltas, formulas and glycan compositions
* The unmodified `unimod.xml` and `PSI-MOD.obo` are bundled (see [Bundled data](#bundled-data))

## Requirement

* Java 8 or higher
* Maven 3.5 or higher

## Build

```sh
$ mvn clean install
```

## Install

Add the dependency to your `pom.xml`.

```xml
<dependency>
  <groupId>ninja.mspp</groupId>
  <artifactId>proformabuilder</artifactId>
  <version>0.1.1</version>
</dependency>
```

The package is published to GitHub Packages, and GlycanMassUtility to the GitLab package registry.
Add both repositories to your `pom.xml`.

```xml
<repositories>
  <repository>
    <id>github</id>
    <url>https://maven.pkg.github.com/masspp/ProFormaBuilder</url>
  </repository>
  <repository>
    <id>GlycanMassUtility-gitlab-maven</id>
    <url>https://gitlab.com/api/v4/projects/43695885/packages/maven</url>
  </repository>
</repositories>
```

GitHub Packages requires authentication even for public packages.
Add a personal access token with the `read:packages` scope to your `~/.m2/settings.xml`.

```xml
<settings>
  <servers>
    <server>
      <id>github</id>
      <username>YOUR_GITHUB_USERNAME</username>
      <password>YOUR_PERSONAL_ACCESS_TOKEN</password>
    </server>
  </servers>
</settings>
```

## Usage

```java
import ninja.mspp.proformabuilder.FragmentIonType;
import ninja.mspp.proformabuilder.ProForma;

public class TestProForma {

	public static void main(String[] args) {
		ProForma pf = ProForma.parse("EM[Oxidation]EVN[Glycan:HexNAc2Hex3]ES[Phospho]PEK/2");

		System.out.println("Sequence: "+pf.getSequence());
		// -> EMEVNESPEK

		System.out.println("Monoisotopic: "+pf.getMonoisotopicMass());
		// -> 2178.79097956439

		System.out.println("Average: "+pf.getAverageMass());
		// -> 2180.05502799842

		System.out.println("m/z (charge from /2): "+pf.getMz());
		// -> 1090.4027662144

		System.out.println("b ions: "+pf.getFragments(FragmentIonType.B, 1));
		// -> [130.0498695210, 277.0852692290, 406.1278623178, ...]
	}
}
```

`ProForma` is an `IMassElement` of GlycanMassUtility, so it can also be used with `IonCloud` and `MassUtils`.
Use `ProForma.parseOrNull` to get `null` for blank input, as the parsers of GlycanMassUtility do.

The charge state of the string gives the ion for `getMz()`.
Following ProForma 2.1, a charge without charge carriers adds protons for a positive charge
but electrons for a negative charge.
Write the charge carriers to get the ions of negative ion mode.

| Ion | ProForma |
|---|---|
| [M+2H]<sup>2+</sup> | `PEPTIDE/2` |
| [M+2Na]<sup>2+</sup> | `PEPTIDE/[Na:z+1^2]` |
| [M+H+Na]<sup>2+</sup> | `PEPTIDE/[Na:z+1,H:z+1]` |
| [M-2H]<sup>2-</sup> | `PEPTIDE/[H-1:z-1^2]` |
| [M+2e]<sup>2-</sup> | `PEPTIDE/-2` |

## ProForma 2.1 Support

| Compliance level | Feature | Support |
|---|---|---|
| Level 1 | Standard amino acids (incl. U, O); Unimod and PSI-MOD names and accessions; mass deltas; N/C-terminal and labile modifications; several modifications per residue; INFO tags | Yes |
| Level 2 | Ambiguous residues B, Z (all possible masses), J, X; prefixed names and mass deltas; mass gaps; formulas with isotopes and negative counts; alternatives `[A\|B]`; unknown positions `?`, `^n`; ranges; sequence ambiguity `(?DQ)` | Yes |
| Level 2 | Groups `#g1` with localisation scores | Parsed and written; not assigned per site in fragment ions |
| Level 2 + top-down | RESID (resolved through PSI-MOD, narrowed down by the residue); names `(>...)`, `(>>...)` | Yes |
| Level 2 + cross-linking | `#XL`, `//`, `#BRANCH` | Parsed and written; masses when the linker is given by Unimod name, formula or mass; no fragment ions |
| Level 2 + cross-linking | XL-MOD names/accessions | No (parsed; resolving the mass raises an exception) |
| Level 2 + glycans | `Glycan:` compositions with the ProForma 2.1 monosaccharide names, braces components | Yes, except Sug, Tri and Tet (parsed; resolving the mass raises an exception) |
| Level 2 + glycans | GNOme | No (parsed; resolving the mass raises an exception) |
| Level 3 | Charged formulas; global isotopes; fixed modifications (incl. `N-term:Q`); chimeric strings; charge and charge carriers; ion notation | Yes |
| Level 3 | Placement controls (`Position:`, `Limit:`, `CoMKP`, `CoMUP`) | Parsed and written only |

Masses (monoisotopic and average; average masses are approximate for mass-only modifications), m/z and fragment ions
(a/b/c/x/y/z and variants, single chain) are calculated. Writing back keeps the modification tags as written
and outputs residues in upper case.

The monosaccharide names are resolved with the dictionaries of GlycanMassUtility.
Sug, Tri and Tet are not defined there: Sugose does not exist as a monosaccharide,
and trioses and tetroses rarely appear in glycans.

## Bundled Data

The unmodified `unimod.xml` (Unimod, Design Science License) and `PSI-MOD.obo` (PSI-MOD, CC BY 4.0) are bundled in
`src/main/resources/ninja/mspp/proformabuilder/cv/`.
Their source, version, retrieval date, SHA-256 and license are listed in `CV_SOURCES.md` in the same directory,
and the attribution is in [NOTICE](NOTICE).
Newer files can be used instead through the system properties `proformabuilder.unimod.xml` and `proformabuilder.psimod.obo`.

## License

BSD 3-Clause License. See [LICENSE](LICENSE).
The bundled data are distributed under their own licenses (see [NOTICE](NOTICE)).

## Release note
* 2026/11/23: released ver 0.1.1
  * Fixed the charge of negative charge carriers (e.g. `/[H-1:z-1^2]`) and negatively charged formulas
    (e.g. `[Formula:Cl:z-1]`), which was 0 or wrong
* 2026/11/23: released ver 0.1.0
  * Parse and write ProForma 2.1 strings, and calculate masses, m/z values and fragment ions
    (based on the ProForma support contributed for GlycanMassUtility)
  * Unimod modifications are looked up by their titles only, not by their full names
  * `ProForma.getMz(int)` removes protons for a negative charge, unlike the charge state `/-2` of the string
  * Added `equals` and `hashCode` to `ProForma`
