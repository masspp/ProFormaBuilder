package ninja.mspp.proformabuilder;

/**
 * Raised for invalid ProForma strings and for modifications that cannot be resolved to a mass.
 * Extends {@link IllegalArgumentException}, the exception used by the other parsers of this library.
 */
public class ProFormaException extends IllegalArgumentException {

	private static final long serialVersionUID = 1L;

	private final int position;

	public ProFormaException(String message) {
		this(message, -1);
	}

	public ProFormaException(String message, int position) {
		super(position >= 0 ? message + " (at position " + position + ")" : message);
		this.position = position;
	}

	public ProFormaException(String message, Throwable cause) {
		super(message, cause);
		this.position = -1;
	}

	/** Character position in the input string, or -1 if not applicable. */
	public int getPosition() {
		return this.position;
	}
}
