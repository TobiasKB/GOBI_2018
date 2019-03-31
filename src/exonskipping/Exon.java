package exonskipping;

public class Exon extends Region {

	private String id;

	public Exon(int start, int end, String identifier) {
        super(start, end);
		this.id = identifier;
    }

	public String get_ID() {
		return this.id;
	}
}