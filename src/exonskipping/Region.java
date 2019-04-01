package exonskipping;

/**
 * @author TKB
 */
public class Region implements Comparable<Region> {

    /**
     *
     */
    public final int start;
    public final int end;
    public final int length;

    public Region(int start, int end) {
        this.start = start;
        this.end = end;
        this.length = end - start;
    }

    public int getEnd() {
        return end;
    }

    public int getStart() {
        return start;
    }

    public String toString() {
        return start + ":" + end;
    }

    @Override
    public int hashCode() {
        return ((start * 104723) % 104729) + ((end * 104717) % 104711);
    }

    /* (non-Javadoc)
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    @Override
    public int compareTo(Region other) {

/*
        if (this.start != other.start) {
            return this.start - other.start;
        } else {
            return this.end - other.end;
        }
*/


        Region otherr = other;
        if (this.start > otherr.start) {
            return 1;
        } else if (this.end < otherr.start) {
            return -1;
        } else {
            return 0;
        }
    }

    public boolean equals(Region other) {
        return this.start == other.start && this.end == other.end;
    }

    public int get_skippedBases(Region other) {
        return other.start - this.end;
    }


    public int getLength() {
        return this.length;
    }

}