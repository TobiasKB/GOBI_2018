package exonskipping;

import utils.TesException;

import java.util.TreeSet;

/**
 * @author TKB
 */
public class RegionVector implements Comparable<RegionVector> {

    /**
     * size: wie viele Regionen wurden verschmolzen
     */
    public RegionVector inverse;
    public TreeSet<Region> regions;
    public int length;
    public int size;
    public int start;
    public int end;


    <T> RegionVector(T regionSet) {
        this.regions = (TreeSet) regionSet;
        this.length = 0;
        this.size = regions.size();
        calculate_length();
    }


    /*
     * Reduziert ein Set von Regionen auf Regionen von einem bestimmten Start bis zu einem gegebenen Stop
     * */
    public <T> RegionVector(T regionSet, int newstart, int newstop) {
        this.regions = new TreeSet<>();
//        this.length = 0;
        this.size = 0;
        try {
            for (Object r : (TreeSet) regionSet) {
                if (((Region) r).start >= newstart && ((Region) r).end < newstop) {
                    regions.add((Region) r);
                    this.start = Math.min(((Region) r).start, this.start);
                    this.end = Math.max(((Region) r).end, this.end);
                }
            }
            if (regions.size() == 0) {
                this.size = 0;
                this.start = newstart;
                this.end = newstop;
                this.length = 0;

            } else {
                this.size = regions.size();
                calculate_length();

            }
        } catch (NullPointerException e) {
            new TesException("Fehler beim einlesen eines RegionVectors" + e);
        }
    }


    private void calculate_length() {
        int temp_length = 0;
        for (Region r : regions) {
            temp_length += r.getLength();
        }
        temp_length += this.size;

        this.length = temp_length;

    }


    public void calculate_invers() {

    }

    /**
     * @param region
     */
    public void merge(Region region) {

    }

    @Override
    public int compareTo(RegionVector o) {
        return 0;
    }
}