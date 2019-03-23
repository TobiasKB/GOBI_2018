package utils;

import java.util.HashSet;
import java.util.Set;

/**
 * @author TKB
 */
public final class HashUtil {

    /**
     *
     */
    private HashUtil() {
    }

    public static <T> T get_nthElement_Set(Set<T> set, int n) {
        if (null != set && n >= 0 && n < set.size()) {
            int count = 0;
            for (T element : set) {
                if (n == count)
                    return element;
                count++;
            }
        }
        return null;
    }

    public static <T> HashSet<T> mergeSet(HashSet<T> a, HashSet<T> b) {
        return new HashSet<T>() {
            {
                addAll(a);
                addAll(b);
            }
        };

    }
}