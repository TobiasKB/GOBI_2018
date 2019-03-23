package utils;

/**
 * @author TKB
 */

public class TesException extends RuntimeException {

    public TesException(String message, Throwable err) {
        super(message, err);

    }

    public TesException(String message) {
        super(message);

    }


}