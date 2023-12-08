public class CodeHex {
    public static void main(String[] args) {
        double lon = 131.27, lat = -17.87;
        byte[] buf = new byte[8];
        long temps;

        /* 编码 */
        temps = (long) (lon * 1e5);
        buf[0] = (byte) (temps >> 24 & 0xFF);
        buf[1] = (byte) (temps >> 16 & 0xFF);
        buf[2] = (byte) (temps >> 8 & 0xFF);
        buf[3] = (byte) (temps & 0xFF);
        temps = (long) (lat * 1e5);
        buf[4] = (byte) (temps >> 24 & 0xFF);
        buf[5] = (byte) (temps >> 16 & 0xFF);
        buf[6] = (byte) (temps >> 8 & 0xFF);
        buf[7] = (byte) (temps & 0xFF);
        for (int i = 0; i < 4; i++)
        {
            System.out.printf("%02X%02X ", buf[2 * i], buf[2 * i + 1]);
        }

        /* 解码 */
        temps  = (buf[0] & 0xFF) << 24 | (buf[1] & 0xFF) << 16 | (buf[2] & 0xFF) << 8 | (buf[3] & 0xFF);
        lon = temps * 1e-5;
        temps  = (buf[4] & 0xFF) << 24 | (buf[5] & 0xFF) << 16 | (buf[6] & 0xFF) << 8 | (buf[7] & 0xFF);
        lat = temps * 1e-5;
        System.out.printf("\nlon: %.2f, lat: %.2f", lon, lat);
    }
}
