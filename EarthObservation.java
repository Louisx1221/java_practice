public class EarthObservation {
    static double RE = 6378.1365, JD2009 = 2454832.5;
    public static void main(String[] args) {
        double t = 497894421.5;
        double jd = JD2009 + t / 86400.0;
        double[] r_eci = {-6121.222491222, -2258.448349219, 2189.751497198};
        double[] r_ecef = eci2ecef(jd, r_eci);
        double[] lla = ecef2geod(r_ecef);
        System.out.printf("Longtitude: %.5f degrees\n", Math.toDegrees(lla[0]));
        System.out.printf("Latitude:   %.5f degrees\n", Math.toDegrees(lla[1]));
        System.out.printf("Altitude:   %.3f km\n", lla[2]);
    }

    public static double[] eci2ecef(double jd, double[] r_eci) {
        /* ECI位置转ECEF */
        double[][] rot_mat = rotation_matrix(jd);
        double[] r_ecef = matrix_times_vector(rot_mat, r_eci);
        return r_ecef;
    }

    public static double[] ecef2geod(double r_ecef[]) {
        /* ECEF位置转大地经纬高 */
        double e2 = 0.00669, rho = Math.sqrt(r_ecef[1] * r_ecef[1] + r_ecef[0] * r_ecef[0]), v = 0;
        double[] lla = new double[3];
        lla[0] = Math.atan2(r_ecef[1], r_ecef[0]);
        lla[1] = Math.atan2(r_ecef[2], rho);
        for (int j = 0; j < 5; j++)
        {
            v = RE / Math.sqrt(1.0 - e2 * Math.sin(lla[1]) * Math.sin(lla[1]));
            lla[1] = Math.atan2(r_ecef[2] + e2 * v * Math.sin(lla[1]), rho);
        }
        lla[2] = r_ecef[0] / Math.cos(lla[0]) / Math.cos(lla[1]) - v;
        return lla;
    }

    public static double[][] rotation_matrix(double jd) {
        /* 地球自转矩阵 */
        double jt = (jd - 2451545.0) / 36525.0;
        double thast_rot = 67310.54841 + (876600 * 3600.0 + 8640184.812866) * jt + 0.0093104 * jt * jt - 0.0000062 * jt * jt * jt;
        thast_rot %= 86400.0;
        thast_rot /= 240.0;
        thast_rot = Math.toRadians(thast_rot);
        double[][] rot = {
            {Math.cos(thast_rot), Math.sin(thast_rot), 0},
            {-Math.sin(thast_rot), Math.cos(thast_rot), 0},
            {0, 0, 1}
        };
        return rot;
    }

    public static double[] matrix_times_vector(double[][] mat, double[] vin) {
        /* 矩阵乘向量 */
        double[] vout = new double[3];
        for (int j = 0; j < 3; j++)
            vout[j] = mat[j][0] * vin[0] + mat[j][1] * vin[1] + mat[j][2] * vin[2];
        return vout;
    }

}