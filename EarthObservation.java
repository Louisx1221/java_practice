public class EarthObservation {
    static double RE = 6378.1365, JD2009 = 2454832.5;
    public static void main(String[] args) {
        TimeSystem.Utc utc = new TimeSystem.Utc(2024, 1, 22, 10, 30, 0.);
        double jd = TimeSystem.utc2jd(utc);
        double[] r_eci = {4806.000777,  3980.971899,    2891.883100};
        double[] v_eci = {-5.444155,    4.415641,       2.969023};
        double[] r_ecef = eci2ecef(jd, r_eci);
        double[] lla = ecef2geod(r_ecef);

        /* 星下点 */
        System.out.print("Ground track\n");
        System.out.printf("Longtitude:\t%.5f \tdegrees\n", Math.toDegrees(lla[0]));
        System.out.printf("Latitude:  \t%.5f \tdegrees\n", Math.toDegrees(lla[1]));
        System.out.printf("Altitude:  \t%.3f \tkm\n", lla[2]);

        /* 欧拉角 */
        double[] lla_tar = {Math.toRadians(121.532), Math.toRadians(25.0478), 0};
        double[] ang = observe_angle(jd, r_eci, v_eci, lla_tar);
        System.out.print("Eular angle\n");
        System.out.printf("Roll angle: \t%.5f \tdegrees\n", Math.toDegrees(ang[0]));
        System.out.printf("Pitch angle:\t%.5f \tdegrees\n", Math.toDegrees(ang[1]));
        System.out.printf("Yaw angle:  \t%.5f \tdegrees\n", Math.toDegrees(ang[2]));
    }

    public static double[] eci2ecef(double jd, double[] r_eci) {
        /* ECI位置转ECEF */
        double[][] rot_mat = rotation_matrix(jd);
        double[] r_ecef = matrix_times_vector(rot_mat, r_eci);
        return r_ecef;
    }

    public static double[] ecef2eci(double jd, double[] r_ecef) {
        /* ECEF位置转 ECI */
        double[][] rot_mat = rotation_matrix(jd);
        rot_mat = matrix_transpose(rot_mat);
        double[] r_eci = matrix_times_vector(rot_mat, r_ecef);
        return r_eci;
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

    public static double[] geod2ecef(double lla[]) {
        /* 大地经纬高转ECEF位置 */
        double lon = lla[0], lat = lla[1], alt = lla[2];
        double e2 = 0.00669;
        double v = RE / Math.sqrt(1 - e2 * Math.sin(lat) * Math.sin(lat));
        double[] r_ecef = new double[3];
        r_ecef[0] = (v + alt) * Math.cos(lat) * Math.cos(lon);
        r_ecef[1] = (v + alt) * Math.cos(lat) * Math.sin(lon);
        r_ecef[2] = ((1.0 - e2) * v + alt) * Math.sin(lat);
        return r_ecef;
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

    public static double[] eci2rtn(double[] r_eci, double[] v_eci, double[] vec_eci) {
        /* ECI向量转RTN(RSW) */
        double[] r_unit = unitize(r_eci);
        double[] n_unit = cross(r_eci, v_eci);
        n_unit = unitize(n_unit);
        double[] t_unit = cross(n_unit, r_unit);
        double[][] rot_mat = horzcat(r_unit, t_unit, n_unit);
        double[] vec_rtn = vector_times_matrix(vec_eci, rot_mat);
        return vec_rtn;
    }

    public static double[] observe_angle(double jd, double[] r_eci, double[] v_eci, double[] lla)
    {
        double[] r_ecef_tar = geod2ecef(lla);
        double[] r_eci_tar = ecef2eci(jd, r_ecef_tar);
        double[] dr_eci = new double[3];
        for (int i = 0; i < 3; i++)
            dr_eci[i] = r_eci[i] - r_eci_tar[i];
        double[] dr_rtn = eci2rtn(r_eci, v_eci, dr_eci);
        double[] dr_body = rtn2body(dr_rtn);
        double[] obs_ang = eular_angle(dr_body);
        return obs_ang;
    }

    public static double[] rtn2body(double[] vec_rtn) {
        /* RTN坐标系转本体坐标系 */
        double[] vec_body = {vec_rtn[1], vec_rtn[2], vec_rtn[0]};
        return vec_body;
    }

    public static double[] eular_angle(double[] vec) {
        /* 欧拉角 */
        double[] ang = new double[3];
        ang[0] = Math.asin(vec[1] / norm(vec)); /* roll */
        ang[1] = -Math.atan2(vec[2], vec[0]) + Math.PI / 2; /* pitch */
        return ang;
    }

    public static double[] matrix_times_vector(double[][] mat, double[] vin) {
        /* 矩阵乘向量 */
        double[] vout = new double[3];
        for (int i = 0; i < 3; i++)
            vout[i] = mat[i][0] * vin[0] + mat[i][1] * vin[1] + mat[i][2] * vin[2];
        return vout;
    }

    public static double[] vector_times_matrix(double[] vin, double[][] mat) {
        /* 向量乘矩阵 */
        double[] vout = new double[3];
        for (int i = 0; i < 3; i++)
            vout[i] = mat[0][i] * vin[0] + mat[1][i] * vin[1] + mat[2][i] * vin[2];
        return vout;
    }

    public static double[][] matrix_transpose(double[][] mat_in){
        /* 矩阵转置 */
        double[][] mat_out = new double[3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                mat_out[i][j] = mat_in[j][i];
        }
        return mat_out;
    }

    public static double[] cross(double vin1[], double vin2[]) {
        /* 向量叉乘 */
        double[] vout = new double[3];
        vout[0] = vin1[1] * vin2[2] - vin1[2] * vin2[1];
        vout[1] = vin1[2] * vin2[0] - vin1[0] * vin2[2];
        vout[2] = vin1[0] * vin2[1] - vin1[1] * vin2[0];
        return vout;
    }

    public static double dot(double vin1[], double vin2[]) {
        /* 向量点乘 */
        double out = 0;
        for (int i = 0; i < 3; i++)
            out += vin1[i] * vin2[i];
        return out;
    }

    public static double norm(double v[]) {
        /* 向量的模 */
        return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    public static double[] unitize(double vin[]) {
        /* 向量单位化 */
        double[] vout = new double[3];
        double vm = norm(vin);
        for (int i = 0; i < 3; i++)
            vout[i] = vin[i] / vm;
        return vout;
    }

    public static double[][] horzcat(double vec0[], double vec1[], double vec2[]) {
        /* 串联矩阵 */
        double[][] mat_out = new double[3][3];
        for (int i = 0; i < 3; i++)
        {
            mat_out[i][0] = vec0[i];
            mat_out[i][1] = vec1[i];
            mat_out[i][2] = vec2[i];
        }
        return mat_out;
    }
}
