public class LinearAlgebra {
    public static double[] matrix_times_vector(double[][] mat, double[] vin) {
        /* 矩阵乘向量 */
        double[] vout = new double[3];
        for (int i = 0; i < 3; i++)
            vout[i] = mat[i][0] * vin[0] + mat[i][1] * vin[1] + mat[i][2] * vin[2];
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
}
