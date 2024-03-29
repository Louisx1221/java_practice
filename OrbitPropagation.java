public class OrbitPropagation {
    static double MU = 398600.4415, RE = 6378.1365, J2 = 1.082627e-3;
    public static void main(String[] args) {
        /* 1. 轨道六根数 */
        double[] r0 = {-6121.222491222, -2258.448349219, 2189.751497198};
        double[] v0 = {-2.776819144, 0.596247651, -7.145828963};
        double[] coe0 = eci2coe(r0, v0);
        for (int i = 0; i < 2; i++)
            System.out.printf("%.6f\t", coe0[i]);
        for (int i = 2; i < 6; i++)
            System.out.printf("%.6f\t", Math.toDegrees(coe0[i]));
        System.out.print("\n");
        // STK: 7029.280916 0.020921 100.805 16.579 160.923 0.177

        /* 2. 轨道递推 */
        double dt = 600;
        double[] coe = new double[6];
        coe = orb_prop_2body(coe0, dt);
        for (int i = 0; i < 2; i++)
            System.out.printf("%.6f\t", coe[i]);
        for (int i = 2; i < 6; i++)
            System.out.printf("%.6f\t", Math.toDegrees(coe[i]));
        System.out.print("\n");
        // STK: 7029.280916 0.020921 100.805 16.579 160.923 38.471

        coe = orb_prop_j2(coe0, dt);
        for (int i = 0; i < 2; i++)
            System.out.printf("%.6f\t", coe[i]);
        for (int i = 2; i < 6; i++)
            System.out.printf("%.6f\t", Math.toDegrees(coe[i]));
        System.out.print("\n");
        // STK: 7029.280916 0.020921 100.805 16.588 160.903 38.448

        /* 3. ECI坐标系位置速度矢量 */
        double[] r = new double[3], v = new double[3];
        coe2eci(coe, r, v);
        for (int i = 0; i < 3; i++)
            System.out.printf("%.6f\t", r[i]);
        for (int i = 0; i < 3; i++)
            System.out.printf("%.6f\t", v[i]);
        System.out.print("\n");
        // STK: -6372.419981 -1448.181730  -2254.818593 1.959712 2.002236 -7.124916
    }

    public static double[] orb_prop_2body(double coe0[], double dt) {
        /* 轨道六根数轨道递推(二体模型) */
        double a, e, f, E, M, n;
        double[] coe = new double[6];

        /* 初始轨道六根数 */
        a = coe0[0];
        e = coe0[1];
        f = coe0[5];
        n = Math.sqrt(MU / a / a / a);
        E = Math.atan2(Math.sin(f) * Math.sqrt(1 - e * e), e + Math.cos(f));
        M = E - e * Math.sin(E);

        /* 终止轨道六根数 */
        M += n * dt;
        E = solve_kepler(M, e);
        f = 2 * Math.atan(Math.sqrt((1 + e) / (1 - e)) * Math.tan(E / 2));
        for (int i = 0; i < 5; i++)
            coe[i] = coe0[i];
        coe[5] = f % Math.TAU;
        return coe;
    }

    public static double[] orb_prop_j2(double coe0[], double dt) {
        /* 轨道六根数轨道递推(线性J2摄动模型) */
        double a, e, f, E, M, n, p, cosi, CJ2;
        double[] coe = new double[6];

        /* 初始轨道六根数 */
        a = coe0[0];
        e = coe0[1];
        f = coe0[5];
        n = Math.sqrt(MU / a / a / a);
        p = a * (1 - e * e);
        E = Math.atan2(Math.sin(f) * Math.sqrt(1 - e * e), e + Math.cos(f));
        M = E - e * Math.sin(E);
        cosi = Math.cos(coe0[2]);

        /* 终止轨道六根数 */
        CJ2 = J2 * n * (RE / p) * (RE / p);
        M += (n + 0.75 * CJ2 * Math.sqrt(1 - e * e) * (3 * cosi * cosi - 1)) * dt;
        E = solve_kepler(M, e);
        f = 2 * Math.atan(Math.sqrt((1 + e) / (1 - e)) * Math.tan(E / 2));
        for (int i = 0; i < 5; i++)
            coe[i] = coe0[i];
        coe[3] += -1.5 * CJ2 * cosi * dt;
        coe[4] += 0.75 * CJ2 * (5 * cosi * cosi - 1)* dt;
        coe[5] = f % Math.TAU;
        for (int i = 3; i < 6; i++)
            coe[i] += (coe[i] < 0) ? Math.TAU : 0;
        return coe;
    }

    public static double solve_kepler(double M, double e) {
        /* 解开普勒方程 */
        int N = 10, i = 0;
        double E = M % Math.TAU, v = 1, tol = 1e-9;

        while (Math.abs(v) > tol && i < N)
        {
            v = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
            E -= v;
            i++;
        }
        return E;
    }

    public static double[] eci2coe(double r[], double v[]) {
        /* ECi位置速度转轨道六根数 */
        double[] h_, r_, ecc_, f_ = new double[3], g_ = new double[3];
        double sma, ecc, inc, raan, argper, tanom, p, q, const1, h, xk, tlon;
        int i;

        h_ = LinearAlgebra.cross(r, v); /* angular momentum vector */
        ecc_ = LinearAlgebra.cross(v, h_); /* eccentricity vector */
        r_ = LinearAlgebra.unitize(r);
        for (i = 0; i < 3; i++)
            ecc_[i] = ecc_[i] / MU - r_[i];

        sma = 1 / (2 / LinearAlgebra.norm(r) - LinearAlgebra.dot(v, v) / MU); /* semi-major axis */

        h_ = LinearAlgebra.unitize(h_);
        p = h_[0] / (1 + h_[2]);
        q = -h_[1] / (1 + h_[2]);

        const1 = 1 / (1 + p * p + q * q);
        f_[0] = const1 * (1 - p * p + q * q);
        f_[1] = const1 * 2 * p * q;
        f_[2] = -const1 * 2 * p;
        g_[0] = const1 * 2 * p * q;
        g_[1] = const1 * (1 + p * p - q * q);
        g_[2] = const1 * 2 * q;
        h = LinearAlgebra.dot(ecc_, g_); /* angular momentum*/
        xk = LinearAlgebra.dot(ecc_, f_);

        ecc = Math.sqrt(h * h + xk * xk); /* eccentricity */
        ecc = Math.clamp(ecc, 0, 0.97);

        inc = 2 * Math.atan(Math.sqrt(p * p + q * q)); /* inclination */
        tlon = Math.atan2(LinearAlgebra.dot(r, g_), LinearAlgebra.dot(r, f_)); /* true longitude */

        /* check for equatorial orbit*/
        if (inc > 0.00000001)
            raan = Math.atan2(p, q); /* right ascension of the ascending node */
        else
            raan = 0;

        /* check for circular orbit*/
        if (ecc > 0.00000001)
            argper = (Math.atan2(h, xk) - raan) % Math.TAU; /* argument of perigee */
        else
            argper = 0;

        tanom = (tlon - raan - argper) % Math.TAU; /* true anomaly */

        /* singular value when h_(3)==-1*/
        if (1 - h_[2] < 0.00000001)
        {
            inc = Math.PI;
            raan = Math.PI;
            argper = (-Math.atan2(h, xk) + raan) % Math.TAU;
        }

        double[] coe = {sma, ecc, inc, raan, argper, tanom};
        for (i = 3; i < 6; i++)
            coe[i] += (coe[i] < 0) ? Math.TAU : 0;
        return coe;
    }

    public static void coe2eci(double coe[], double r[], double v[])
    {
        /* 轨道六根数转ECI位置速度 */
        double sma, ecc, inc, raan, argper, tanom, slr, rm, arglat, sarglat, carglat, c4, c5, c6;
        double sinc, cinc, sraan, craan;

        sma = coe[0];
        ecc = coe[1];
        inc = coe[2];
        raan = coe[3];
        argper = coe[4];
        tanom = coe[5];

        slr = sma * (1.0 - ecc * ecc);
        rm = slr / (1.0 + ecc * Math.cos(tanom));
        arglat = argper + tanom;
        sarglat = Math.sin(arglat);
        carglat = Math.cos(arglat);

        c4 = Math.sqrt(MU / slr);
        c5 = ecc * Math.cos(argper) + carglat;
        c6 = ecc * Math.sin(argper) + sarglat;
        sinc = Math.sin(inc);
        cinc = Math.cos(inc);
        sraan = Math.sin(raan);
        craan = Math.cos(raan);

        /* position vector and velocity vector*/
        r[0] = rm * (craan * carglat - sraan * sarglat * cinc);
        r[1] = rm * (sraan * carglat + craan * sarglat * cinc);
        r[2] = rm * sarglat * sinc;
        v[0] = -c4 * (c6 * craan + c5 * sraan * cinc);
        v[1] = -c4 * (c6 * sraan - c5 * craan * cinc);
        v[2] = c4 * c5 * sinc;
    }
}
