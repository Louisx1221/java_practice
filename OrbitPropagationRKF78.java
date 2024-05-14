public class OrbitPropagationRKF78 {
    /* constant */
    static double MU = 398600.4415, RE = 6378.1365, J2 = 1.082627e-3;
    static double[] ch = {0., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 0., 41. / 840., 41. / 840.};
    static double[] alpha = {0., 2. / 27., 1. / 9., 1. / 6., 5. / 12., 0.5, 5. / 6., 1. / 6., 2. / 3., 1. / 3., 1., 0., 1.};
    static double[][] beta = {
        {0.,            0.,         0.,         0.,             0.,             0.,             0.,             0.,         0.,         0.,         0., 0., },
        {2. / 27.,      0.,         0.,         0.,             0.,             0.,             0.,             0.,         0.,         0.,         0., 0., },
        {1. / 36.,      1. / 12.,   0.,         0.,             0.,             0.,             0.,             0.,         0.,         0.,         0., 0., },
        {1. / 24.,      0.,         1. / 8.,    0.,             0.,             0.,             0.,             0.,         0.,         0.,         0., 0., },
        {5. / 12.,      0.,         -25. / 16., 25. / 16.,      0.,             0.,             0.,             0.,         0.,         0.,         0., 0., },
        {0.05,          0.,         0.,         0.25,           0.2,            0.,             0.,             0.,         0.,         0.,         0., 0., },
        {-25. / 108.,   0.,         0.,         125. / 108.,    -65. / 27.,     125. / 54.,     0.,             0.,         0.,         0.,         0., 0., },
        {31. / 300.,    0.,         0.,         0.,             61. / 225.,     -2. / 9.,       13. / 900.,     0.,         0.,         0.,         0., 0., },
        {2.,            0.,         0.,         -53. / 6.,      704. / 45.,     -107. / 9.,     67. / 90.,      3.,         0.,         0.,         0., 0., },
        {-91. / 108.,   0.,         0.,         23. / 108.,     -976. / 135.,   311. / 54.,     -19. / 60.,     17. / 6.,   -1. / 12.,  0.,         0., 0., },
        {2383. / 4100., 0.,         0.,         -341. / 164.,   4496. / 1025.,  -301. / 82.,    2133. / 4100.,  45. / 82.,  45. / 164., 18. / 41.,  0., 0., },
        {3. / 205.,     0.,         0.,         0.,             0.,             -6. / 41.,      -3. / 205.,     -3. / 41.,  3. / 41.,   6. / 41.,   0., 0., },
        {-1777. / 4100.,0.,         0.,         -341. / 164.,   4496. / 1025.,  -289. / 82.,    2193. / 4100.,  51. / 82.,  33. / 164., 12. / 41.,  0., 1., },
    };

    @FunctionalInterface
    interface DynamicalModel {
        double[] differentialEquation (double[] x);
    }

    /* 二体模型 */
    public static DynamicalModel orbPropTwoBody = x -> {
        double rm = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        double r3 = rm * rm * rm;

        double[] dx = new double[6];
        dx[0] = x[3];
        dx[1] = x[4];
        dx[2] = x[5];
        dx[3] = -MU * x[0] / r3;
        dx[4] = -MU * x[1] / r3;
        dx[5] = -MU * x[2] / r3;

        return dx;
    };

    /* J2摄动模型 */
    public static DynamicalModel orbPropJ2Pert = x -> {
        double rm = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        double r2 = rm * rm, r3 = rm * rm * rm, rz2 = x[2] * x[2], re2 = RE * RE;

        double[] dx = new double[6];
        dx[0] = x[3];
        dx[1] = x[4];
        dx[2] = x[5];
        dx[3] = -MU * x[0] / r3 * (1 + (1.5 * J2 * re2 / r2 * (1 - 5 * rz2 / r2)));
        dx[4] = -MU * x[1] / r3 * (1 + (1.5 * J2 * re2 / r2 * (1 - 5 * rz2 / r2)));
        dx[5] = -MU * x[2] / r3 * (1 + (1.5 * J2 * re2 / r2 * (3 - 5 * rz2 / r2)));

        return dx;
    };

    public static double[] RKF78(DynamicalModel dyn_mdl, double[] x0, double ti, double tf, double h, int neq) {
        /**
         * solve first ODE system Runge-Kutta-Fehlberg 7(8) method
         * input
         * dyn_mdl.deq = name of function which defines the system of differential equations
         * neq   = number of differential equations
         * ti    = initial simulation time
         * tf    = final simulation time
         * h     = initial guess for integration step size
         * tetol = truncation error tolerance (non-dimensional)
         * x     = integration vector at time = ti
         * output
         * xout  = integration vector at time = tf
         */

        double dt, twrk, tetol = 1e-12, xerr, ter, tol, tconst, sdt;
        double[] xwrk = new double[7], xdot = new double[7];
        double[][] f = new double[7][13];
        int i, j, k;

        double[] x = new double[neq];
        for (i = 0; i < neq; i++)
            x[i] = x0[i];

        //compute integration "direction"
        sdt = Math.signum(tf - ti);
        dt = Math.abs(h) * sdt;

        while (true)
        {
            // load "working" time and integration vector
            twrk = ti;
            for (i = 0; i < neq; i++)
                xwrk[i] = x[i];
            // check for last dt
            if (Math.abs(dt) > Math.abs(tf - ti))
                dt = tf - ti;
            if (Math.abs(ti - tf) < 1e-8)
            {
                break;
            }
            // evaluate equations of motion
            xdot = dyn_mdl.differentialEquation(x);
            for (i = 0; i < neq; i++)
                f[i][0] = xdot[i];
            // compute solution
            for (k = 1; k < 13; k++)
            {
                for (i = 0; i < neq; i++)
                {
                    x[i] = 0.;
                    for (j = 0; j < k; j++)
                        x[i] += beta[k][j] * f[i][j];
                    x[i] *= dt;
                    x[i] += xwrk[i];
                }
                ti = twrk + alpha[k] * dt;
                xdot = dyn_mdl.differentialEquation(x);
                for (j = 0; j < neq; j++)
                    f[j][k] = xdot[j];
            }
            for (i = 0; i < neq; i++)
            {
                x[i] = 0.;
                for (j = 0; j < 13; j++)
                    x[i] += ch[j] * f[i][j];
                x[i] *= dt;
                x[i] += xwrk[i];
            }

            // truncation error calculations
            xerr = tetol;
            for (i = 0; i < neq; i++)
            {
                ter = Math.abs((f[i][0] + f[i][10] - f[i][11] - f[i][12]) * ch[11] * dt);
                tol = Math.abs(x[i]) * tetol + tetol;
                tconst = ter / tol;
                if (tconst > xerr)
                    xerr = tconst;
            }

            // compute new step size
            //dt = 0.8 * dt * (1 / xerr) ^ (1 / 8);
            dt *= 0.8 * Math.pow(1. / xerr, 1. / 8.);
            if (xerr > 1)
            {
                // reject current step
                ti = twrk;
                for (i = 0; i < neq; i++)
                    x[i] = xwrk[i];
            }
        }

        return x;
    }

    public static double[][] orbitPropagation(double[] r0, double[] v0, double dt, double tspan, int model) {
        /**
         * r0[3]    = ECI位置向量(km)
         * v0[3]    = ECI速度向量(km/s)
         * dt       = 递推时长(s)
         * tspan    = 时间间隔(s)
         * model    = 0(二体模型), 1(J2摄动模型)
         */

        double[] x = {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]};
        int N = (int)Math.floor(dt / tspan);
        double[][] xout = new double[N + 1][7];

        for (int i = 0; i < N; i++)
        {
            xout[i][0] = i * tspan;
            for (int j = 0; j < 6; j++)
                xout[i][j + 1] = x[j];
            if (model == 1)
                x = RKF78(orbPropJ2Pert, x, 0, tspan, 1e-2, 6);
            else
                x = RKF78(orbPropTwoBody, x, 0, tspan, 1e-2, 6);
        }
        xout[N][0] = N * tspan;
        for (int j = 0; j < 6; j++)
            xout[N][j + 1] = x[j];

        return xout;
    }

    public static void main(String[] args) {
        // 初始轨道参数
        double[] r0 = {-6121.222491222, -2258.448349219, 2189.751497198};   // km
        double[] v0 = {-2.776819144, 0.596247651, -7.145828963};            // km/s

        // 递推时长及时间间隔
        double dt = 6000, tspan = 1; // 递推时长, 时间间隔(s)

        double[][] trv_mat = orbitPropagation(r0, v0, dt, tspan, 1);

        int N = trv_mat.length - 1;
        double tf = trv_mat[N][0];
        double[] rf = {trv_mat[N][1], trv_mat[N][2], trv_mat[N][3]};
        double[] vf = {trv_mat[N][4], trv_mat[N][5], trv_mat[N][6]};

        System.out.printf("time: %.2f s\t", tf);
        System.out.printf("position: [%.3f, %.3f, %.3f] km\t", rf[0], rf[1], rf[2]);
        System.out.printf("velocity: [%.3f, %.3f, %.3f] km/s\n", vf[0], vf[1], vf[2]);

    }
}
