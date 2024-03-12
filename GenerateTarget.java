import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

public class GenerateTarget {
    public static void main(String[] args) throws IOException {
        /* 输入 */
        int tar_num = 100; /* 目标数量 */
        double rho = 0; /* 随机范围(°) */
        double t0 = 494395200; /* 星上时 */
        double[] r0 = {6878012.849e-3, -37902.952e-3, -16468.482e-3}; /* 位置(km) */
        double[] v0 = {44.818744e-3, 6235.960167e-3, 4366.073205e-3}; /* 速度(km/s) */
        double ti = t0 + 300; /* 起始时刻 */
        double tf = t0 + 900; /* 终止时刻 */
        double dt = tf - ti;


        /* 生成目标点 */
        double[] coe0 = OrbitPropagation.eci2coe(r0, v0);
        double t = 0., jd = 0.;
        double[] coe = new double[6], r = new double[3], v = new double[3], lla = new double[3];
        /* 存储 */
        Writer fw = new FileWriter("./目标参数.csv");
        BufferedWriter bw =new BufferedWriter(fw);
        bw.write("编号,经度,纬度");
        bw.newLine();
        for (int i = 0; i < tar_num; i++){
            t = ti + dt * i / tar_num;
            jd = EarthObservation.JD2009 + t / 86400.;
            TimeSystem.Utc utc = TimeSystem.jd2utc(jd);
            System.out.printf("%d\t%d\t%d\t%d\t%d\t%.2f\n", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
            coe = OrbitPropagation.orb_prop_j2(coe0, t - t0);
            OrbitPropagation.coe2eci(coe, r, v);
            r = EarthObservation.eci2ecef(jd, r);
            lla = EarthObservation.ecef2geod(r);
            lla[0] += Math.toRadians((Math.random() * 2. - 1.) * rho);
            lla[1] += Math.toRadians((Math.random() * 2. - 1.) * rho);
            lla[2] = 0.;
            bw.write(String.valueOf(i + 1));
            bw.write(",");
            bw.write(String.valueOf(Math.toDegrees(lla[0])));
            bw.write(",");
            bw.write(String.valueOf(Math.toDegrees(lla[1])));
            bw.newLine();
        }
        bw.close();
        fw.close();
    }

}
