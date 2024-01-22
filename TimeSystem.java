public class TimeSystem {
    static double JD2009 = 2454832.5;
    public static void main(String[] args) {
        /* 1. UTC时间 */
        double t0 = 497894421.5;
        double jd0 = JD2009 + t0 / 86400.0;
        Utc utc0 = jd2utc(jd0);
        System.out.printf("%d\t%d\t%d\t%d\t%d\t%.2f\n", utc0.year, utc0.month, utc0.day, utc0.hour, utc0.minute, utc0.second);
        // 2024    10      11      16      0       21.50

        /* 2. 儒略日 */
        double dt = 600;
        double jd = utc2jd(utc0) + dt / 86400.0;
        Utc utc = jd2utc(jd);
        System.out.printf("%d\t%d\t%d\t%d\t%d\t%.2f\n", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
        // 2024    10      11      16      10      21.50

        Utc utcx = new Utc(2024, 1, 22, 10, 30, 0.);
        double jdx = utc2jd(utcx);
        System.err.println(jdx);
    }

    public static class Utc {
        int year, month, day, hour, minute;
        double second;
        Utc(int year, int month, int day, int hour, int minute, double second) {
            this.year = year;
            this.month = month;
            this.day = day;
            this.hour = hour;
            this.minute = minute;
            this.second = second;
        }
    }

    public static Utc jd2utc(double jd) {
        /* 儒略日转UTC时间 */
        double z, a, b, c, d, e, fday, alpha;
        int year, month, day, hour, minute;
        double second;

        z = Math.floor(jd + 0.5);
        fday = jd + 0.5 - z;

        if (fday < 0)
        {
            fday += 1;
            z -= 1;
        }

        if (z < 2299161)
            a = z;
        else
        {
            alpha = Math.floor((z - 1867216.25) / 36524.25);
            a = z + 1 + alpha - Math.floor(alpha / 4.0);
        }

        b = a + 1524.0;
        c = Math.floor((b - 122.1) / 365.25);
        d = Math.floor(365.25 * c);
        e = Math.floor((b - d) / 30.6001);

        day = (int)Math.floor(b - d - Math.floor(30.6001 * e) + fday);
        hour = (int)Math.floor(fday * 24);
        minute = (int)Math.floor((fday * 24 - hour) * 60);
        second = ((fday * 24 - hour) * 60 - minute) * 60;

        if (e < 14)
            month = (int) (e - 1.0);
        else
            month = (int) (e - 13.0);

        if (month > 2)
            year = (int) (c - 4716.0);
        else
            year = (int) (c - 4715.0);

        return new Utc(year, month, day, hour, minute, second);
    }

    public static double utc2jd(Utc utc) {
        /* UTC时间转儒略日 */
        double jd = 367.0 * utc.year
                    - Math.floor((7 * (utc.year + Math.floor((utc.month + 9) / 12.0))) * 0.25)
                    + Math.floor(275 * utc.month / 9) + utc.day + 1721013.5
                    + ((utc.second / 60.0 + utc.minute) / 60.0 + utc.hour) / 24.0;
        return jd;
    }
}
