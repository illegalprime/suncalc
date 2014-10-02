package edu.gatech.meden3.arealdaysnight.util;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

/**
 * Created by mde on 10/1/14.
 */
public class CelestialCalculations {

    // sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas
    // date/time constants and conversions
    private static final int dayMs = 1000 * 60 * 60 * 24;
    private static final int J1970 = 2440588;
    private static final int J2000 = 2451545;
    private static final double J0 = 0.0009;
    // obliquity of the Earth
    private static final double rad = Math.PI / 180;
    private static final double  e  = rad * 23.4397;

    private static final int DECLINATION     = 0;
    private static final int RIGHT_ASCENSION = 1;
    private static final int DISTANCE        = 2;

    private static final int AZIMUTH  = 0;
    private static final int ALTITUDE = 1;

    private static final int FRACTION = 0;
    private static final int PHASE    = 1;
    private static final int ANGLE    = 2;

    private static final int SOLAR_NOON = 0;
    private static final int NADIR      = 1;

    private static final int sunrise = 0,       sunset = 1;
    private static final int sunriseEnd = 0,    sunsetStart = 1;
    private static final int dawn = 0,          dusk = 1;
    private static final int nauticalDawn = 0,  nauticalDusk = 1;
    private static final int nightEnd = 0,      night = 1;
    private static final int goldenHourEnd = 0, goldenHour = 1;

    private static final ArrayList<String[]> timesString = new ArrayList<String[]>();
    private static final ArrayList<Double>   timesDouble = new ArrayList<Double>();


    public static double toJulian(Date date) {
        return date.getTime() / dayMs - 0.5 + J1970;
    }

    private static Date fromJulian(double j)  {
        return new Date((long) (j + 0.5 - J1970) * dayMs);
    }

    private static double toDays(Date date) {
        return toJulian(date) - J2000;
    }

    // general calculations for position
    private static double rightAscension(double l, double b) {
        return Math.atan((Math.sin(l) * Math.cos(e) - Math.tan(b) * Math.sin(e)) / Math.cos(l));
    }

    private static double declination(double l, double b) {
        return Math.asin(Math.sin(b) * Math.cos(e) + Math.cos(b) * Math.sin(e) * Math.sin(l));
    }

    private static double azimuth(double H, double phi, double dec) {
        return Math.atan(Math.sin(H) / (Math.cos(H) * Math.sin(phi) - Math.tan(dec) * Math.cos(phi)));
    }

    private static double altitude(double H, double phi, double dec) {
        return Math.asin(Math.sin(phi) * Math.sin(dec) + Math.cos(phi) * Math.cos(dec) * Math.cos(H));
    }

    private static double siderealTime(double d, double lw) {
        return rad * (280.16 + 360.9856235 * d) - lw;
    }

    // general sun calculations
    private static double solarMeanAnomaly(double d) {
        return rad * (357.5291 + 0.98560028 * d);
    }

    private static double eclipticLongitude(double M) {
        double C = rad * (1.9148 * Math.sin(M) + 0.02 * Math.sin(2 * M) + 0.0003 * Math.sin(3 * M)), // equation of center
        P = rad * 102.9372; // perihelion of the Earth

        return M + C + P + Math.PI;
    }

    private static double[] sunCoords(double d) {
        double M = solarMeanAnomaly(d);
        double L = eclipticLongitude(M);

        double[] value = new double[2];
        value[DECLINATION]     = new Double(declination(L, 0));
        value[RIGHT_ASCENSION] = new Double(rightAscension(L, 0));
        return value;
    }


    // calculates sun position for a given date and latitude/longitude
    private static double[] getPosition(Date date, double lat, double lng) {
        double   lw  = rad * -lng;
        double   phi = rad * lat;
        double   d   = toDays(date);
        double[] c   = sunCoords(d);
        double   H   = siderealTime(d, lw) - c[RIGHT_ASCENSION];
        double[] out = new double[2];

        out[AZIMUTH]  = azimuth(H, phi, c[DECLINATION]);
        out[ALTITUDE] = altitude(H, phi, c[DECLINATION]);
        return out;
    }

    private static void init() {
        // sun times configuration (angle, morning name, evening name)
        String[] data1 = {"sunrise", "sunset"};
        String[] data2 = {"sunriseEnd", "sunsetStart"};
        String[] data3 = {"dawn", "dusk"};
        String[] data4 = {"nauticalDawn", "nauticalDusk"};
        String[] data5 = {"nightEnd", "night"};
        String[] data6 = {"goldenHourEnd", "goldenHour"};
        timesString.add(data1);
        timesString.add(data2);
        timesString.add(data3);
        timesString.add(data4);
        timesString.add(data5);
        timesString.add(data6);

        timesDouble.add(-0.833);
        timesDouble.add(-0.3);
        timesDouble.add(-6.0);
        timesDouble.add(-12.0);
        timesDouble.add(-18.0);
        timesDouble.add(6.0);
    }

    // adds a custom time to the times config
    public static void addTime (double angle, String riseName, String setName) {
        String[] data = {riseName, setName};
        timesString.add(data);
        timesDouble.add(angle);
    }


    // calculations for sun times
    private static double julianCycle(double d, double lw) {
        return Math.round(d - J0 - lw / (2 * Math.PI));
    }

    private static double approxTransit(double Ht, double lw, double n) {
        return J0 + (Ht + lw) / (2 * Math.PI) + n;
    }

    private static double solarTransitJ(double ds, double M, double L)  {
        return J2000 + ds + 0.0053 * Math.sin(M) - 0.0069 * Math.sin(2 * L);
    }

    private static double hourAngle(double h, double phi, double d) {
        return Math.acos((Math.sin(h) - Math.sin(phi) * Math.sin(d)) / (Math.cos(phi) * Math.cos(d)));
    }

    // returns set time for the given sun altitude
    private static double getSetJ(double h, double lw, double phi, double dec, double n, double M, double L) {
        double w = hourAngle(h, phi, dec);
        double a = approxTransit(w, lw, n);
        return solarTransitJ(a, M, L);
    }


    // calculates sun times for a given date and latitude/longitude
    private static HashMap<String, Object> getTimes(Date date, double lat, double lng) {

        double lw = rad * -lng;
        double phi = rad * lat;

        double d = toDays(date);
        double n = julianCycle(d, lw);
        double ds = approxTransit(0, lw, n);

        double M = solarMeanAnomaly(ds);
        double L = eclipticLongitude(M);
        double dec = declination(L, 0);

        double Jnoon = solarTransitJ(ds, M, L);
        double Jset, Jrise;

        HashMap<String, Object> result = new HashMap<String, Object>();
        result.put("solarNoon", fromJulian(Jnoon));
        result.put("solarNoon", fromJulian(Jnoon - 0.5));

        for (int strs = 0, angls = 0, lenStrs = timesString.size(), lenAngls = timesDouble.size();
             strs < lenStrs && angls < lenAngls; ++strs, ++angls) {

            Jset = getSetJ(timesDouble.get(angls) * rad, lw, phi, dec, n, M, L);
            Jrise = Jnoon - (Jset - Jnoon);

            result.put(timesString.get(strs)[0], fromJulian(Jrise));
            result.put(timesString.get(strs)[1], fromJulian(Jset));
        }

        return result;
    }


    // moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas
    private static double[] moonCoords(double d) { // geocentric ecliptic coordinates of the moon

        double L = rad * (218.316 + 13.176396 * d); // ecliptic longitude
        double M = rad * (134.963 + 13.064993 * d); // mean anomaly
        double F = rad * (93.272 + 13.229350 * d);  // mean distance
        double l  = L + rad * 6.289 * Math.sin(M); // longitude
        double b  = rad * 5.128 * Math.sin(F);     // latitude
        double dt = 385001 - 20905 * Math.cos(M);  // distance to the moon in km

        double[] data = new double[3];
        data[RIGHT_ASCENSION] = rightAscension(l, b);
        data[DECLINATION]     = declination(l, b);
        data[DISTANCE]        =  dt;
        return data;
    }

    private static double[] getMoonPosition(Date date, double lat, double lng) {

        double lw  = rad * -lng;
        double phi = rad * lat;
        double d   = toDays(date);

        double[] c = moonCoords(d);
        double H = siderealTime(d, lw) - c[RIGHT_ASCENSION];
        double h = altitude(H, phi, c[DECLINATION]);

        // altitude correction for refraction
        h = h + rad * 0.017 / Math.tan(h + rad * 10.26 / (h + rad * 5.10));

        double[] data = new double[3];
        data[AZIMUTH] = azimuth(H, phi, c[DECLINATION]);
        data[ALTITUDE] = h;
        data[DISTANCE] = c[DISTANCE];
        return data;
    }

    // calculations for illumination parameters of the moon,
    // based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and
    // Chapter 48 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
    private static double[] getMoonIllumination(Date date) {
        double d = toDays(date);
        double[] s = sunCoords(d);
        double[] m = moonCoords(d);

        long sdist = 149598000; // distance from Earth to Sun in km

        double    phi = Math.acos(Math.sin(s[DECLINATION]) * Math.sin(m[DECLINATION]) + Math.cos(s[DECLINATION]) * Math.cos(m[DECLINATION]) * Math.cos(s[RIGHT_ASCENSION] - m[RIGHT_ASCENSION]));
        double    inc = Math.atan((sdist * Math.sin(phi)) / (m[DISTANCE] - sdist * Math.cos(phi)));
        double    angle = Math.atan((Math.cos(s[DECLINATION]) * Math.sin(s[RIGHT_ASCENSION] - m[RIGHT_ASCENSION])) / (Math.sin(s[DECLINATION]) * Math.cos(m[DECLINATION]) - Math.cos(s[DECLINATION]) * Math.sin(m[DECLINATION]) * Math.cos(s[RIGHT_ASCENSION] - m[RIGHT_ASCENSION])));


        double[] data = new double[3];
        data[FRACTION] = (1 + Math.cos(inc)) / 2;
        data[PHASE] = 0.5 + 0.5 * inc * (angle < 0 ? -1 : 1) / Math.PI;
        data[ANGLE] = angle;
        return data;
    }
}
