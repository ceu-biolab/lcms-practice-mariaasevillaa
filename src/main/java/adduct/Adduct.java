package adduct;
import java.util.regex.*;
public class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {
        int multimer = extractMultimer(adduct); //(\\d+)M
        int charge = extractCharge(adduct);// (\\d+)([+-]$)
        Double massToSearch;
        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (adductMass == null) {
            adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if (adductMass == null) {
            throw new IllegalArgumentException("Invalid adduct: " + adduct);
        }
        massToSearch = (mz * charge) + adductMass / multimer;
        return massToSearch;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).

        /*
        if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge

        if adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer

        return monoisotopicMass;

         */

    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        Double mz;

        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (adductMass == null) {
            adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if(adductMass == null) {
            throw new IllegalArgumentException("Invalid adduct: " + adduct);
        }
        mz = ((monoisotopicMass * multimer) - adductMass) / charge;

        return mz;
    }
    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass   Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }
    private static int extractMultimer (String adduct){
        Pattern pMultimer = Pattern.compile("(\\d+)M"); //"([0-9]+)M")
        Matcher mMultimer = pMultimer.matcher(adduct);

        if (mMultimer.find()) {
            return Integer.parseInt(mMultimer.group(1));
        } else {
            return 1;
        }
    }

    private static int extractCharge (String adduct){
        Pattern pCharge = Pattern.compile("(\\d+)([+-]$)"); //"([+-]?\\d+)"
        Matcher mCharge = pCharge.matcher(adduct);

        if (mCharge.find()) {
            return Integer.parseInt(mCharge.group(1));
        } else {
            return 1;
        }
    }
}





