using MathNet.Numerics.Statistics;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Resources;
using Chemistry;
using MassSpectrometry.Deconvolution.Parameters;

namespace MassSpectrometry
{
    public class ChargeStateLadderMatch
    {
        public ChargeStateLadder TheoreticalLadder { get; set; }
        public List<MatchedPeak> PeakList { get; set; }
        public double EnvelopeScore { get; set; }
        public double Score { get; set; }
        public int MassIndex { get; set; }
        public double MonoisotopicMass => MonoGuesses.Average();
        public double StdDev => MonoGuesses.StandardDeviation();
        public double SequentialChargeStateScore { get; set; }
        public List<double> MonoGuesses { get; set; }
        public List<double> MonoErrors { get; set; }
        public double PercentageMzValsMatched { get; set; }

        public ChargeStateLadderMatch(List<double> matchingMzPeaks,
            List<double> intensitiesOfMatching, List<double> chargesOfMatchingPeaks,
            ChargeStateLadder theoreticalLadder, double ppmMatchErrorTolerance)
        {
            MonoGuesses = new List<double>();
            MonoErrors = new List<double>();
            TheoreticalLadder = theoreticalLadder;

            List<MatchedPeak> peakList = new();
            for (int i = 0; i < matchingMzPeaks.Count; i++)
            {
                var peak = new MatchedPeak(matchingMzPeaks[i], intensitiesOfMatching[i],
                    chargesOfMatchingPeaks[i]);
                peak.CalculateMatchError(TheoreticalLadder);
                if (Math.Abs(peak.MatchError) <= ppmMatchErrorTolerance)
                {
                    peakList.Add(peak);
                }
            }
            PeakList = peakList;
        }

        /// <summary>
        /// Calculates the average spacing between distinct charge states. If the charge states are sequential,
        /// then the charges will be close to -1. If the charge states correspond to a high harmonic, then the values of the charges states are going to be
        /// less than -1.
        /// </summary>
        internal void CalculateChargeStateScore()
        {

            var chargesList = PeakList
                .Select(i => i.Charge)
                .Zip(PeakList.Select(i => i.Charge).Skip(1), (x, y) => y - x)
                // Where clause removes any diffs that are less than 0.1, which would indicate that the 
                // peak matching tolerance was broad enough to grab more than one peak. However, because the 
                // ppm tolerance for peak matching is so low, the difference between consecutive peaks should remain 
                // approximately the same. 
                .Where(i => Math.Abs(i) > 0.1)
                .ToList();
            if (chargesList.Count() > 1)
            {
                var numberUniqueChargesWithCorrectChargeState = chargesList
                    .Select(i => Math.Abs(i - 1d))
                    .Count(i => i <= 0.1);

                SequentialChargeStateScore = (double)numberUniqueChargesWithCorrectChargeState / (double)chargesList.Count;
                return;
            }

            SequentialChargeStateScore = -10000;
        }
        /// <summary>
        /// The envelope score determines the coefficient of best fit between the (mz,intensity) values of the
        /// matched peaks to a second order polynomial. This will eliminate peaks that are primarily made of noise values and
        /// peaks that don't have enough points to calculate the score. 
        /// </summary>
        internal void CalculateEnvelopeScore()
        {
            if (PeakList.Select(i => i.Mz).Count() < 3)
            {
                EnvelopeScore = 0d;
                return;
            }
            double[] chargesArray = PeakList.Select(i => i.Charge).ToArray();
            double[] intensitiesArray = PeakList.Select(i => i.Intensity).ToArray();

            double[] coefficients =
                Fit.Polynomial(chargesArray, intensitiesArray, 2);
            double c = coefficients[0];
            double b = coefficients[1];
            double a = coefficients[2];

            // calculate theoretical polynomial to get R^2. 
            double[] theoreticalPolynom = chargesArray
                .Select(x => c + b * x + a * x * x)
                .ToArray();
            double sum1 = 0;
            double sum2 = 0;
            double yMean = intensitiesArray.Mean();

            for (int i = 0; i < chargesArray.Length; i++)
            {
                sum1 += Math.Pow(intensitiesArray[i] - theoreticalPolynom[i], 2);
                sum2 += Math.Pow(intensitiesArray[i] - yMean, 2);
            }

            EnvelopeScore = 1 - sum1 / sum2;
        }
        /// <summary>
        /// Calculates the percentage of matched mz values in an experimental spectrum compared to a set of theoretical mz values. 
        /// </summary>
        /// <param name="ladderMatch"></param>
        /// <returns></returns>
        internal void CompareTheoreticalNumberChargeStatesVsActual()
        {
            int integerUniqueChargeValuesLength = PeakList.Select(i => i.Charge)
                .Select(i => (int)Math.Round(i))
                .Distinct()
                .Count();

            PercentageMzValsMatched = (double)integerUniqueChargeValuesLength / (double)TheoreticalLadder.MzVals.Length;
        }
    }

    public class MatchedPeak
    {
        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double Charge { get; set; }
        public double MatchError { get; set; }
        public double NeutralMass => CalculateNeutralMass();
        public MatchedPeak(double mz, double intensity, double charge)
        {
            Mz = mz;
            Intensity = intensity;
            Charge = charge;
        }

        public void CalculateMatchError(ChargeStateLadder theoreticalLadder)
        {
            int index = DeChargeDeconvolutionAlgorithm.GetBucket(theoreticalLadder.MzVals, Mz);
            MatchError = Math.Abs((Mz - theoreticalLadder.MzVals[index]) / Mz * 1E6);
        }

        private double CalculateNeutralMass()
        {
            return Mz.ToMass((int)Math.Round(Charge));
        }
    }
}

