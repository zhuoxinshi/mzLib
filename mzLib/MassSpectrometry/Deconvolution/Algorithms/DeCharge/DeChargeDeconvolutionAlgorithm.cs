using Chemistry;
using MassSpectrometry.Deconvolution.Parameters;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using static MassSpectrometry.IsoDecAlgorithm;

namespace MassSpectrometry
{
    internal class DeChargeDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        private DeChargeDeconvolutionParameters DeChargeParams { get; set; }
        internal DeChargeDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {
            DeChargeParams = deconParameters as DeChargeDeconvolutionParameters;
        }

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            return DeconvolutePrivateFast(spectrum,
            range, DeChargeParams.EnvelopeThreshold)
            .OrderByDescending(i => i.Score)
            .DistinctBy(i => new
            {
                i.MonoisotopicMass,
                i.Charge,
                i.Score
            })
            .ToList();
        }

        /// <summary>
        /// The function that actually does the work in deconvolution, optimized for speed. 
        /// </summary>
        /// <param name="scan"></param>
        /// <param name="deconvolutionRange"></param>
        /// <param name="spectralSimMatchThresh"></param>
        /// <returns></returns>
        internal IEnumerable<IsotopicEnvelope> DeconvolutePrivateFast(MzSpectrum scan, MzRange deconvolutionRange, double spectralSimMatchThresh)
        {
            // I think the isotopic envelope has to be a List<IsotopicEnvelope>. 
            // And charge states will get added to the isotopic envelope object. 
            // 
            ConcurrentDictionary<double, List<IsotopicEnvelope>> ieHashSet = new(new DoubleEqualityComparer());


            var output = PreFilterMzVals(scan.XArray,
                scan.YArray, DeChargeParams.MinCharge,
                DeChargeParams.MaxCharge,
                DeChargeParams.MinimumMassDa,
                DeChargeParams.MaximumMassDa,
                DeChargeParams.DeltaMass,
                DeChargeParams.PreFilterDeconvolutionType);

            var ladder = CreateChargeStateLadders(output, DeChargeParams.MinCharge, DeChargeParams.MaxCharge,
                scan.FirstX!.Value, scan.LastX!.Value);

            Parallel.ForEach(ladder, (m) =>
            {
                var index = MatchChargeStateLadder(scan, m);
                //match each Mz in m to the Mz values in scan, return a list of index

                var ladderMatch = TransformToChargeStateLadderMatch(index, scan, m, DeChargeParams.PeakMatchPpmTolerance);
                //create a ChargeStateLadderMatch with matched Mz values between scan and m, and the calculated charges
                //m is the theoretical ladder for the ChargeStateLadderMatch

                var successfulMatch = ScoreChargeStateLadderMatch(ladderMatch, scan);

                if (successfulMatch)
                {
                    FindIsotopicEnvelopes(ladderMatch!, scan, deconvolutionRange, ieHashSet, DeChargeParams.EnvelopeThreshold);
                }
            });

            return ieHashSet.Values.SelectMany(i => i);
        }

        /// <summary>
        /// Given a charge state ladder match and the original data, adds 
        /// </summary>
        /// <param name="match"></param>
        /// <param name="scan"></param>
        /// <param name="range"></param>
        /// <param name="ieHashSet"></param>
        /// <param name="minimumThreshold"></param>
        internal void FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum spectrum, MzRange range,
            ConcurrentDictionary<double, List<IsotopicEnvelope>> ieHashSet, double minimumThreshold)
        {
            //TODO: Needs to be refactored to remove side effects. 
            List<double> chargesList = match.PeakList.Select(i => i.Charge).ToList();
            List<double> mzList = match.PeakList.Select(i => i.Mz).ToList();

            //For each charge state in the scan
            for (int i = 0; i < chargesList.Count; i++)
            {
                int charge = (int)Math.Round(chargesList[i]);
                if (range.Contains(mzList[i]))
                {
                    double[] neutralXArray = spectrum.XArray.Select(j => j.ToMass(charge)).ToArray();
                    //take the experimental mz values and convert them to neutral mass? why convert all the Mz values using the same charge?
                    //??

                    MzSpectrum neutralMassSpectrum = new(neutralXArray, spectrum.YArray, true);

                    double maxMassToTake = AverageResidueModel.GetAllTheoreticalMasses(match.MassIndex).Max();
                    double minMassToTake = AverageResidueModel.GetAllTheoreticalMasses(match.MassIndex).Min();

                    MzRange newRange = new(minMassToTake, maxMassToTake);

                    var envelope = FillIsotopicEnvelopeByBounds(match, neutralMassSpectrum, newRange, charge);
                    //make an isotopic envelope using the monoisotopic mass from the ChargeStateLaddarMatch, and peaks from the neutral mass spectrum, given the range of masses

                    if (envelope != null)
                    {
                        // need to perform a scoring if the envelope only consists of low resolution data. 
                        RescoreIsotopicEnvelope(envelope);
                        if (envelope.Score >= minimumThreshold)
                        {
                            if (ieHashSet.ContainsKey(match.MonoisotopicMass))
                            {
                                if (ieHashSet.TryGetValue(match.MonoisotopicMass, out var tempList))
                                {
                                    tempList.Add(envelope);
                                }
                            }
                            else
                            {
                                ieHashSet.TryAdd(match.MonoisotopicMass, new List<IsotopicEnvelope> { envelope });
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Generates an sequence of values where the n+1 value is a the 2 * ppm tolerance away from the previous value. 
        /// </summary>
        /// <param name="minMass"></param>
        /// <param name="maxMass"></param>
        /// <param name="ppmTolerance"></param>
        /// <returns></returns>
        internal static double[] GenerateMassesIndexMann(double minMass, double maxMass, double ppmTolerance = 15)
        {
            List<double> masses = new();
            // put mass in the center of the ppm error range. 

            double value = minMass;
            double epsilon = ppmTolerance / 1e6;
            while (value < maxMass)
            {
                double newValue = value * (epsilon + 1) / (1 - epsilon);
                masses.Add(newValue);
                value = newValue;
            }
            return masses.ToArray();
        }

        internal static double[] GenerateMassesIndexMult(double minMass, double maxMass, double mzSpacing)
        {
            List<double> masses = new();
            double value = minMass;
            while (value < maxMass)
            {
                double newValue = value + mzSpacing;
                masses.Add(newValue);
                value = newValue;
            }
            return masses.ToArray();
        }

        /// <summary>
        /// Calculates theoretical charge state ladders. 
        /// </summary>
        /// <param name="indexOfMaxIntensityPeak">Index derived from the original, full scan. </param>
        /// <param name="mzValueArray"></param>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="minMzValAllowed"></param>
        /// <param name="maxMzValAllowed"></param>
        /// <returns></returns>
        internal static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(int indexOfMaxIntensityPeak,
            double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed)
        {
            double mzVal = mzValueArray[indexOfMaxIntensityPeak];
            for (int i = minCharge; i <= maxCharge; i++)
            {
                double tempMass = mzVal.ToMass(i);
                List<(int charge, double mz)> tempLadder = new();
                for (int j = maxCharge; j > minCharge; j--)
                {
                    double mz = tempMass.ToMz(j);
                    if (mz <= maxMzValAllowed && mz >= minMzValAllowed)
                    {
                        tempLadder.Add((j, tempMass.ToMz(j)));
                    }
                }
                // clean up the ladder ;
                yield return new ChargeStateLadder(tempMass, tempLadder.Select(k => k.mz).ToArray());
            }
        }

        internal static IEnumerable<double> PreFilterMzVals(double[] mzVals, double[] intensityArray, int minCharge,
            int maxCharge, double minMass, double maxMass, double spacing,
            PreFilterDeconvolutionType deconType = PreFilterDeconvolutionType.Mann)
        {
            ConcurrentDictionary<double, int> concurrentHashSet = new(new DoubleEqualityComparer());
            ConcurrentDictionary<int, double> cumulativeIntensityDict = new();



            //double medianCounts = ArrayStatistics.Mean(cumulativeIntensities.Where(i => i >= 1.1).ToArray());
            // return indexes where counts are bigger than the median. 
            HashSet<int> indexToKeepList = new();


            double[] masses = CreateMassesArray(deconType, minMass, maxMass, spacing);
            double[] cumulativeIntensityArray = new double[masses.Length];

            switch (deconType)
            {
                case PreFilterDeconvolutionType.Mann:
                    for (int i = 0; i < masses.Length; i++)
                    {
                        cumulativeIntensityDict.TryAdd(i, 0d);
                    }
                    MannDeconvolution(mzVals, intensityArray, masses, minCharge, maxCharge, minMass, maxMass,
                        concurrentHashSet, cumulativeIntensityDict, out cumulativeIntensityArray);
                    break;

                case PreFilterDeconvolutionType.Multiplicative:
                    for (int i = 0; i < masses.Length; i++)
                    {
                        cumulativeIntensityDict.TryAdd(i, 1d);
                    }
                    MultiplicativeCorrAlgoDecon(mzVals, intensityArray, masses, minCharge, maxCharge, minMass, maxMass,
                        concurrentHashSet, cumulativeIntensityDict, out cumulativeIntensityArray);
                    break;
            }


            var meanVariance = ArrayStatistics.QuantileInplace(cumulativeIntensityArray.Where(z => z > 0).ToArray(), 0.99);
            double threshold = meanVariance;


            for (int i = 0; i < cumulativeIntensityArray.Length; i++)
            {
                if (cumulativeIntensityArray[i] >= threshold)
                {
                    indexToKeepList.Add(i);
                }
            }
            // Define the threshold value as mean + 1.5 * standard deviation of the intensity values


            return concurrentHashSet
                .Where(i => indexToKeepList.Contains((i.Value)))
                .GroupBy(x => x.Value, new DoubleEqualityComparer())
                .ToDictionary(t => t.Key,
                    t => t.Select(r => r.Key).Average())
                .Select(i => i.Value);
        }

        internal static double[] CreateMassesArray(PreFilterDeconvolutionType deconType, double minMass, double maxMass, double spacing)
        {
            var masses = deconType switch
            {
                PreFilterDeconvolutionType.Mann => GenerateMassesIndexMann(minMass, maxMass, spacing),
                PreFilterDeconvolutionType.Multiplicative => GenerateMassesIndexMult(minMass, maxMass, spacing)
            };
            return masses;
        }

        internal static void MannDeconvolution(double[] mzVals, double[] intensityArray, double[] masses, int minCharge, int maxCharge,
            double minMass, double maxMass, ConcurrentDictionary<double, int> concurrentHashSet,
            ConcurrentDictionary<int, double> cumulativeIntensityDict, out double[] neutralMassIntensityArray)
        {
            Parallel.For(0, mzVals.Length, j =>
            {
                for (int i = maxCharge; i >= minCharge; i--)
                {
                    var testMass = mzVals[j].ToMass(i);
                    if (testMass > maxMass || testMass < minMass)
                    {
                        continue;
                    }

                    int index = GetBucket(masses, testMass);
                    concurrentHashSet.TryAdd(testMass, index);
                    while (true)
                    {
                        if (cumulativeIntensityDict.TryUpdate(index, cumulativeIntensityDict[index] + intensityArray[j],
                                cumulativeIntensityDict[index])) break;
                    }
                }
            });

            neutralMassIntensityArray = cumulativeIntensityDict.OrderBy(i => i.Key)
                .Select(i => i.Value).ToArray();
        }

        internal static void MultiplicativeCorrAlgoDecon(double[] mzVals, double[] intensityArray, double[] masses, int minCharge, int maxCharge,
            double minMass, double maxMass, IDictionary<double, int> concurrentHashSet,
            IDictionary<int, double> cumulativeIntensityDict, out double[] neutralMassIntensityArray)
        {
            double scanRmsIntensity = intensityArray.RootMeanSquare();

            int[] massesLengthTrackingArray = new int[masses.Length];

            for (int j = 0; j < mzVals.Length; j++)
            {
                for (int i = maxCharge; i >= minCharge; i--)
                {
                    var testMass = mzVals[j].ToMass(i);
                    if (testMass > maxMass || testMass < minMass)
                    {
                        continue;
                    }

                    int index = GetBucket(masses, testMass);
                    concurrentHashSet.TryAdd(testMass, index);

                    cumulativeIntensityDict[index] *= intensityArray[j];
                    massesLengthTrackingArray[index]++;
                }
            }

            // apply the rms correction, which is sum of intensities / rms^n, where n is the number of points in the calculation
            var orderedIntensityDict = cumulativeIntensityDict
                .OrderBy(i => i.Key)
                .ToDictionary(i => i.Key, i => i.Value);
            for (int i = 0; i < massesLengthTrackingArray.Length; i++)
            {
                orderedIntensityDict[i] /= Math.Pow(scanRmsIntensity, massesLengthTrackingArray[i]);
            }

            neutralMassIntensityArray = orderedIntensityDict.Values.Select(i => i - 1d).ToArray();
        }

        internal static int GetBucket(double[] array, double value)
        {
            int index = Array.BinarySearch(array, value);
            if (index < 0)
            {
                index = ~index;
            }

            if (index >= array.Length)
            {
                index = array.Length - 1;
            }

            if (index != 0 && array[index] - value >
                value - array[index - 1])
            {
                index--;
            }
            return index;
        }
        internal static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(IEnumerable<double> neutralMass,
            int minCharge, int maxCharge, double minMzAllowed, double maxMzAllowed)
        {
            foreach (var mass in neutralMass)
            {
                List<(int charge, double mz)> tempLadder = new();
                for (int j = maxCharge; j > minCharge; j--)
                {
                    double mz = mass.ToMz(j);
                    if (mz <= maxMzAllowed && mz >= minMzAllowed)
                    {
                        tempLadder.Add((j, mass.ToMz(j)));
                    }
                }
                yield return new ChargeStateLadder(mass, tempLadder.Select(k => k.mz).ToArray());
            }
        }

        internal static List<int> MatchChargeStateLadder(MzSpectrum scan, ChargeStateLadder ladder)
        {
            ConcurrentBag<int> output = new();
            foreach (var t in ladder.MzVals)
            {
                int index = GetBucket(scan.XArray, t);
                output.Add(index);
            }

            return output.ToList();
        }

        internal static ChargeStateLadderMatch TransformToChargeStateLadderMatch(List<int> ladderToIndexMap,
            MzSpectrum scan, ChargeStateLadder ladder, double ppmMatchTolerance)
        {
            List<double> listMzVals = new();
            List<double> listIntVals = new();
            for (int i = 0; i < ladderToIndexMap.Count; i++)
            {
                listMzVals.Add(scan.XArray[ladderToIndexMap[i]]);
                listIntVals.Add(scan.YArray[ladderToIndexMap[i]]);
            }
            List<double> chargesOfMatchingPeaks = listMzVals.Select(i => ladder.Mass / i).ToList();
            //should the neutral mass (ladder.Mass) be used here to calculate charge?

            return new ChargeStateLadderMatch(listMzVals, listIntVals,
                chargesOfMatchingPeaks, ladder, ppmMatchTolerance);
        }

        /// <summary>
        /// Performs spectral similarity calculation between the theoretical isotopic envelopes and the experimental data.
        /// </summary>
        /// <param name="envelope"></param>
        internal void RescoreIsotopicEnvelope(IsotopicEnvelope envelope)
        {
            if (envelope == null) return;

            int massIndex = AverageResidueModel.GetMostIntenseMassIndex(envelope.MostAbundantObservedIsotopicMass);
            double[] theoreticalMasses = AverageResidueModel.GetAllTheoreticalMasses(massIndex);
            double[] theoreticalIntensities = AverageResidueModel.GetAllTheoreticalIntensities(massIndex);
            double diff = envelope.MonoisotopicMass + Averagine.DiffToMonoisotopic[massIndex] - theoreticalMasses[0];
            double[] theoreticalMzs = AverageResidueModel.GetAllTheoreticalMasses(massIndex).Select(i => (i + diff).ToMz(envelope.Charge)).ToArray();
            double[] normalizedTheoreticalIntensities = theoreticalIntensities.Select(i => i / theoreticalIntensities.Max()).ToArray();

            var theoreticalSpectrum = new MzSpectrum(theoreticalMzs, normalizedTheoreticalIntensities, true)
                .FilterByY(0.001, 1.0);
            var spectrum0 = new MzSpectrum(theoreticalSpectrum.Select(i => (i.Mz)).ToArray(),
                theoreticalSpectrum.Select(i => i.Intensity).ToArray(), true);

            MzSpectrum experimentalSpectrum = new MzSpectrum(envelope.Peaks.Select(i => i.mz.ToMz(envelope.Charge)).ToArray(),
                envelope.Peaks.Select(i => i.intensity).ToArray(), true);

            SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, spectrum0,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, DeChargeParams.PeakMatchPpmTolerance,
                true);

            double? score = similarity.SpectralContrastAngle();
            if (score.HasValue)
            {
                envelope.Rescore(score.Value);
            }
            else
            {
                envelope.Rescore(0);
            }
        }

        internal bool ScoreChargeStateLadderMatch(ChargeStateLadderMatch match, MzSpectrum scan)
        {
            GetMassIndex(match);

            match.CompareTheoreticalNumberChargeStatesVsActual();
            if (match.PercentageMzValsMatched < DeChargeParams.PercentageMatchedThresh) return false;

            match.CalculateChargeStateScore();
            if (match.SequentialChargeStateScore < DeChargeParams.SequentialChargeStateDiff) return false;

            match.CalculateEnvelopeScore();
            if (match.EnvelopeScore < DeChargeParams.EnvelopeScoreThresh) return false;

            CalculateMonoMass(match);

            return true;
        }

        internal void CalculateMonoMass(ChargeStateLadderMatch match)
        {
            var diff = AverageResidueModel.GetDiffToMonoisotopic(match.MassIndex);
            match.MonoGuesses = match.PeakList
                .Select(i => i.NeutralMass - diff)
                .ToList();
        }

        /// <summary>
        /// Adds the mass index to the ChargeStateLadderMatch object. Must be run from the deconvolution to have access to the precalculated averagine isotopic envelopes. 
        /// </summary>
        /// <param name="match"></param>
        /// <param name="scan"></param>
        public void GetMassIndex(ChargeStateLadderMatch match)
        {
            int massIndex = AverageResidueModel.GetMostIntenseMassIndex(match.TheoreticalLadder.Mass);
            match.MassIndex = massIndex;
        }

        /// <summary>
        /// get the range of the theoretical isotopic envelopes and then pulls all the peaks within that range from the original data. 
        /// </summary>
        /// <param name="match"></param>
        /// <param name="scan"></param>
        /// <param name="isolationRange"></param>
        /// <param name="chargeState"></param>
        /// <returns></returns>
        public IsotopicEnvelope FillIsotopicEnvelopeByBounds(ChargeStateLadderMatch match, MzSpectrum scan, MzRange isolationRange, int chargeState)
        {
            List<(double, double)> listOfPeaks = scan.Extract(isolationRange).Select(i => (i.Mz, i.Intensity)).ToList();
            double totalIntensity = listOfPeaks.Sum(i => i.Item2);
            if (listOfPeaks.Any())
            {
                return new IsotopicEnvelope(listOfPeaks, match.MonoisotopicMass, chargeState,
                                       totalIntensity, match.StdDev);
            }

            return null;
        }
    }

    public class DoubleEqualityComparer : IEqualityComparer<double>
    {
        public bool Equals(double a, double b)
        {
            return Math.Round(a, 2) == Math.Round(b, 2);
        }

        public int GetHashCode(double value)
        {
            return Math.Round(value, 2).GetHashCode();
        }
    }
}
