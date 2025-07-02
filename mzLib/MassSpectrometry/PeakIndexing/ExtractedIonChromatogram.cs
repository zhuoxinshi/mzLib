using MathNet.Numerics.Interpolation;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    /// <summary>
    /// A generic XIC class for all IIndexedPeak objects (mz peak, isotopic envelope, etc.) that can be traced across retention time.
    /// </summary>
    public class ExtractedIonChromatogram
    {
        public List<IIndexedPeak> Peaks { get; set; }

        public double ApexRT;
        public int ApexScanIndex;
        public double AveragedM => AverageM();
        public (double, double)[] XYData { get; set; }
        public double[] NormalizedPeakIntensities { get; set; }
        public double StartRT { get; set; }
        public double EndRT { get; set; }
        public int StartScanIndex { get; set; }
        public int EndScanIndex { get; set; }
        public double AverageM()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            double averagedM = 0;
            foreach (var peak in Peaks)
            {
                double weight = peak.Intensity / sumIntensity;
                averagedM += weight * peak.M;
            }
            return averagedM;
        }

        public ExtractedIonChromatogram(List<IIndexedPeak> peaks)
        {
            Peaks = peaks;
            ApexRT = Peaks.OrderByDescending(p => p.Intensity).First().RetentionTime;
            StartRT = Peaks.Min(p => p.RetentionTime);
            EndRT = Peaks.Max(p => p.RetentionTime);
            ApexScanIndex = Peaks.OrderByDescending(p => p.Intensity).First().ZeroBasedScanIndex;
            StartScanIndex = Peaks.Min(p => p.ZeroBasedScanIndex);
            EndScanIndex = Peaks.Max(p => p.ZeroBasedScanIndex);
        }

        public void SetNormalizedPeakIntensities()
        {
            double sumIntensity = Peaks.Sum(p => p.Intensity);
            NormalizedPeakIntensities = Peaks.Select(p => p.Intensity / sumIntensity * 100).ToArray();
        }

        public static List<ExtractedIonChromatogram> GetAllXICsFromIndexedPeaks(IndexingEngine<IIndexedPeak> indexingEngine, Tolerance peakFindingTolerance, int maxMissedScanAllowed, double maxRTRange, int numPeakThreshold = 2)
        {
            var xics = new List<ExtractedIonChromatogram>();
            var sortedPeaks = indexingEngine.IndexedPeaks.SelectMany(peaks => peaks).OrderBy(p => p.Intensity).ToList();
            foreach(var peak in sortedPeaks)
            {
                if (peak.XIC == null)
                {
                    var peakList = indexingEngine.GetXic(peak.M, peak.RetentionTime, peakFindingTolerance, maxMissedScanAllowed, maxRTRange);
                    if (peakList.Count >= numPeakThreshold)
                    {
                        var newXIC = new ExtractedIonChromatogram(peakList);
                        foreach(var p in peakList)
                        {
                            p.XIC = newXIC; 
                        }
                        xics.Add(newXIC);
                    }
                }
            }
            return xics;
        }
    }
}
