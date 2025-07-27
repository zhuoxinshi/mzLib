using MathNet.Numerics.RootFinding;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.PeakIndexing
{
    public class ExtractedMassChromatogram : ExtractedIonChromatogram
    {
        public double MonoIsotopicMass { get; set; } 
        public int Charge { get; set; }

        public ExtractedMassChromatogram(List<IIndexedPeak> peaks)
            : base(peaks)
        {
            MonoIsotopicMass = peaks.MaxBy(p => p.Intensity).M;
        }
    }
}
