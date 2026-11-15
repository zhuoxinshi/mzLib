using ISD;
using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.IsdTests
{
    /// <summary>
    /// Regenerates the mzML files behind the HCD_YBC ISD analysis from the instrument .raw files.
    /// Explicit: they need the local .raw files and write about 15 GB.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestISD
    {
        // the ten HCD_YBC ISD60-80-100 .raw files
        private const string RawFolder = @"E:\ISD Project\071626_ISD\Quant\HCD_YBC";
        private const string OutputFolder = @"E:\Proteomics_software\TopPIC\toppic-windows-1.8.1\Data\ISD_quant\HCD_YBC\ISD\Reproduce_test";

        private static List<string> IsdRawFiles()
        {
            return Directory.GetFiles(RawFolder, "*.raw").Where(f => f.Contains("ISD60-80-100")).ToList();
        }

        /// <summary>
        /// _labelCorrected.mzML, the TopFD input for TopFD_labelCorrected_redo
        /// (originally test ISD0725, which wrote next to the .raw files).
        /// </summary>
        [Test]
        [Explicit("needs the local HCD_YBC .raw files")]
        public static void WriteLabelCorrectedMzml()
        {
            var outFolder = Path.Combine(OutputFolder, "labelCorrected");
            Directory.CreateDirectory(outFolder);
            foreach (var file in IsdRawFiles())
            {
                LabelCorrection.ISD60_80_100_LabelCorrectionFromRaw(file, outFolder);
            }
        }

        /// <summary>
        /// _ms1/_isd60/_isd80/_isd100.mzML, the TopFD inputs for SeparateISD\TopFD
        /// (originally test SeparateFilesForTopFD).
        /// </summary>
        [Test]
        [Explicit("needs the local HCD_YBC .raw files")]
        public static void WriteSidSeparatedMzml()
        {
            var outFolder = Path.Combine(OutputFolder, "SeparateISD");
            Directory.CreateDirectory(outFolder);
            foreach (var file in IsdRawFiles())
            {
                LabelCorrection.SeparateFilesForTopFD(file, outFolder);
            }
        }
    }
}
