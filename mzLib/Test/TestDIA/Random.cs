using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.TestDIA
{
    public class Random
    {
        [Test]
        public static void TestFlashLFQSetting()
        {
            var lcmsPeak = new ChromatographicPeak(null, false, null);
        }

        [Test]
        public static void TestDecon()
        {
            string dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            var dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            var ms1scans = dataFile.GetMsDataScans().Where(s => s.MsnOrder == 1).ToList();
            DeconvolutionParameters parameters = new ClassicDeconvolutionParameters(1, 5, 4, 3); ;
            var envelopes = Deconvoluter.Deconvolute(ms1scans[1], parameters);
;        }
    }
}
