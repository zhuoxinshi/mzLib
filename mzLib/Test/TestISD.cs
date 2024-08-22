using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    internal class TestISD
    {
        [Test]
        public static void SnipMzMlForISD()
        {
            string origDataFile = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample1_2uL_ISD.mzML";
            FilteringParams filter = new FilteringParams(200, 0.01, 1, null, false, false, true);
            var reader = MsDataFileReader.GetDataFile(origDataFile);
            reader.LoadAllStaticData(filter, 1);

            var scans = reader.GetAllScansList();
            int startScan = 1449;
            int endScan = 1552;
            var scansToKeep = scans.Where(x => x.OneBasedScanNumber >= startScan && x.OneBasedScanNumber <= endScan).ToList();

            List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)> scanNumbers = new List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)>();

            foreach (var scan in scansToKeep)
            {
                if (scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    scanNumbers.Add((scan.OneBasedScanNumber, scan.OneBasedPrecursorScanNumber.Value));
                }
                else
                {
                    scanNumbers.Add((scan.OneBasedScanNumber, null));
                }
            }

            Dictionary<int, int> scanNumberMap = new Dictionary<int, int>();

            foreach (var scanNumber in scanNumbers)
            {
                if (!scanNumberMap.ContainsKey(scanNumber.oneBasedScanNumber) && (scanNumber.oneBasedScanNumber - startScan + 1) > 0)
                {
                    scanNumberMap.Add(scanNumber.oneBasedScanNumber, scanNumber.oneBasedScanNumber - startScan + 1);
                }
                if (scanNumber.oneBasedPrecursorScanNumber.HasValue && !scanNumberMap.ContainsKey(scanNumber.oneBasedPrecursorScanNumber.Value) && (scanNumber.oneBasedPrecursorScanNumber.Value - startScan + 1) > 0)
                {
                    scanNumberMap.Add(scanNumber.oneBasedPrecursorScanNumber.Value, scanNumber.oneBasedPrecursorScanNumber.Value - startScan + 1);
                }
            }
            List<MsDataScan> scansForTheNewFile = new List<MsDataScan>();


            foreach (var scanNumber in scanNumbers)
            {
                MsDataScan scan = scansToKeep.First(x => x.OneBasedScanNumber == scanNumber.oneBasedScanNumber);

                int? newOneBasedPrecursorScanNumber = null;
                if (scan.OneBasedPrecursorScanNumber.HasValue && scanNumberMap.ContainsKey(scan.OneBasedPrecursorScanNumber.Value))
                {
                    newOneBasedPrecursorScanNumber = scanNumberMap[scan.OneBasedPrecursorScanNumber.Value];
                }
                MsDataScan newDataScan = new MsDataScan(
                    scan.MassSpectrum,
                    scanNumberMap[scan.OneBasedScanNumber],
                    scan.MsnOrder,
                    scan.IsCentroid,
                    scan.Polarity,
                    scan.RetentionTime,
                    scan.ScanWindowRange,
                    scan.ScanFilter,
                    scan.MzAnalyzer,
                    scan.TotalIonCurrent,
                    scan.InjectionTime,
                    scan.NoiseData,
                    scan.NativeId.Replace(scan.OneBasedScanNumber.ToString(), scanNumberMap[scan.OneBasedScanNumber].ToString()),
                    scan.SelectedIonMZ,
                    scan.SelectedIonChargeStateGuess,
                    scan.SelectedIonIntensity,
                    scan.IsolationMz,
                    scan.IsolationWidth,
                    scan.DissociationType,
                    newOneBasedPrecursorScanNumber,
                    scan.SelectedIonMonoisotopicGuessMz,
                    scan.HcdEnergy
                );
                scansForTheNewFile.Add(newDataScan);
            }
            string replace = "_RT" + Math.Round(scans.Where(s => s.OneBasedScanNumber == startScan).First().RetentionTime, 2) + "-" +
                Math.Round(scans.Where(s => s.OneBasedScanNumber == endScan).First().RetentionTime, 2) + ".mzML";
            string outPath = origDataFile.Replace(".mzML", replace).ToString();

            SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
                reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.FileChecksumType, reader.SourceFile.Uri, reader.SourceFile.Id, reader.SourceFile.FileName);


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scansForTheNewFile.ToArray(), sourceFile), outPath, false);

        }

        [Test]
        public static void SnipMzMlForDDA()
        {
            string origDataFile = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_DDA.raw";
            FilteringParams filter = new FilteringParams(200, 0.01, 1, null, false, false, true);
            var reader = MsDataFileReader.GetDataFile(origDataFile);
            reader.LoadAllStaticData(filter, 1);

            var scans = reader.GetAllScansList();
            int startScan = 1449;
            int endScan = 2156;
            var scansToKeep = scans.Where(x => x.OneBasedScanNumber >= startScan && x.OneBasedScanNumber <= endScan).ToList();

            List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)> scanNumbers = new List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)>();

            foreach (var scan in scansToKeep)
            {
                if (scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    scanNumbers.Add((scan.OneBasedScanNumber, scan.OneBasedPrecursorScanNumber.Value));
                }
                else
                {
                    scanNumbers.Add((scan.OneBasedScanNumber, null));
                }
            }

            Dictionary<int, int> scanNumberMap = new Dictionary<int, int>();

            foreach (var scanNumber in scanNumbers)
            {
                if (!scanNumberMap.ContainsKey(scanNumber.oneBasedScanNumber) && (scanNumber.oneBasedScanNumber - startScan + 1) > 0)
                {
                    scanNumberMap.Add(scanNumber.oneBasedScanNumber, scanNumber.oneBasedScanNumber - startScan + 1);
                }
                if (scanNumber.oneBasedPrecursorScanNumber.HasValue && !scanNumberMap.ContainsKey(scanNumber.oneBasedPrecursorScanNumber.Value) && (scanNumber.oneBasedPrecursorScanNumber.Value - startScan + 1) > 0)
                {
                    scanNumberMap.Add(scanNumber.oneBasedPrecursorScanNumber.Value, scanNumber.oneBasedPrecursorScanNumber.Value - startScan + 1);
                }
            }
            List<MsDataScan> scansForTheNewFile = new List<MsDataScan>();


            foreach (var scanNumber in scanNumbers)
            {
                MsDataScan scan = scansToKeep.First(x => x.OneBasedScanNumber == scanNumber.oneBasedScanNumber);

                int? newOneBasedPrecursorScanNumber = null;
                if (scan.OneBasedPrecursorScanNumber.HasValue && scanNumberMap.ContainsKey(scan.OneBasedPrecursorScanNumber.Value))
                {
                    newOneBasedPrecursorScanNumber = scanNumberMap[scan.OneBasedPrecursorScanNumber.Value];
                }
                MsDataScan newDataScan = new MsDataScan(
                    scan.MassSpectrum,
                    scanNumberMap[scan.OneBasedScanNumber],
                    scan.MsnOrder,
                    scan.IsCentroid,
                    scan.Polarity,
                    scan.RetentionTime,
                    scan.ScanWindowRange,
                    scan.ScanFilter,
                    scan.MzAnalyzer,
                    scan.TotalIonCurrent,
                    scan.InjectionTime,
                    scan.NoiseData,
                    scan.NativeId.Replace(scan.OneBasedScanNumber.ToString(), scanNumberMap[scan.OneBasedScanNumber].ToString()),
                    scan.SelectedIonMZ,
                    scan.SelectedIonChargeStateGuess,
                    scan.SelectedIonIntensity,
                    scan.IsolationMz,
                    scan.IsolationWidth,
                    scan.DissociationType,
                    newOneBasedPrecursorScanNumber,
                    scan.SelectedIonMonoisotopicGuessMz,
                    scan.HcdEnergy
                );
                scansForTheNewFile.Add(newDataScan);
            }
            string replace = "_RT" + Math.Round(scans.Where(s => s.OneBasedScanNumber == startScan).First().RetentionTime, 2) + "-" +
                Math.Round(scans.Where(s => s.OneBasedScanNumber == endScan).First().RetentionTime, 2) + ".mzML";
            string outPath = origDataFile.Replace(".raw", replace).ToString();

            SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
                reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.FileChecksumType, reader.SourceFile.Uri, reader.SourceFile.Id, reader.SourceFile.FileName);


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scansForTheNewFile.ToArray(), sourceFile), outPath, false);

        }
    }
}
