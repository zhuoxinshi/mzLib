using MassSpectrometry;
using Readers;

namespace ISD
{
    /// <summary>
    /// Prepares in-source dissociation (ISD) Thermo .raw files for TopFD.
    ///
    /// The ISD method repeats a four-scan cycle: a survey scan at source-induced dissociation energy
    /// sid=15, then fragmentation scans at sid=60, 80 and 100. The instrument records all four as MS1
    /// scans; the "sid=" token in the scan filter is the only thing that tells them apart.
    /// </summary>
    public static class LabelCorrection
    {
        /// <summary>
        /// Returns the scans acquired at one source energy, picked by the "sid=" token in the scan filter.
        /// </summary>
        public static IEnumerable<MsDataScan> GetISDScans(this List<MsDataScan> scans, double ISDVoltage)
        {
            return scans.Where(i => i.ScanFilter.Contains($"sid={ISDVoltage}"));
        }

        /// <summary>
        /// Relabels every sid=60/80/100 scan as an MS2 of the sid=15 survey scan that opens its cycle
        /// (1, 2 and 3 scans earlier). The isolation window is set to the survey scan's full m/z range.
        /// All scans are written to "[raw file name]_labelCorrected.mzML", next to the .raw file unless
        /// an output folder is given.
        /// </summary>
        public static void ISD60_80_100_LabelCorrectionFromRaw(string inputPath, string? outputFolder = null)
        {
            var file = MsDataFileReader.GetDataFile(inputPath);
            var scansFull = file.GetAllScansList();
            foreach (MsDataScan scan in scansFull)
            {
                if (scan.ScanFilter.Contains("sid=60"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 1;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var ms1scan = scansFull.Where(s => s.OneBasedScanNumber == precursorScanNumber).First();
                    var isolationWidth = ms1scan.ScanWindowRange.Maximum - ms1scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(ms1scan.ScanWindowRange.Minimum + isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = ms1scan.ScanWindowRange.Minimum + isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=80"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 2;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var ms1scan = scansFull.Where(s => s.OneBasedScanNumber == precursorScanNumber).First();
                    var isolationWidth = ms1scan.ScanWindowRange.Maximum - ms1scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(ms1scan.ScanWindowRange.Minimum + isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = ms1scan.ScanWindowRange.Minimum + isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=100"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 3;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var ms1scan = scansFull.Where(s => s.OneBasedScanNumber == precursorScanNumber).First();
                    var isolationWidth = ms1scan.ScanWindowRange.Maximum - ms1scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(ms1scan.ScanWindowRange.Minimum + isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = ms1scan.ScanWindowRange.Minimum + isolationWidth / 2;
                }
            }
            var outPath = inputPath.Replace(".raw", "_labelCorrected.mzML");
            if (outputFolder != null)
            {
                outPath = Path.Combine(outputFolder, Path.GetFileName(outPath));
            }
            SourceFile genericSourceFile3 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFileCombined = new GenericMsDataFile(scansFull.ToArray(), genericSourceFile3);
            msFileCombined.ExportAsMzML(outPath, false);
        }

        /// <summary>
        /// Splits one ISD .raw file by source energy into four files for separate TopFD runs:
        /// "_ms1" (sid=15 survey scans), "_isd60", "_isd80" and "_isd100". Scans are written unchanged,
        /// so every file holds MS1 scans only.
        /// </summary>
        public static void SeparateFilesForTopFD(string rawFilePath, string outputFolder)
        {
            var reader = new ThermoRawFileReader(rawFilePath);
            var scans = reader.GetAllScansList();
            var ms1Scans = scans.GetISDScans(15).ToArray();
            var isdScans60 = scans.GetISDScans(60).ToArray();
            var isdScans80 = scans.GetISDScans(80).ToArray();
            var isdScans100 = scans.GetISDScans(100).ToArray();

            var ms1outPath = System.IO.Path.Combine(outputFolder, reader.SourceFile.FileName.Replace(".raw", "_ms1.mzML"));
            var isd60outPath = System.IO.Path.Combine(outputFolder, reader.SourceFile.FileName.Replace(".raw", "_isd60.mzML"));
            var isd80outPath = System.IO.Path.Combine(outputFolder, reader.SourceFile.FileName.Replace(".raw", "_isd80.mzML"));
            var isd100outPath = System.IO.Path.Combine(outputFolder, reader.SourceFile.FileName.Replace(".raw", "_isd100.mzML"));
            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFileCombined = new GenericMsDataFile(ms1Scans, genericSourceFile);
            msFileCombined.ExportAsMzML(ms1outPath, false);
            GenericMsDataFile isd60File = new GenericMsDataFile(isdScans60, genericSourceFile);
            isd60File.ExportAsMzML(isd60outPath, false);
            GenericMsDataFile isd80File = new GenericMsDataFile(isdScans80, genericSourceFile);
            isd80File.ExportAsMzML(isd80outPath, false);
            GenericMsDataFile isd100File = new GenericMsDataFile(isdScans100, genericSourceFile);
            isd100File.ExportAsMzML(isd100outPath, false);
        }
    }
}
