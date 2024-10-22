using System.Diagnostics;
using System.Text;
using Easy.Common.Interfaces;
using MassSpectrometry;
using Readers;


namespace ISD
{
    public static class LabelCorrection
    {
        public static IEnumerable<MsDataScan> GetISDScans(this List<MsDataScan> scans, double ISDVoltage)
        {
            return scans.Where(i => i.ScanFilter.Contains($"sid={ISDVoltage}"));
        }
        public static IEnumerable<MsDataScan> GetMs1Scans(this List<MsDataScan> scans)
        {
            return scans.Where(i => i.ScanFilter.Contains("sid=15"));
        }
        public static IEnumerable<MsDataScan> InterleaveScans(this IEnumerable<MsDataScan> ms1s,
            IEnumerable<MsDataScan> ms2s)
        {
            return ms1s.Zip(ms2s, (f, s) => new[] { f, s }).SelectMany(f => f);
        }

        public static IEnumerable<MsDataScan> UpdateMs2MetaData(this IEnumerable<MsDataScan> scansList)
        {
            foreach (MsDataScan ms in scansList)
            {
                if (!ms.ScanFilter.Contains("sid=15"))
                {
                    var isolationWidth = ms.ScanWindowRange.Maximum - ms.ScanWindowRange.Minimum;
                    ms.MsnOrder = 2;
                    ms.SetIsolationMz(isolationWidth / 2);
                    ms.IsolationWidth = isolationWidth;
                    ms.SelectedIonMZ = isolationWidth / 2;
                }

                yield return ms;
            }
            

        }

        /*
        public static IEnumerable<MsDataScan> UpdateScanStringMetaData(this IEnumerable<MsDataScan> scans)
        {
            int oneBasedScanNumber = 1;
            string nativeBase = "controllerType=0 controllerNumber=1 scan=";
            foreach (MsDataScan scan in scans)
            {
                scan.SetOneBasedScanNumber(oneBasedScanNumber);
                int precursorScanNumber = scan.OneBasedScanNumber - 1;
                scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                scan.SetNativeID(nativeBase + oneBasedScanNumber);
                StringBuilder sb = new StringBuilder();
                sb.AppendJoin("", "+ ", "NSI ", "Full ms ", "[", scan.ScanWindowRange.Minimum, "-",
                    scan.ScanWindowRange.Maximum, "]");
                scan.ScanFilter = sb.ToString();
                oneBasedScanNumber++;
                yield return scan;
            }
        }
        */

        public static IEnumerable<MsDataScan> UpdateIsdScanMetaData(this IEnumerable<MsDataScan> scans)
        {
            int oneBasedScanNumber = 1;
            string nativeBase = "controllerType=0 controllerNumber=1 scan=";
            foreach (MsDataScan scan in scans)
            {
                scan.SetOneBasedScanNumber(oneBasedScanNumber);
                if (!scan.ScanFilter.Contains("sid=15"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 1;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                }
                    scan.SetOneBasedScanNumber(oneBasedScanNumber);
                    scan.SetNativeID(nativeBase + oneBasedScanNumber);
                    StringBuilder sb = new StringBuilder();
                    //sb.AppendJoin("", "+ ", "NSI ", "Full ms ", "[", scan.ScanWindowRange.Minimum, "-",
                    //    scan.ScanWindowRange.Maximum, "]");
                    //scan.ScanFilter = sb.ToString();
                    oneBasedScanNumber++;
                yield return scan;
            }
        }

        //public static IEnumerable<MsDataScan> UpdateMetaDataForAll(this IEnumerable<MsDataScan> scans)
        //{
        //    foreach (MsDataScan scan in scans)
        //    {
        //        if (!scan.ScanFilter.Contains("sid=15"))
        //        {
        //            int precursorScanNumber = scan.OneBasedScanNumber - 1;
        //            scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
        //            var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
        //            scan.MsnOrder = 2;
        //            scan.SetIsolationMz(isolationWidth / 2);
        //            scan.IsolationWidth = isolationWidth;
        //            scan.SelectedIonMZ = isolationWidth / 2;
        //        }
        //        yield return scan;
        //    }
        //}

        public static void MultipleVoltageFileConversion(string rawFilePath, string outputPath)
        {
            var file = new ThermoRawFileReader(rawFilePath);
            var scansFull = file.GetAllScansList();
            var ms1scan = scansFull.Where(s => s.ScanFilter.Contains("sid=15")).FirstOrDefault();
            var isolationWidth = ms1scan.ScanWindowRange.Maximum - ms1scan.ScanWindowRange.Minimum;
            var isolationMz = isolationWidth / 2 + ms1scan.ScanWindowRange.Minimum;
            foreach (MsDataScan scan in scansFull)
            {
                if (scan.ScanFilter.Contains("sid=60"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 1;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationMz);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationMz;
                }
                if (scan.ScanFilter.Contains("sid=80"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 2;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationMz);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationMz;
                }
                if (scan.ScanFilter.Contains("sid=100"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 3;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationMz);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationMz;
                }
            }
            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFileCombined = new GenericMsDataFile(scansFull.ToArray(), genericSourceFile);
            msFileCombined.ExportAsMzML(outputPath, false);
        }

    }
}
