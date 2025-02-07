using Chemistry;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;
using Readers;
using ISD;
using MassSpectrometry;
using UsefulProteomicsDatabases.Generated;
using static System.Net.WebRequestMethods;
using System.Text.RegularExpressions;
using System.Security.Cryptography.X509Certificates;
using System.Runtime.CompilerServices;
using System.Windows.Shapes;
using Nett;
using SpectralAveraging;
using ISD;
using System.Windows.Media.Animation;

namespace Test.IsdTests
{
    internal class DataLoading
    {

        [Test]
        public void TestImportExport2()
        {
            string path = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.raw";
            string path2 = @"E:\ISD Project\ISD_240606\06-07-24_mix_1pmol_5uL_ISD.mzML";

            var file = new ThermoRawFileReader(path);
            var scansFull = file.GetAllScansList();

            var ms1Scans = scansFull.GetMs1Scans();
            var voltage = 100;
            var isdScans = scansFull.GetISDScans(voltage).ToList();
            var interleaved = ms1Scans.InterleaveScans(isdScans).ToList();

            var results = interleaved.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();

            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile = new GenericMsDataFile(results, genericSourceFile);
            msFile.ExportAsMzML(path2, false);
        }

        [Test]
        public void LoadInAllScansInFile()
        {
            string path = @"E:\ISD Project\ISD_240305\HighResOrbiMS1_ISD100.mzML";
            var file = MsDataFileReader.GetDataFile(path);
            var scans = file.GetAllScansList();
        }


        [Test]
        [TestCase(@"E:\ISD Project\ISD_240305\FileList.txt")]
        public void ConvertEntireList(string directoryPath)
        {
            using var streamReader = new StreamReader(directoryPath);

            string? name;
            while ((name = streamReader.ReadLine()) is not null)
            {
                name = name.Trim();
                string folderName = @"ISD_240305\";
                string outputPath = @"E:\ISD Project\" + folderName + name + ".mzML";

                string inputPath = @"E:\ISD Project\" + folderName + name +".raw";
                var file = new ThermoRawFileReader(inputPath);
                var scansFull = file.GetAllScansList();

                var ms1Scans = scansFull.GetMs1Scans();
                Regex v = new Regex(@"\d{2,3}$");
                Match m = v.Match(name);
                int voltage = int.Parse(m.Value);
                var isdScans = scansFull.GetISDScans(voltage).ToList();
                var interleaved = ms1Scans.InterleaveScans(isdScans).ToList();

                var results = interleaved.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();

                SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                    null, null, null);
                GenericMsDataFile msFile = new GenericMsDataFile(results, genericSourceFile);
                msFile.ExportAsMzML(outputPath, false);
            }


        }

        [Test]
        public void TestFileReading()
        {
            string path3 = @"E:\ISD Project\ISD_240305\HighResOrbiMS1_ISD75.mzML";
            var reader = MsDataFileReader.GetDataFile(path3);
            reader.LoadAllStaticData();
        }

        [Test]
        public void FileConversion0418()
        {
            string inputPath = @"E:\ISD Project\ISD_240606\06-06-24_mix_higherConc_5uL_ISD.raw";
            string path75 = @"E:\ISD Project\ISD_240606\06-06-24_mix_higherConc_5uL_ISD_75V.mzML";
            
            var file = new ThermoRawFileReader(inputPath);
            var scansFull = file.GetAllScansList();
            var ms1Scans = scansFull.GetMs1Scans();
            var isdScans75 = scansFull.GetISDScans(75).ToList();
            var interleaved75 = ms1Scans.InterleaveScans(isdScans75).ToList();
            var results75 = interleaved75.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile75 = new GenericMsDataFile(results75, genericSourceFile);
            msFile75.ExportAsMzML(path75, false);

            string path100 = @"E:\ISD Project\ISD_240606\06-06-24_mix_higherConc_5uL_ISD_100V.mzML";
            var file2 = new ThermoRawFileReader(inputPath);
            var scansFull2 = file2.GetAllScansList();
            var ms1Scans2 = scansFull2.GetMs1Scans();
            var isdScans100 = scansFull2.GetISDScans(100).ToList();
            var interleaved100 = ms1Scans2.InterleaveScans(isdScans100).ToList();
            var results100 = interleaved100.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();

            SourceFile genericSourceFile2 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile100 = new GenericMsDataFile(results100, genericSourceFile2);
            msFile100.ExportAsMzML(path100, false);
        }

        [Test]
        public void ConvertFileList()
        {
            string fileListPath = @"E:\ISD Project\ISD_240606\FileList2.txt";

            string[] names = System.IO.File.ReadAllLines(fileListPath);

            foreach (string name in names)
            {
                string folderName = @"ISD_240606\";
                string outputPath = @"E:\ISD Project\" + folderName + name + ".mzML";

                string inputPath = @"E:\ISD Project\" + folderName + name + ".raw";
                var file = new ThermoRawFileReader(inputPath);
                var scansFull = file.GetAllScansList();

                var ms1Scans = scansFull.GetMs1Scans();
                int voltage = 100;
                var isdScans = scansFull.GetISDScans(voltage).ToList();
                var interleaved = ms1Scans.InterleaveScans(isdScans).ToList();

                var results = interleaved.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();

                SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                    null, null, null);
                GenericMsDataFile msFile = new GenericMsDataFile(results, genericSourceFile);
                msFile.ExportAsMzML(outputPath, false);
            }


        }

        [Test]
        public void FileConversionForMultipleVoltages()
        {
            string inputPath = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1.raw";
            string path40only = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1_40only.mzML";
            string path50only = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1_50only.mzML";
            string path60only = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1_60only.mzML";
            string path80only = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1_80only.mzML";
            string path100only = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionB_orbiMS1_ISD40-50-60-80-100_micro1_100only.mzML";

            var file = new ThermoRawFileReader(inputPath);
            var scansFull = file.GetAllScansList();
            var ms1Scans = scansFull.GetMs1Scans();
            var isdScans60 = scansFull.GetISDScans(60).ToList();
            var interleaved60 = ms1Scans.InterleaveScans(isdScans60).ToList();
            var results60 = interleaved60.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile60 = new GenericMsDataFile(results60, genericSourceFile);
            msFile60.ExportAsMzML(path60only, false);

            var file2 = new ThermoRawFileReader(inputPath);
            var scansFull2 = file2.GetAllScansList();
            var ms1Scans2 = scansFull2.GetMs1Scans();
            var isdScans100 = scansFull2.GetISDScans(100).ToList();
            var interleaved100 = ms1Scans2.InterleaveScans(isdScans100).ToList();
            var results100 = interleaved100.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile2 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile100 = new GenericMsDataFile(results100, genericSourceFile2);
            msFile100.ExportAsMzML(path100only, false);

            var file3 = new ThermoRawFileReader(inputPath);
            var scansFull3 = file3.GetAllScansList();
            var ms1Scans3 = scansFull3.GetMs1Scans();
            var isdScans80 = scansFull3.GetISDScans(80).ToList();
            var interleaved80 = ms1Scans3.InterleaveScans(isdScans80).ToList();
            var results80 = interleaved80.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile3 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile80 = new GenericMsDataFile(results80, genericSourceFile3);
            msFile80.ExportAsMzML(path80only, false);

            var file4 = new ThermoRawFileReader(inputPath);
            var scansFull4 = file4.GetAllScansList();
            var ms1Scans4 = scansFull4.GetMs1Scans();
            var isdScans40 = scansFull4.GetISDScans(40).ToList();
            var interleaved40 = ms1Scans4.InterleaveScans(isdScans40).ToList();
            var results40 = interleaved40.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile4 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile40 = new GenericMsDataFile(results40, genericSourceFile4);
            msFile40.ExportAsMzML(path40only, false);

            var file5 = new ThermoRawFileReader(inputPath);
            var scansFull5 = file5.GetAllScansList();
            var ms1Scans5 = scansFull5.GetMs1Scans();
            var isdScans50 = scansFull5.GetISDScans(50).ToList();
            var interleaved50 = ms1Scans5.InterleaveScans(isdScans50).ToList();
            var results50 = interleaved50.UpdateMs2MetaData().UpdateIsdScanMetaData().ToArray();
            SourceFile genericSourceFile5 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile50 = new GenericMsDataFile(results50, genericSourceFile5);
            msFile50.ExportAsMzML(path50only, false);
        }

        [Test]
        public void FileConversionForMultipleVoltagesCombined()
        {
            string inputPath = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_averagedAll_0.5mzstep.mzML";
            string pathCombined = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_averagedAll_0.5mzstep_relabled.mzML";
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
            SourceFile genericSourceFile3 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFileCombined = new GenericMsDataFile(scansFull.ToArray(), genericSourceFile3);
            msFileCombined.ExportAsMzML(pathCombined, false);
        }

        [Test]
        public static void Convert0927()
        {
            string fileListPath = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionABCDE_orbiMS1_tech-rep.txt";

            string[] names = System.IO.File.ReadAllLines(fileListPath);

            foreach (string name in names)
            {
                string folderName = @"ISD_240927\";
                string outputPath = @"E:\ISD Project\" + folderName + name + ".mzML";

                string inputPath = @"E:\ISD Project\" + folderName + name + ".raw";
            }
        }

        [Test]
        public void FileConversion5voltages()
        {
            string inputPath = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionD_orbiMS1_ISD40-50-60-80-100_micro4.raw";
            string pathCombined = @"E:\ISD Project\ISD_240812\test.mzML";
            var file = MsDataFileReader.GetDataFile(inputPath);
            var scansFull = file.GetAllScansList();
            foreach (MsDataScan scan in scansFull)
            {
                if (scan.ScanFilter.Contains("sid=40"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 1;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=50"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 2;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=60"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 3;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=80"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 4;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationWidth / 2;
                }
                if (scan.ScanFilter.Contains("sid=100"))
                {
                    int precursorScanNumber = scan.OneBasedScanNumber - 5;
                    scan.SetOneBasedPrecursorScanNumber(precursorScanNumber);
                    var isolationWidth = scan.ScanWindowRange.Maximum - scan.ScanWindowRange.Minimum;
                    scan.MsnOrder = 2;
                    scan.SetIsolationMz(isolationWidth / 2);
                    scan.IsolationWidth = isolationWidth;
                    scan.SelectedIonMZ = isolationWidth / 2;
                }
            }
            SourceFile genericSourceFile3 = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFileCombined = new GenericMsDataFile(scansFull.ToArray(), genericSourceFile3);
            msFileCombined.ExportAsMzML(pathCombined, false);
        }

        [Test]
        public static void CountNumPeaks()
        {
            string inputPath = @"E:\ISD Project\ISD_241001\10-03-24_PEPPI_FractionD_orbiMS1_ISD40-50-60-80-100_micro4.raw";
            var file = new ThermoRawFileReader(inputPath);
            var scansFull = file.GetAllScansList();
            var isdScans40 = scansFull.GetISDScans(40).ToList();
            var isdScans50 = scansFull.GetISDScans(50).ToList();
            var isdScans60 = scansFull.GetISDScans(60).ToList();
            var isdScans80 = scansFull.GetISDScans(80).ToList();
            var isdScans100 = scansFull.GetISDScans(100).ToList();
        }

        [Test]
        public static void TestTopDIADecon()
        {
            string file = @"E:\DIA\TopDIA\20231117_DIA_720_800_rep2_RT46.72-48.71.mzML";
            var reader = MsDataFileReader.GetDataFile(file);
            var ms1scans = reader.GetMS1Scans();
            var deconParams = new ClassicDeconvolutionParameters(6, 60, 20, 3, Polarity.Positive);
            foreach(var scan in ms1scans)
            {
                var envelope = Deconvoluter.Deconvolute(scan, deconParams, new MzLibUtil.MzRange(720,733)).ToArray();
            }
        }

        [Test]
        public static void SeparateFilesForTopFD()
        {
            string file = @"E:\ISD Project\ISD_240812\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100.raw";
            var reader = new ThermoRawFileReader(file);
            var scans = reader.GetAllScansList();
            var ms1Scans = scans.GetISDScans(15).ToArray();
            var isdScans60 = scans.GetISDScans(60).ToArray();
            var isdScans80 = scans.GetISDScans(80).ToArray();
            var isdScans100 = scans.GetISDScans(100).ToArray();
            var ms1outPath = @"E:\toppic-windows-1.7.4\ISD_yeast_search\separateVoltageFile\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_MS1.mzML";
            var isd60outPath = @"E:\toppic-windows-1.7.4\ISD_yeast_search\separateVoltageFile\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_60.mzML";
            var isd80outPath = @"E:\toppic-windows-1.7.4\ISD_yeast_search\separateVoltageFile\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_80.mzML";
            var isd100outPath = @"E:\toppic-windows-1.7.4\ISD_yeast_search\separateVoltageFile\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_100.mzML";
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

        [Test]
        public static void TestReadFeatureFile()
        {
            var file = @"E:\toppic-windows-1.6.5\ISD\test1\id_08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_60_ms1.feature";
            var featureFile = new Ms1FeatureFile(file);
            //var feature = featureFile.First();
            featureFile.LoadResults();
        }

        [Test]
        public static void Test0927()
        {
            var path1 = @"E:\ISD Project\ISD_250128\01-28-25_td-ISD_PEPPI-YB_105min_ISD60-80-100_120k_micro1.raw";
            var path2 = @"E:\ISD Project\ISD_250128\01-28-25_td-ISD_PEPPI-YB_105min_ISD60-80-100_60k_micro1.raw";
            var path3 = @"E:\ISD Project\ISD_250128\01-28-25_td-DIA_PEPPI-YB_105min_50mz_21-23-25HCD_AGC1e6_200ms.raw";
            var path4 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionC_orbiMS1_ISD60.raw";
            var path5 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionC_orbiMS1_ISD80.raw";
            var path6 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionC_orbiMS1_ISD100.raw";
            var path7 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionD_orbiMS1_ISD60.raw";
            var path8 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionD_orbiMS1_ISD80.raw";
            var path9 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionD_orbiMS1_ISD100.raw";
            var path10 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionA_orbiMS1_ISD60.raw";
            var path11 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionA_orbiMS1_ISD80.raw";
            var path12 = @"E:\ISD Project\ISD_240927\09-27-24_PEPPI_FractionA_orbiMS1_ISD100.raw";
            var paths60 = new List<string> { path1, path4, path5, path7, path10 };
            var paths80 = new List<string> { path2, path5, path8, path11 };
            var paths100 = new List<string> { path3, path6, path9, path12 };
            foreach (var path in paths60)
            {
                LabelCorrection.ISDSingleVoltageLabelCorrectionFromRaw(path, 60);
            }
            foreach(var path in paths80)
            {
                LabelCorrection.ISDSingleVoltageLabelCorrectionFromRaw(path, 80);
            }
            foreach(var path in paths100)
            {
                LabelCorrection.ISDSingleVoltageLabelCorrectionFromRaw(path, 100);
            }
        }

        [Test]
        public static void Test250128()
        {
            var path1 = @"E:\ISD Project\ISD_250128\02-03-25_td-ISD_PEPPI-YC_105min_ISD60-80-100_0-1ug.raw";
            var path2 = @"E:\ISD Project\ISD_250128\02-03-25_td-ISD_PEPPI-YC_105min_ISD60-80-100_0-05ug.raw";
            var path3 = @"E:\ISD Project\ISD_250128\02-03-25_td-ISD_PEPPI-YC_105min_ISD60-80-100_0-2ug.raw";
            var path4 = @"E:\ISD Project\ISD_250128\02-03-25_td-ISD_PEPPI-YC_105min_ISD60-80-100_0-5ug.raw";
            var paths = new List<string> { path1, path2, path3, path4};
            foreach (var path in paths)
            {
                LabelCorrection.ISD60_80_100_LabelCorrectionFromRaw(path);
            }
            LabelCorrection.ISDSingleVoltageLabelCorrectionFromRaw(path4, 60);
        }
    }
}
