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

namespace Test.IsdTests
{
    internal class DataLoading
    {

        [Test]
        public void TestImportExport2()
        {
            string path = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample13_10uL_5pmol_ISD_IT.raw";
            string path2 = @"E:\ISD Project\ISD_240606\06-11-24_mix_sample13_10uL_5pmol_ISD_IT.mzML";

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

    }
}
