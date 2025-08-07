using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using MassSpectrometry;
using System.IO;
using System.Text.RegularExpressions;
using System.Windows.Interop;

namespace Test
{
    public class nothing
    {
        [Test]
        public static void RNA()
        {
            var path1 = @"E:\Temp\20250612_21mer#1-Am.raw";
            var path2 = @"E:\Temp\20250612_22mer#1-m6Am.raw";
            var path3 = @"E:\Temp\20250612_22mer#1-unmod.raw";
            var path4 = @"E:\Temp\20250612_RNA-Mix_10V.raw";
            var path5 = @"E:\Temp\20250612_RNA-Mix_30V.raw";
            var path6 = @"E:\Temp\20250612-22mer#1-Am.raw";
            var path7 = @"E:\Temp\20250612-22mer#1-m6A.raw";
            var path8 = @"E:\Temp\20250616_21mer#1-Am.raw";
            var paths = new List<string> { path1, path2, path3, path4, path5, path6, path7, path8 };
            var dictionary = new Dictionary<string, MsDataScan[]>();
            foreach (var path in paths)
            {
                var msData = MsDataFileReader.GetDataFile(path);
                var scans = msData.GetAllScansList().Where(S => S.MsnOrder == 2).ToArray();
                dictionary[Path.GetFileNameWithoutExtension(path)] = scans;
            }
            var osmPath = @"E:\Temp\AllOSMs.osmtsv";
            var allOsms = SpectrumMatchTsvReader.ReadOsmTsv(osmPath, out List<string> warnings).ToList();

            var cidEnergys = new List<double>();
            string pattern = $@"cid(\d+)";
            foreach (var osm in allOsms)
            {
                var ms2Scan = dictionary[osm.FileName].First(s => s.OneBasedScanNumber == osm.Ms2ScanNumber);
                var match = Regex.Match(ms2Scan.ScanFilter, pattern);
                double voltage = double.Parse(match.Groups[1].Value);
                cidEnergys.Add(voltage);
            }

            // Write to TSV file with one column
            var outPath = @"E:\Temp\cidEnergys.tsv";
            using (StreamWriter writer = new StreamWriter(outPath)) 
            { writer.WriteLine("CID"); // Header
            foreach (var value in cidEnergys) { 
                    writer.WriteLine(value); 
                } }
        }

        [Test]
        public static void Mod()
        {
            var allmods = ModificationConverter.AllKnownMods;
            var substitutionMods = ModificationConverter.AllKnownMods.Where(m => m.ModificationType.Contains("substitution")).ToList();
        }
    }
}
