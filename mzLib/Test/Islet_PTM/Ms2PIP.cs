using System;
using System.Diagnostics;
using Microsoft.Win32;
using System.Windows;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using MzLibUtil;
using CsvHelper.Configuration;
using CsvHelper;
using System.Globalization;
using System.Text.RegularExpressions;
using MassSpectrometry.MzSpectra;
using NUnit.Framework;
using Omics;
using Readers.SpectralLibrary;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Omics;
using System.Data.Entity.Core.Common.CommandTrees.ExpressionBuilder;
using System.IO;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;
using Chemistry;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;
using MassSpectrometry;

namespace Test
{
    public class Ms2PIP
    {
        public static string ParseSequenceForMs2PIP(string fullSequence)
        {
            var sb = new StringBuilder();
            var mods = SpectrumMatchFromTsv.ParseModifications(fullSequence);
            for (int i = 0; i < fullSequence.Length; i++)
            {
                char c = fullSequence[i];
                sb.Append(c);
                if (mods.ContainsKey(i))
                {
                    var modName = mods[i].Split(':')[1].Split(' ')[0];
                    var parsedModName = DeepLC.ParseModNameForDeepLC(modName);
                    sb.Append($"[{parsedModName}]");
                }
            }
            return sb.ToString();
        }

        [Test]
        public static void TestMs2PIPinput()
        {
            var allPeptideLow_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesLow_file = new PsmFromTsvFile(allPeptideLow_path);
            var allPeptidesLow = allPeptidesLow_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");
            //var allPeptidesHigh_path = @"E:\Islets\Brian_data\HighResTMT\noCali_LFgptmdFilterPruned\Task1-SearchTask\AllPeptides.psmtsv";
            //var allPeptidesHigh_file = new PsmFromTsvFile(allPeptidesHigh_path);
            //var allPeptidesHigh = allPeptidesHigh_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            //var allPeptidesHighNoMod = Ptm_tmt.GetUnmodifiedPeptides(allPeptidesHigh);
            //var allPeptidesWithModsHigh = Ptm_tmt.GetModifiedPeptides(allPeptidesHigh);

            var outPath = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\TargetMod_input_ms2pip.peprec";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine("spec_id modifications peptide charge");
                foreach (var peptide in Ptm_tmt.GetModifiedPeptides(allPeptidesLow))
                {
                    var parsedMods = DeepLC.ParseModsForDeepLC(peptide.FullSequence);
                    if (parsedMods == null) continue;
                    var outString = new List<string> { $"scan{peptide.Ms2ScanNumber}", parsedMods, peptide.BaseSeq, peptide.ChargeState.ToString() };
                    writer.WriteLine(string.Join(" ", outString));
                }
            }
        }

        [Test]
        public static void TestSpectralSimilarity()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            var allPeptides_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides_file = new PsmFromTsvFile(allPeptides_path);
            var allPeptides = allPeptides_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");
            var allPeptidesNoMod = allPeptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")) && !p.FullSequence.Contains("|"));// && p.Score >= 10
            var allPeptidesWithMods = allPeptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)))).ToList();

            var raw_directory = @"E:\Islets\Brian_data\Real_islets\Frxn";
            var raw_paths = Directory.GetFiles(raw_directory, "*.raw");//, SearchOption.AllDirectories
            var rawFilesDict = raw_paths.ToDictionary(r => r, r => MsDataFileReader.GetDataFile(r));

            var libraryPath = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\TargetMod_modelCID-TMT_ms2pip.msp";
            var spectralLibrary = new SpectralLibrary(new List<string> { libraryPath });
            var allLibrarySepctra = spectralLibrary.Results;
            var similarities = new List<double>();

            var similarityOutPath = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\SpectralSimilarity_TargetMod_modelCID-TMT_ms2pip_0.2Da.tsv";
            using (StreamWriter writer = new StreamWriter(similarityOutPath))
            {
                writer.WriteLine("FileName\tRawScan\tPeptide");
                foreach (var peptide in Ptm_tmt.GetModifiedPeptides(allPeptides))
                {
                    var rawFile = rawFilesDict.FirstOrDefault(kvp => kvp.Key.Contains(peptide.FileName)).Value;
                    var rawScan = rawFile.GetOneBasedScan(peptide.Ms2ScanNumber);
                    var precursorTol = new PpmTolerance(20);
                    var productTol = new AbsoluteTolerance(0.2);
                    var librarySpectrum = spectralLibrary.Results.FirstOrDefault(s => IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(s.Sequence) == peptide.BaseSeq && s.ChargeState == peptide.ChargeState && precursorTol.Within(s.PrecursorMz.ToMass(s.ChargeState), peptide.MonoisotopicMass));//Math.Round(s.PrecursorMz, 0) == Math.Round(peptide.PrecursorMass.ToMz(peptide.ChargeState), 0)
                    if (librarySpectrum != null)
                    {
                        var similarity = new SpectralSimilarity(rawScan.MassSpectrum, librarySpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, 20, false, tol: productTol);
                        similarities.Add(similarity.CosineSimilarity().Value);
                        var outString = new List<string> { $"{rawFile.SourceFile.FileName}",rawScan.OneBasedScanNumber.ToString(), peptide.FullSequence };
                        writer.WriteLine(string.Join("\t", outString));
                    }
                }
            }
        }

        [Test]
        public static void TestPeptideIntensity()
        {
            var allPeptideLow_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesLow_file = new PsmFromTsvFile(allPeptideLow_path);
            var allPeptidesLow = allPeptidesLow_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");

            var pairs = new List<(double ms1, double ms3)> { };
            foreach(var peptide in allPeptidesLow)
            {
                pairs.Add((peptide.PrecursorIntensity.Value, peptide.Intensities.Average()));
            }
        }

        [Test]
        public static void TestMs2PIP()
        {
            string fullSeq = "[Multiplex Label:TMT18 on X]GQAGPEGAAP[Common Biological:Hydroxylation on P]APEEDK[Multiplex Label:TMT18 on K]";
            var parsed = ParseSequenceForMs2PIP(fullSeq);
        }

        public static void RunPredictBatchCLI(string pythonExePath, string inputPsmPath, string psmFileType = null, string outputName = null, string outputFormat = "msp", bool addRetentionTime = false, bool addIonMobility = false, string model = "HCD", string modelDir = null, int? processes = null)
        {
            // Build command line arguments
            string args = $"-m ms2pip predict-batch \"{inputPsmPath}\"";

            if (!string.IsNullOrEmpty(psmFileType)) args += $" --psm-filetype {psmFileType}";

            if (!string.IsNullOrEmpty(outputName)) args += $" --output-name \"{outputName}\"";

            if (!string.IsNullOrEmpty(outputFormat)) args += $" --output-format {outputFormat}";

            if (addRetentionTime) args += " --add-retention-time";


            if (addIonMobility) args += " --add-ion-mobility";

            if (!string.IsNullOrEmpty(model)) args += $" --model {model}";

            if (!string.IsNullOrEmpty(modelDir)) args += $" --model-dir \"{modelDir}\"";

            if (processes.HasValue) args += $" --processes {processes.Value}";

            // Setup process
            var psi = new ProcessStartInfo
            {
                FileName = pythonExePath,
                Arguments = args,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };

            Console.WriteLine($"Running command: {pythonExePath} {args}");

            using (var process = Process.Start(psi))
            {
                string stdout = process.StandardOutput.ReadToEnd();
                string stderr = process.StandardError.ReadToEnd();
                process.WaitForExit();

                Console.WriteLine("Standard Output:\n" + stdout);
                if (!string.IsNullOrWhiteSpace(stderr))
                    Console.WriteLine("Standard Error:\n" + stderr);
            }
        }

        public static void CheckAndRunMs2Pip(string inputPsmPath, string pythonPath = null, string psmFileType = null, string outputName = null, string outputFormat = "msp", bool addRetentionTime = false, bool addIonMobility = false, string model = "HCD", string modelDir = null, int? processes = null)
        {
            if (pythonPath == null)
                pythonPath = FindPythonExe();
            CheckMs2PipInstalled(pythonPath);

            RunPredictBatchCLI(pythonPath, inputPsmPath, psmFileType, outputName, outputFormat,
                            addRetentionTime, addIonMobility, model, modelDir, processes);
        }

        private static string FindPythonExe()
        {
            // Try PATH
            try
            {
                var psi = new ProcessStartInfo
                {
                    FileName = "where",
                    Arguments = "python",
                    RedirectStandardOutput = true,
                    UseShellExecute = false,
                    CreateNoWindow = true
                };
                using (var proc = Process.Start(psi))
                {
                    string firstLine = proc.StandardOutput.ReadLine();
                    if (File.Exists(firstLine)) return firstLine;
                }
            }
            catch { }

            // Try common locations
            string[] guesses =
            {
            $@"C:\Users\{Environment.UserName}\AppData\Local\Programs\Python\Python311\python.exe",
            $@"C:\Users\{Environment.UserName}\AppData\Local\Programs\Python\Python312\python.exe",
            $@"C:\Users\{Environment.UserName}\Anaconda3\python.exe"
        };
            foreach (var path in guesses)
                if (File.Exists(path)) return path;

            throw new MzLibException("Python was not found. Please install Python from https://www.python.org/downloads/ and try again.");
        }

        private static bool CheckMs2PipInstalled(string pythonPath)
        {
            var psi = new ProcessStartInfo
            {
                FileName = pythonPath,
                Arguments = "-m ms2pip --help",
                RedirectStandardError = true,
                RedirectStandardOutput = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };

            try
            {
                using (var proc = Process.Start(psi))
                {
                    string output = proc.StandardOutput.ReadToEnd();
                    string error = proc.StandardError.ReadToEnd();
                    proc.WaitForExit();
                    return output.Contains("Usage:") || error.Contains("Usage:");
                }
            }
            catch
            {
                throw new MzLibException("MS²PIP is not installed in your Python environment.\n\nPlease run the following in Command Prompt:\n\npip install ms2pip");
            }
        }
    }


}

