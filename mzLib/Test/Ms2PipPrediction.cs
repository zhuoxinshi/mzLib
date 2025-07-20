using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Python.Runtime;
using NUnit.Framework;
using System.Diagnostics;
using Microsoft.Win32;
using System.Windows;
using System.IO;
using MzLibUtil;
using CsvHelper.Configuration;
using CsvHelper;
using Proteomics.PSM;
using Readers;
using static Nett.TomlObjectFactory;
using System.Globalization;
using System.Text.RegularExpressions;

namespace Test
{
    public class Ms2PipPrediction
    {
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

        public static void CheckAndRunMs2Pip(string inputPsmPath, string pythonPath = null, string psmFileType = null, string outputName = null, string outputFormat = "msp", bool addRetentionTime = false,bool addIonMobility = false, string model = "HCD", string modelDir = null, int? processes = null)
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

        [Test]
        public static void TestMs2Pip()
        {
            var inputPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_filtered_sub_inputTest.tsv";
            var outPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\test4.msp";
            CheckAndRunMs2Pip(inputPath, outputName: outPath);
        }
    }

    public class Ms2PipInput
    {
        public string peptidoform { get; set; }
        public string spectrum_id { get; set; }

        public Ms2PipInput(PsmFromTsv psmTsv)
        {
            string updatedFullSeq = psmTsv.FullSequence;
            if (psmTsv.FullSequence.Contains("substitution"))
            {
                updatedFullSeq = PsmValidations.ParseSubstitutedFullSequence(psmTsv.FullSequence);
            }
            peptidoform = ParseModsForMs2PipInput(updatedFullSeq) + "/" + psmTsv.PrecursorCharge;
            spectrum_id = "scan" + psmTsv.Ms2ScanNumber.ToString() + ": " + psmTsv.FullSequence;
        }
        public Ms2PipInput() { }

        public static string ParseModsForMs2PipInput(string fullSequence)
        {
            // Regex matches [anything:ModificationName on X]
            return Regex.Replace(
                fullSequence,
                @"\[[^\[\]:]*:[ ]*([A-Za-z]+)[^\[\]]*\]",
                m => $"[{m.Groups[1].Value}]"
            );
        }
    }

    public class Ms2PipInputFile : ResultFile<Ms2PipInput>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public Ms2PipInputFile() : base() { }
        public Ms2PipInputFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<Ms2PipInput>().ToList();
        }

        public string FullFileName { get; set; }
        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<Ms2PipInput>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }
}
