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

namespace Proteomics
{
    public class Ms2PIP
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
