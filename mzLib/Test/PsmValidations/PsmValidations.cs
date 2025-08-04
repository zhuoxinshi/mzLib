using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics.RetentionTimePrediction.Chronologer;
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownProteomics;
using System.Text.RegularExpressions;
using Omics;
using MathNet.Numerics.Statistics;
using System.IO;
using Readers;
using Proteomics;
using Plotly;
using Plotly.NET;
using Readers.SpectralLibrary;
using MassSpectrometry;
using Chemistry;
using MassSpectrometry.MzSpectra;
using Omics.SpectrumMatch;

namespace Test
{
    public class PsmValidations
    {
        [Test]
        public static void OneSubOnly()
        {
            var psmFilePath_sub_1614 = @"E:\Aneuploidy\DDA\071525\1614_E1-8_calied-generalGPTMD+1NAsub_noTrunc\Task2-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_PSMs.psmtsv";
            var psmtsv_sub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_sub_1614, out List<string> warnings).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var psmFilePath_noSub_1614 = @"E:\Aneuploidy\DDA\071525\1614_E1-8_cali-generalGPTMD_noTrunc\Task3-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_PSMs.psmtsv";
            var psmtsv_noSub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_noSub_1614, out List<string> warnings2).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();

            //candidate PSMs that have a substitution
            var psmtsvWithSub_1614 = psmtsv_sub_1614.Where(p => p.FullSequence.Contains("substitution")).ToList();

            //filter out PSMs where the predicted RT is considered as an outlier
            var rtFilteredPsms_1614 = FilterPsmTsvFromPredictedRT(psmtsvWithSub_1614, 1, 1.68, out List<(int, string, double, float)> filteredPredictions);
            //PlotPredictedRt(filteredPredictions).Show();

            //filter out PSMs that can be explained by other PTMs
            var filteredPsm_sub_1614 = FilterUniquePsmTsv(rtFilteredPsms_1614, psmtsv_noSub_1614);
            var filteredPep_sub_1614 = filteredPsm_sub_1614.GroupBy(p => p.FullSequence).Select(g => g.First()).ToList();
            var pepToWrite = filteredPep_sub_1614.Where(p => SpectrumMatchFromTsv.ParseModifications(Ms2PipInput.ParseSubstitutedFullSequence(p.FullSequence)).Values.SelectMany(v => v).All(mod => mod == "Common Fixed:Carbamidomethyl on C" || mod == "Common Variable:Oxidation on M")).ToList();

            //write Ms2Pip input file for spectral prediction
            var libraryOutPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_HCDch2_oneSub_ms2pip.msp";
            var inputFilePath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_oneSub_ms2PipInput.tsv";
            if (!File.Exists(libraryOutPath))
            {
                WriteMs2PipInputFileFromPsmTsv(pepToWrite, inputFilePath);
                Ms2PIP.CheckAndRunMs2Pip(inputFilePath, null, null, libraryOutPath, "msp", false, false, "HCDch2", null);
            }

            //Calculate spectal similarity
            var pathList = new List<string> { libraryOutPath };
            var library = new SpectralLibrary(pathList);
            var librarySpectra = library.GetAllLibrarySpectra().ToList();
            var rawPath = @"E:\Aneuploidy\DDA\071525\07-15-25_1614-R1-Q_E1+5-calib.mzML";
            var rawFile = MsDataFileReader.GetDataFile(rawPath);
            var ms2Scans = rawFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var filteredPsms = new List<PsmFromTsv>();
            foreach (var psmTsv in filteredPsm_sub_1614)
            {
                var substitutedFullSeq = Ms2PipInput.ParseSubstitutedFullSequence(psmTsv.FullSequence);
                if (library.TryGetSpectrum(substitutedFullSeq, psmTsv.PrecursorCharge, out LibrarySpectrum libSpectrum))
                {
                    var rawScan = ms2Scans.FirstOrDefault(s => s.OneBasedScanNumber == psmTsv.Ms2ScanNumber);
                    var similarity = new SpectralSimilarity(rawScan.MassSpectrum, libSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, 20, false);
                    if (similarity.CosineSimilarity() >= 0.7) 
                    {
                        filteredPsms.Add(psmTsv);
                    }
                }
            }

            var fileteredScanOutPath = @"E:\Aneuploidy\DDA\071525\1614_oneSub_denovo.mzML";
            WriteOutScansForDenovo(filteredPsms, rawPath, fileteredScanOutPath);
        }

        [Test]
        public static void TestNormal()
        {
            var psmFilePath_sub_1614 = @"E:\Aneuploidy\DDA\071525\1614_E1-8_calied-generalGPTMD+1NAsub_noTrunc\Task2-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_PSMs.psmtsv";
            var psmtsv_sub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_sub_1614, out List<string> warnings).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var psmFilePath_noSub_1614 = @"E:\Aneuploidy\DDA\071525\1614_E1-8_cali-generalGPTMD_noTrunc\Task3-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_PSMs.psmtsv";
            var psmtsv_noSub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_noSub_1614, out List<string> warnings2).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();

            var psmsWithMod = psmtsv_sub_1614.Where(p => p.FullSequence.Contains("[")).ToList();
            var testPep = psmsWithMod.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.SelectMany(v => v).All(mod => mod == "Common Fixed:Carbamidomethyl on C" || mod == "Common Variable:Oxidation on M")).GroupBy(p => p.FullSequence).Select(g => g.First()).ToList();
            var rawPath = @"E:\Aneuploidy\DDA\071525\07-15-25_1614-R1-Q_E1+5-calib.mzML";

            var fileteredScanOutPath = @"E:\Aneuploidy\DDA\071525\1614_testNormal_denovo.mzML";
            WriteOutScansForDenovo(testPep, rawPath, fileteredScanOutPath);
        }
        [Test]
        public static void CasanovoOneSubOnly()
        {
            var path = @"E:\Aneuploidy\DDA\071525\casanovo_20250731231227.mztab";
            var casanovoResults = CasanovoResultFile.ReadInCasanovoResults(path.Replace(".mzML", ".mztab"));
            var filteredResults = casanovoResults.Where(r => IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(Ms2PipInput.ParseSubstitutedFullSequence(r.FullSequenceFromMM)) == r.DenovoSequence).ToList();
        }

        [Test]
        public static void TestWrittenFileReading()
        {
            var path = @"E:\Aneuploidy\DDA\071525\1614_oneSub_denovo.mzML";
            var rawFile = MsDataFileReader.GetDataFile(path);
            var scans = rawFile.GetAllScansList().ToList();
        }

        [Test]
        public static void TrunctatedPeptides()
        {
            var pepFilePath_sub_1614 = @"E:\Aneuploidy\DDA\062525\1614_E1-8_cali-gptmd(bioMods+1NAsub)-xml+trunc\Task3-SearchTask\Individual File Results\06-26-25_1614-R1-Q_E1+5-calib_Peptides.psmtsv";
            var pep_sub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(pepFilePath_sub_1614, out List<string> warnings).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var psmFilePath_noSub_1614_pep = @"E:\Aneuploidy\DDA\062525\1614_E1-8_calied-gptmd(bioMods)-xml+trunc\Task2-SearchTask\Individual File Results\06-26-25_1614-R1-Q_E1+5-calib_Peptides.psmtsv";
            var pep_noSub_1614 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_noSub_1614_pep, out List<string> warnings2).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();

            var pepFilePath_sub_1611 = @"E:\Aneuploidy\DDA\062525\1611_E1-8_cali-gptmd(bioMods+1NAsub)-xml+trunc\Task3-SearchTask\Individual File Results\06-25-25_1611-R1-Q_E1+5-calib_Peptides.psmtsv";
            var pep_sub_1611 = SpectrumMatchTsvReader.ReadPsmTsv(pepFilePath_sub_1611, out List<string> warnings3).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var psmFilePath_noSub_1611_pep = @"E:\Aneuploidy\DDA\062525\1611_E1-8_calied-gptmd(bioMods)-xml+trunc\Task2-SearchTask\Individual File Results\06-25-25_1611-R1-Q_E1+5-calib_Peptides.psmtsv";
            var pep_noSub_1611 = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath_noSub_1611_pep, out List<string> warnings4).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();

            var pep_sub_1614_trunc = pep_sub_1614.Where(p => p.BaseSeq.Last() != 'K' && p.BaseSeq.Last() != 'R').ToList();
            var pep_noSub_1614_trunc = pep_noSub_1614.Where(p => p.BaseSeq.Last() != 'K' && p.BaseSeq.Last() != 'R').ToList();
            var pep_sub_1611_trunc = pep_sub_1611.Where(p => p.BaseSeq.Last() != 'K' && p.BaseSeq.Last() != 'R').ToList();
            var pep_noSub_1611_trunc = pep_noSub_1611.Where(p => p.BaseSeq.Last() != 'K' && p.BaseSeq.Last() != 'R').ToList();
            //psm.BaseSeq.Last() != 'K' && psm.BaseSeq.Last() != 'R'
        }

        public static List<PsmFromTsv> FilterUniquePsmTsv(List<PsmFromTsv> psmsWithSub, List<PsmFromTsv> psmsWithoutSub)
        {
            var filteredPsms = new List<PsmFromTsv>();
            foreach (var psm in psmsWithSub)
            {
                if (psm.FullSequence.Contains("|"))
                {
                    continue;
                }
                if (psmsWithoutSub.Any(p => p.Ms2ScanNumber == psm.Ms2ScanNumber && p.PrecursorCharge == psm.PrecursorCharge && Math.Round(p.PrecursorMass, 2) == Math.Round(psm.PrecursorMass, 2)))
                {
                    continue;
                }
                filteredPsms.Add(psm);
            }
            return filteredPsms;
        }

        public static List<PsmFromTsv> FilterPsmTsvFromPredictedRT(List<PsmFromTsv> psms, double rtWindow, double zScoreCutOff, out List<(int, string, double, float)> filteredPredictions)
        {
            var filteredPsms = new List<PsmFromTsv>();
            var predictions = new List<(int, string, double, float)>();
            filteredPredictions = new List<(int, string, double, float)>();
            foreach (var psm in psms)
            {
                if (psm.FullSequence.Contains("|"))
                {
                    continue;
                }
                var allMods = SpectrumMatchFromTsv.ParseModifications(psm.FullSequence).Values.SelectMany(m => m).ToList();
                if (allMods.Any(m => m.Contains("substitution")))
                {
                    string modifiedFullSeq = Ms2PipInput.ParseSubstitutedFullSequence(psm.FullSequence);
                    string modifiedBaseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(modifiedFullSeq);
                    var predictedRt = ChronologerEstimator.PredictRetentionTime(modifiedBaseSeq, modifiedFullSeq);
                    predictions.Add((psm.Ms2ScanNumber, psm.FullSequence, psm.RetentionTime, predictedRt));
                }
                else
                {
                    var predictedRt = ChronologerEstimator.PredictRetentionTime(psm.BaseSeq, psm.FullSequence);
                    predictions.Add((psm.Ms2ScanNumber, psm.FullSequence, psm.RetentionTime, predictedRt));
                }
            }
            foreach (var prediction in predictions)
            {
                var localPredictions = predictions.Where(p => Math.Abs(p.Item3 - prediction.Item3) <= rtWindow).ToList();
                var zScore = (prediction.Item4 - localPredictions.Select(p => p.Item4).Average()) / localPredictions.Select(p => p.Item4).StandardDeviation();
                if (Math.Abs(zScore) < zScoreCutOff)
                {
                    filteredPsms.Add(psms.FirstOrDefault(p => p.Ms2ScanNumber == prediction.Item1 && p.FullSequence == prediction.Item2));
                    filteredPredictions.Add(prediction);
                }
            }
            return filteredPsms;
        }

        public static void WriteOutScansForDenovo(List<PsmFromTsv> psmTsvs, string rawFilePath, string outPath)
        {
            var rawFile = MsDataFileReader.GetDataFile(rawFilePath);
            var sourceFile = rawFile.GetSourceFile();
            var ms2Scans = rawFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var filteredScans = new List<MsDataScan>();
            int oneBasedNumber = 1;
            foreach(var psmTsv in psmTsvs)
            {
                var rawScan = ms2Scans.Where(s => s.OneBasedScanNumber == psmTsv.Ms2ScanNumber).FirstOrDefault();
                var newScan = new MsDataScan(rawScan.MassSpectrum, oneBasedNumber, 2, true, Polarity.Positive, rawScan.RetentionTime, rawScan.ScanWindowRange, rawScan.ScanFilter, rawScan.MzAnalyzer, rawScan.TotalIonCurrent, rawScan.InjectionTime, rawScan.NoiseData, "scan: " + psmTsv.Ms2ScanNumber + " "+ psmTsv.FullSequence, psmTsv.PrecursorMz, psmTsv.PrecursorCharge, psmTsv.PrecursorIntensity.Value, psmTsv.PrecursorMz, rawScan.IsolationWidth, rawScan.DissociationType, null, psmTsv.PrecursorMass.ToMz(psmTsv.PrecursorCharge));
                filteredScans.Add(newScan);
                oneBasedNumber++;
            }
            var dataFile = new GenericMsDataFile(filteredScans.ToArray(), sourceFile);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(dataFile, outPath, true);
        }

        public static GenericChart PlotPredictedRt(List<(int, string, double, float)> predictions)
        {
            var scatter = Chart2D.Chart.Point<double, float, string>(
                                x: predictions.Select(p => p.Item3),
                                y: predictions.Select(p => p.Item4)).WithMarkerStyle(Color: Color.fromString("blue"));
            return scatter;
        }

        [Test]
        public static void TestParsing()
        {
            //var seq = "S[1 nucleotide substitution:S->C on S]QEELDEMGAPIDYLTPIVADADAGHGGLTAVFK";
            //var modifiedSequence = ParseSubstitutedBaseSequence(seq);
            //var mod4 = ParseSubstitutedFullSequence(seq);
            var seq2 = "IVTEDC[Common Fixed:Carbamidomethyl on C]F[1 nucleotide substitution:F->Y on F]LQIDQSAITGESLAAEK";
            //var modifiedSequence2 = ParseSubstitutedBaseSequence(seq2);
            var mod3 = Ms2PipInput.ParseSubstitutedFullSequence(seq2);
        }

        [Test]
        public static void PredictSpectra()
        {
            var libraryOutPath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_noMod_HCD_predictions.msp";
            var inputFilePath = @"E:\Aneuploidy\DDA\062525\RtPredictionResults\1614_noMod_ms2PipInput.tsv";
            Ms2PIP.CheckAndRunMs2Pip(inputFilePath, null, null, libraryOutPath, "msp", false, false, "HCD", null);
        }

        //public static (double, double) BuildCalibrationCurveFromUnmodifiedPeptides(List<PsmFromTsv> peptides)
        //{
        //    var calibration_1614 = new List<(double, float)>();
        //    foreach (var pep in peptides)
        //    {
        //        var prediction = ChronologerEstimator.PredictRetentionTime(pep.BaseSeq, pep.FullSequence);
        //        calibration_1614.Add((pep.RetentionTime.Value, prediction));
        //    }
        //    var (intercept, slope) = Fit.Line(calibration_1614.Select(p => p.Item1).ToArray(), calibration_1614.Select(p => (double)p.Item2).ToArray());
        //    return (intercept, slope);
        //}

        public static void WriteMs2PipInputFileFromPsmTsv(List<PsmFromTsv> psmTsvList, string outPath)
        {
            var inputs = new List<Ms2PipInput>();
            foreach (var psm in psmTsvList)
            {
                if (psm.FullSequence.Contains("|"))
                {
                    continue;
                }
                var seq = new Ms2PipInput(psm);
                inputs.Add(seq);
            }
            var seqFile = new Ms2PipInputFile { Results = inputs };
            seqFile.WriteResults(outPath);
        }

        

        [Test]
        public static void TestCasanovoFileReading()
        {
            var path = @"E:\Aneuploidy\DDA\071525\casanovo_20250731231227.mztab";
            foreach (var line in File.ReadLines(path))
            {
                // Skip metadata lines (those starting with 'MTD', 'COM', etc.)
                if (line.StartsWith("MTD") || line.StartsWith("COM") || line.StartsWith("PSH") || string.IsNullOrWhiteSpace(line))
                    continue;

                // Only keep data table rows (e.g., PSM section starts with 'PSH' or 'PSM')
                if (line.StartsWith("PSM"))
                {
                    var columns = line.Split('\t');
                    var denovoResult = new CasanovoResult(line);
                }
            }
        }
    }
}
