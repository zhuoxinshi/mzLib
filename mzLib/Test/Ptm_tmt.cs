using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using System.Windows.Markup;
using Omics;
using System.Data.Entity.Core.Common.CommandTrees.ExpressionBuilder;
using System.IO;
using System.Text.RegularExpressions;
using System.Windows.Shapes;
using Easy.Common.Extensions;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using static UsefulProteomicsDatabases.ProteinDbRetriever;

namespace Test
{
    public class Ptm_tmt
    {
        [Test]
        public static void TMT_BestPsm()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMT = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var allPeptidesTMTWithMods = allPeptidesTMT.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain"));

            var model = new Prosit2020iRTTMT();
            var rtPredictionOutPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\TMT_BestPsm_RtPrediction.tsv";
            var rtColumns = new List<string> { "FullSequence", "Mod", "ObservedRT", "PredictedRT" };
            using (StreamWriter writer = new StreamWriter(rtPredictionOutPath))
            {
                writer.WriteLine(string.Join("\t", rtColumns));
                foreach (var peptide in allPeptidesTMT)
                {
                    var predictionInput = new List<RetentionTimePredictionInput> { new RetentionTimePredictionInput (peptide.FullSequence.Replace("X", "N-terminus")) };
                    if (predictionInput.First().SequenceWarning != null)
                    {
                        continue;
                    }
                    var predictedRt = model.Predict(predictionInput).First().PredictedRetentionTime;
                    var mod = "Unmodified";
                    if (!SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable")))
                    {
                        mod = "GptmdMod";
                        if (peptide.FullSequence.Contains("Biological") || peptide.FullSequence.Contains("Uniprot"))
                        {
                            mod = "BioMod";
                            if (peptide.FullSequence.Contains("Phospho"))
                            {
                                mod = "Phospho";
                            }
                        }
                    } 
                    var outString = new List<string> { peptide.FullSequence, mod, peptide.RetentionTime.ToString(), predictedRt.ToString() };
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
            //string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\ModifiedPeptidesTMT.tsv";
            //TmtPair.WriteResults(allPeptidesTMTWithMods, outPath);
        }

        [Test]
        public static void TMT_SumPsms()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPSMs.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMTWithMods = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && (SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain"))).GroupBy(p => p.FullSequence);

            string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\ModifiedPeptidesTMT_PSMs.tsv";
            //TmtPair.WriteResults(allPeptidesTMTWithMods, outPath);

        }

        [Test]
        public static void TransferIds()
        {
            string fullSequence = "AAGGAGAQVGGSISSGSSASSVTVTR";
            string converted = TmtPair.AddTmtLabelsToSequence(fullSequence);
            string fullSeq2 = "EVIITPNSAWGGEGSLGC[Common Fixed:Carbamidomethyl on C]GIGYGYLHR";
            string converted2 = TmtPair.AddTmtLabelsToSequence(fullSeq2);
            string fullSeq3 = "FGIVTSSAGTGTTEDTEAK[]";
            string converted3 = TmtPair.AddTmtLabelsToSequence(fullSeq3);
            string fullSeq4 = "HAEGTFTSDVSS[Common Biological:Phosphorylation on S]YLEGQAAK";
            string converted4 = TmtPair.AddTmtLabelsToSequence(fullSeq4);
            string removed1 = TmtPair.RemoveTmtLabels(converted);
            string removed2 = TmtPair.RemoveTmtLabels(converted2);
            string removed3 = TmtPair.RemoveTmtLabels(converted3);
            string removed4 = TmtPair.RemoveTmtLabels(converted4);

            var fullSeq = "[Multiplex Label:TMT18 on X]GQAGPEGAAP[Common Biological:Hydroxylation on P]APEEDK[Multiplex Label:TMT18 on K]";
            var parsed = SpectrumMatchFromTsv.ParseModifications(fullSeq);
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable" };
            var contains = parsed.Values.Any(v => !notInteresting.Any(key => v.Contains(key)));
        }

        [Test]
        public static void TMT_rt()
        {
            string allPeptidesTMT_noMod_path = @"E:\Islets\Brian_data\Real_islets\Frxn\F2-11_search-cali-search\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_noMod_file = new PsmFromTsvFile(allPeptidesTMT_noMod_path);
            var allPeptidesTMT_noMod = allPeptidesTMT_noMod_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && !p.FullSequence.Contains("|"));
            var observedRts = allPeptidesTMT_noMod.Select(p => p.RetentionTime).ToArray();
            var model = new Prosit2020iRTTMT();
            var inputs = allPeptidesTMT_noMod.Select(p => new RetentionTimePredictionInput(p.FullSequence.Replace("Multiplex Label:TMT18", "Common Fixed:TMTpro").Replace("X", "N-terminus"))).ToList();//.Replace("X", "N-terminus")
            var predictions = model.Predict(inputs).Select(p => p.PredictedRetentionTime).ToArray();

            var model_noTMT = new Prosit2019iRT();
            var input_noTMT = allPeptidesTMT_noMod.Select(p => new RetentionTimePredictionInput(TmtPair.RemoveTmtLabels(p.FullSequence))).ToList();
            var predictions_noTMT = model_noTMT.Predict(input_noTMT).Select(p => p.PredictedRetentionTime).ToArray();

            var outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\F2-11_search-cali-search\Task1-SearchTask\TMT_rt_predictions.tsv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", new List<string> { "PredictedRt_TMT", "PredictedRt_noTMT", "Diff", "ObservedRt"}));
                for (int i = 0; i < predictions.Length; i++)
                {
                    if (predictions[i] == null || predictions_noTMT[i] == null) continue;
                    writer.WriteLine(string.Join("\t", new List<string> { predictions[i].ToString(), predictions_noTMT[i].ToString(), (predictions[i] - predictions_noTMT[i]).ToString(), observedRts[i].ToString() }));
                }
            }
        }
    }

    public class TmtPair
    {
        public PsmFromTsv BestPsm_LF { get; set; }
        public List<PsmFromTsv> AllPsmsTMT { get; set; }
        public double[] ReporterIonIntensities { get; set; }

        public TmtPair(PsmFromTsv bestPsm_LF, IEnumerable<PsmFromTsv> psms = null)
        {
            BestPsm_LF = bestPsm_LF;
            AllPsmsTMT = psms == null ? new List<PsmFromTsv>() : psms.ToList();
        }

        public void FilterPsms(int fractionTol, double rtTol, double similarityTol)
        {
            string fraction_LF = Regex.Match(BestPsm_LF.FileNameWithoutExtension, @"Frxn\d+").Value;
            for (int i = AllPsmsTMT.Count - 1; i >= 0; i--)
            {
                var psm = AllPsmsTMT[i];
                var fraction_tmt = Regex.Match(psm.FileNameWithoutExtension, @"Frxn\d+").Value;
                if (Math.Abs(int.Parse(fraction_LF.Substring(4)) - int.Parse(fraction_tmt.Substring(4))) > fractionTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }

                //RT filter
                var predictedRt = 0;
                if (Math.Abs(predictedRt - psm.RetentionTime) > rtTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }

                //predicted spectrum
                var spectralSimilarity = 0;
                if (spectralSimilarity < similarityTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }
            }
        }

        public static void WriteResults(IEnumerable<PsmFromTsv> allPsms, string outPath)
        {
            var modsToexclude = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                var peptides = allPsms.GroupBy(p => p.FullSequence);
                foreach (var peptide in peptides)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.First().FullSequence).Values.Where(m => !modsToexclude.Any(x => m.Contains(x))));
                    var outString = new List<string> { peptide.First().BaseSequence, peptide.First().FullSequence, mods, peptide.First().ProteinAccession, peptide.First().ProteinName, peptide.First().GeneName, peptide.First().Description };

                    var reporterIonIntensities = new double[peptide.First().Intensities.Count()];
                    foreach (var psm in peptide)
                    {
                        for (int i = 0; i < peptide.First().Intensities.Count(); i++)
                        {
                            reporterIonIntensities[i] += psm.Intensities[i];
                        }
                    }
                    outString.AddRange(reporterIonIntensities.Select(i => i.ToString()));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        public static string AddTmtLabelsToSequence(string fullSequence)
        {
            var sb = new StringBuilder();

            if (fullSequence.StartsWith("["))
            {
                // Already modified N-term
                int endBracket = fullSequence.IndexOf(']');
                sb.Append(fullSequence.Substring(0, endBracket + 1));
                fullSequence = fullSequence.Substring(endBracket + 1);
            }
            else
            {
                // Add TMT to N-term
                sb.Append("[Multiplex Label:TMT18 on X]");
            }

            // Process the rest
            for (int i = 0; i < fullSequence.Length; i++)
            {
                char c = fullSequence[i];
                sb.Append(c);

                if (c == 'K')
                {
                    // Check if next char is a modification
                    if (i + 1 < fullSequence.Length && fullSequence[i + 1] == '[')
                    {
                        // Already modified, do nothing
                        continue;
                    }
                    else
                    {
                        sb.Append("[Multiplex Label:TMT18 on K]");
                    }
                }
            }

            return sb.ToString();
        }

        public static string RemoveTmtLabels(string fullSequence)
        {
            // Remove N-term TMT label if present
            fullSequence = Regex.Replace(
                fullSequence,
                @"^\[Multiplex Label:TMT18 on X\]",
                "",
                RegexOptions.Compiled);

            // Remove TMT label on K (after K)
            fullSequence = Regex.Replace(
                fullSequence,
                @"K\[Multiplex Label:TMT18 on K\]",
                "K",
                RegexOptions.Compiled);

            return fullSequence;
        }
    }
}
