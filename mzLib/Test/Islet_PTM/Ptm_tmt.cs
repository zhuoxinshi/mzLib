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

namespace Test
{
    public class Ptm_tmt
    {
        [Test]
        public static void HighResTMT()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            var allPeptidesHigh_path = @"E:\Islets\Brian_data\HighResTMT\noCali_LFgptmdFilterPruned\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesHigh_file = new PsmFromTsvFile(allPeptidesHigh_path);
            var allPeptidesHigh = allPeptidesHigh_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var allPeptidesNoModHigh = allPeptidesHigh.Where(p => !SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")) && !p.FullSequence.Contains("|"));
            var allPeptidesWithModsHigh = allPeptidesHigh.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain")).ToList();

            var allPeptidesLow_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesLow_file = new PsmFromTsvFile(allPeptidesLow_path);
            var allPeptidesLow = allPeptidesLow_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            //var allPeptidesWithModsLow = allPeptidesLow.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain")).ToList();

            var rtPairs = new List<(double high, double low)>();
            var highResOverlap = new List<PsmFromTsv>();
            foreach (var peptide in allPeptidesWithModsHigh)
            {
                var lowResPeptide = allPeptidesLow.FirstOrDefault(p => p.FullSequence == peptide.FullSequence);
                if (lowResPeptide != null)
                {
                    rtPairs.Add((peptide.RetentionTime, lowResPeptide.RetentionTime));
                    highResOverlap.Add(lowResPeptide);
                }
            }

            var lowResOutPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_LFgptmdFilterPrunedDb\ModifiedPeptidesTMT_HighResOverlap_BestPsm.tsv";
            TmtPair.WriteResults(highResOverlap, lowResOutPath);
        }

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

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPSMs.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPsms = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");

            string outPath = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Peptides_PsmSum.tsv";
            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                var peptides = allPsms.GroupBy(p => p.FullSequence);
                foreach (var peptide in peptides)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.First().FullSequence).Values.Where(m => !notInteresting.Any(x => m.Contains(x))));
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

            var outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\F2-11_search-cali-search\Task1-SearchTask\TMT_deepLC.csv";
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


        public static string ModType(string fullSequence)
        {
            var mod = "Unmodified";
            if (!SpectrumMatchFromTsv.ParseModifications(fullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")))
            {
                mod = "GptmdMod";
                if (fullSequence.Contains("Biological") || fullSequence.Contains("Uniprot"))
                {
                    mod = "BioMod";
                    if (fullSequence.IndexOf("Phospho", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        mod = "Phospho";
                    }
                }
            }
            return mod;
        }

        [Test]
        public static void TestMs2PIP()
        {
            string fullSeq = "[Multiplex Label:TMT18 on X]GQAGPEGAAP[Common Biological:Hydroxylation on P]APEEDK[Multiplex Label:TMT18 on K]";
            var parsed = ParseSequenceForMs2PIP(fullSeq);
        }

        [Test]
        public static void TestDeepLC()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMT = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.MissedCleavage == "0");
            //var allPeptidesTMTWithMods = allPeptidesTMT.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
            var highScoreNoMod = allPeptidesTMT.Where(p => ModType(p.FullSequence) =="Unmodified" && p.Score >= 20);

            var outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\DeepLC_highScoreNoMod_input.csv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "seq", "modifications", "ObservedRt", "ModType" }));
                foreach (var peptide in highScoreNoMod)
                {
                    var parsedMods = ParseModsForDeepLC(peptide.FullSequence, peptide.BaseSeq);
                    var modType = ModType(peptide.FullSequence);
                    writer.WriteLine(string.Join(",", new List<string> { peptide.BaseSeq, parsedMods, Math.Round(peptide.RetentionTime, 2).ToString(), modType }));
                }
            }
        }

        public static string ParseModsForDeepLC(string fullSequence, string baseSequence = null)
        {
            var mods = SpectrumMatchFromTsv.ParseModifications(fullSequence);
            if (mods.Count == 0) return "";
            var sb = new StringBuilder();
            foreach (var mod in mods)
            {
                var modName = mod.Value.Split(':')[1].Split(' ')[0];
                var parsedModName = ParseModNameForDeepLC(modName);
                if (parsedModName == null) return null;
                var modPosition = mod.Key;
                //if (modPosition == baseSequence.Length)
                //{
                //    modPosition = -1; 
                //}
                sb.Append($"{modPosition}|{parsedModName}|");
            }
            return sb.ToString().TrimEnd('|');
        }

        public static string ParseSequenceForMs2PIP(string fullSequence )
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
                    var parsedModName = ParseModNameForDeepLC(modName);
                    sb.Append($"[{parsedModName}]");
                }
            }
            return sb.ToString();
        }

        public static IEnumerable<PsmFromTsv> GetUnmodifiedPeptides(IEnumerable<PsmFromTsv> psms)
        {
            return psms.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")) && !p.FullSequence.Contains("|"));
        }

        public static IEnumerable<PsmFromTsv> GetModifiedPeptides(IEnumerable<PsmFromTsv> psms)
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            return psms.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
        }

        public static string ParseModNameForDeepLC(string modName)
        {
            // C#
            switch (modName)
            {
                case string s when s.Contains("TMT18"):
                    return "TMTpro";
                case string s when s.Contains("Carbamidomethyl"):
                    return "Carbamidomethyl";
                case string s when s.IndexOf("Phospho", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Phospho";
                case string s when s.IndexOf("Hydroxy", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Hydroxylation";
                case string s when s.IndexOf("Acetyl", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Acetyl";
                case string s when s.Contains("Citrullination"):
                    return "Deamidated";
                case string s when s.IndexOf("Sodium", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Sodiated";
                case string s when s.IndexOf("Methyl", StringComparison.OrdinalIgnoreCase) >= 0 && s.Contains("-"):
                    return "Methyl";
                case string s when s.IndexOf("dimethyl", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Dimethyl";
                case string s when s.IndexOf("trimethyl", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Trimethyl";
                case string s when s.IndexOf("Succinyl", StringComparison.OrdinalIgnoreCase) >= 0:
                    return "Succinyl";
                default:
                    return null;
            }

        }

        [Test]
        public static void unlabeledNoModRts()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            var lf_path = @"E:\Islets\Brian_data\LF_ptm\calied_avged_Gptmd_filters_secondPass\Task2-SearchTask\AllPeptides.psmtsv";
            var lf_file = new PsmFromTsvFile(lf_path);
            var lf_peptides = lf_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T").ToList();
            var lf_noMod = lf_peptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable")) && !p.FullSequence.Contains("|")).ToList();
            var lf_mod = lf_peptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)))).ToList();

            var tmtHigh_path = @"E:\Islets\Brian_data\HighResTMT\noCali_LFgptmdFilterPruned\Task1-SearchTask\AllPeptides.psmtsv";
            var tmtHigh_file = new PsmFromTsvFile(tmtHigh_path);
            var tmtHigh_peptides = tmtHigh_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T").ToList();
            var tmtHigh_noMod = tmtHigh_peptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")) && !p.FullSequence.Contains("|")).ToList();
            var tmtHigh_mod = tmtHigh_peptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)))).ToList();

            var outPath = @"E:\Islets\Brian_data\HighResTMT\noCali_LFgptmdFilterPruned\RT_comparisons.tsv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", new List<string> { "LF_RT", "TMT_RT", "NumberOfTags", "Mod"}));
                foreach (var lfPeptide in lf_noMod)
                {
                    var tmtHighPeptide = tmtHigh_noMod.FirstOrDefault(p => p.FullSequence == TmtPair.AddTmtLabelsToSequence(lfPeptide.FullSequence));
                    if (tmtHighPeptide != null)
                    {
                        var outString = new List<string> { lfPeptide.RetentionTime.ToString(), tmtHighPeptide.RetentionTime.ToString(), SpectrumMatchFromTsv.ParseModifications(tmtHighPeptide.FullSequence).Values.Count(v => v.Contains("TMT")).ToString(), "noMod" };
                        writer.WriteLine(string.Join("\t", outString));
                    }
                }
                foreach (var lfPeptide in lf_mod)
                {
                    var tmtHighPeptide = tmtHigh_mod.FirstOrDefault(p => p.FullSequence == TmtPair.AddTmtLabelsToSequence(lfPeptide.FullSequence));
                    if (tmtHighPeptide != null)
                    {
                        var outString = new List<string> { lfPeptide.RetentionTime.ToString(), tmtHighPeptide.RetentionTime.ToString(), SpectrumMatchFromTsv.ParseModifications(tmtHighPeptide.FullSequence).Values.Count(v => v.Contains("TMT")).ToString(), "Mod" };
                        writer.WriteLine(string.Join("\t", outString));
                    }
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
