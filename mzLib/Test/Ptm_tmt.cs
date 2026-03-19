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

namespace Test
{
    public class Ptm_tmt
    {

        [Test]
        public static void FindIDs()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesLF_path = @"E:\Islets\Brian_data\LF_ptm\gptmd_secondPass\Task2-SearchTask\AllPeptides.psmtsv";
            var allPeptidesLF_file = new PsmFromTsvFile(allPeptidesLF_path);
            var allPeptidesLFWithMods = allPeptidesLF_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && (SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)) || p.Description.Contains("chain")))).Take(20).ToList();

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPSMs.psmtsv";
            var allPsmTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPsms = allPsmTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T").ToList();
            //var allPsmTMTWithMods = allPsmTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)))).ToList();

            var tmtPairs = new List<TmtPair>();
            foreach (var peptide in allPeptidesLFWithMods)
            {
                string converted = TmtPair.ConvertSequence(peptide.FullSequence);
                var matchTMT = allPsms.Where(p => p.FullSequence == converted);
                if (matchTMT.Count() > 0)
                {
                    var newPair = new TmtPair(peptide, matchTMT);
                    tmtPairs.Add(newPair);
                }
            }
            tmtPairs.Sort((a, b) => b.AllPsmsTMT.Count().CompareTo(a.AllPsmsTMT.Count()));

            string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\TMT_Pairs.tsv";
            TmtPair.WriteResults(tmtPairs, outPath);
        }

        [Test]
        public static void FindIDsTMTonly()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMTWithMods = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && (SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain")));

            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\ModifiedPeptidesTMT.tsv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                foreach (var peptide in allPeptidesTMTWithMods)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.Where(m => !notInteresting.Any(x => m.Contains(x))));
                    var outString = new List<string> { peptide.BaseSequence, peptide.FullSequence, mods, peptide.ProteinAccession, peptide.ProteinName, peptide.GeneName, peptide.Description };
                    outString.AddRange(peptide.ReporterIonIntensities.Select(i => i.ToString()));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        [Test]
        public static void FindIDsTMTPsms()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\Task1-SearchTask\AllPSMs.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMTWithMods = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && (SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain"))).GroupBy(p => p.FullSequence);

            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\ModifiedPeptidesTMT_PSMs.tsv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                foreach (var peptide in allPeptidesTMTWithMods)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.First().FullSequence).Values.Where(m => !notInteresting.Any(x => m.Contains(x))));
                    var outString = new List<string> { peptide.First().BaseSequence, peptide.First().FullSequence, mods, peptide.First().ProteinAccession, peptide.First().ProteinName, peptide.First().GeneName, peptide.First().Description };

                    var reporterIonIntensities = new double[18];
                    foreach (var psm in peptide)
                    {
                        for (int i = 0; i < 18; i++)
                        {
                            reporterIonIntensities[i] += psm.ReporterIonIntensities[i];
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
            string converted = TmtPair.ConvertSequence(fullSequence);
            string fullSeq2 = "EVIITPNSAWGGEGSLGC[Common Fixed:Carbamidomethyl on C]GIGYGYLHR";
            string converted2 = TmtPair.ConvertSequence(fullSeq2);
            string fullSeq3 = "FGIVTSSAGTGTTEDTEAK[]";
            string converted3 = TmtPair.ConvertSequence(fullSeq3);
            string fullSeq4 = "HAEGTFTSDVSS[Common Biological:Phosphorylation on S]YLEGQAAK";
            string converted4 = TmtPair.ConvertSequence(fullSeq4);

            var fullSeq = "[Multiplex Label:TMT18 on X]GQAGPEGAAP[Common Biological:Hydroxylation on P]APEEDK[Multiplex Label:TMT18 on K]";
            var parsed = SpectrumMatchFromTsv.ParseModifications(fullSeq);
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable" };
            var contains = parsed.Values.Any(v => !notInteresting.Any(key => v.Contains(key)));
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
            AllPsmsTMT = psms == null? new List<PsmFromTsv>() : psms.ToList();
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

        public void AggregateReporterIonIntensities()
        {
            ReporterIonIntensities = new double[18];
            foreach (var psm in AllPsmsTMT)
            {
                for (int i = 0; i < 18; i++)
                {
                    ReporterIonIntensities[i] += psm.ReporterIonIntensities[i];
                }
            }
        }
        public static void WriteResults(IEnumerable<TmtPair> pairs, string outPath)
        {
            var columns = new List<string> { "BaseSequence", "FullSequence", "Description", "Protein Accession", "Protein Name", "Fraction_LF", "Fraction_TMT", "RT_LF", "RT_TMT", "PsmCount_TMT" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));   
                foreach (var pair in pairs)
                {
                    string fraction_LF = Regex.Match(pair.BestPsm_LF.FileNameWithoutExtension, @"Frxn\d+").Value;
                    string fraction_TMT = String.Join(",", pair.AllPsmsTMT.Select(p => Regex.Match(p.FileNameWithoutExtension, @"Frxn\d+").Value));
                    string RT_TMT = String.Join(",",pair.AllPsmsTMT.Select(p => Math.Round(p.RetentionTime, 2).ToString()));
                    var outString = new List<string> { pair.BestPsm_LF.BaseSequence, pair.BestPsm_LF.FullSequence, pair.BestPsm_LF.Description, pair.BestPsm_LF.ProteinAccession, pair.BestPsm_LF.ProteinName, fraction_LF, fraction_TMT, pair.BestPsm_LF.RetentionTime.ToString(), RT_TMT, pair.AllPsmsTMT.Count().ToString() };
                    pair.AggregateReporterIonIntensities();
                    outString.AddRange(pair.ReporterIonIntensities.Select(i => i.ToString()));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        public static string ConvertSequence(string fullSequence)
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
    }
}
