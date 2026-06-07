using NUnit.Framework;
using Omics;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Proteomics.AminoAcidPolymer;
using Readers;
using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Test
{
    public class DeepLC
    {

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
        public static void TestDeepLCInput()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            string allPeptides_path = @"E:\Islets\Brian_data\PTM\LF_ptm\calied_avged_Gptmd_filters_secondPass\Task2-SearchTask\AllPeptides.psmtsv";
            var allPeptides_file = new PsmFromTsvFile(allPeptides_path);
            var allPeptides = allPeptides_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");
            //var allPeptidesTMTWithMods = allPeptidesTMT.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
            //var highScoreNoMod = allPeptides.Where(p => Ptm_tmt.ModType(p.FullSequence) == "Unmodified" && p.Score >= 20);

            var outPath = @"E:\Islets\Brian_data\PTM\Validations\LF_deepLC_input.csv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "seq", "modifications", "ObservedRt", "ModType" }));
                foreach (var peptide in allPeptides)
                {
                    var parsedMods = ParseModsForDeepLC(peptide.FullSequence, peptide.BaseSeq);
                    var modType = Ptm_tmt.ModType(peptide.FullSequence);
                    writer.WriteLine(string.Join(",", new List<string> { peptide.BaseSeq, parsedMods, Math.Round(peptide.RetentionTime, 2).ToString(), modType }));
                }
            }
        }

        [Test]
        public static void CalibrationPeptides()
        {
            string allPeptides_path = @"E:\Islets\Brian_data\PTM\LF_ptm\calied_avged_Gptmd_filters_secondPass\Task2-SearchTask\AllPeptides.psmtsv";
            var allPeptides_file = new PsmFromTsvFile(allPeptides_path);
            var allPeptides = allPeptides_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");
            var highScoreNoMod = allPeptides.Where(p => Ptm_tmt.ModType(p.FullSequence) == "Unmodified" && p.Score >= 20);

            var outPath = @"E:\Islets\Brian_data\PTM\Validations\LF_deepLC_highScoreNoMod_calibration.csv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "seq", "modifications", "ObservedRt", "ModType" }));
                foreach (var peptide in highScoreNoMod)
                {
                    var parsedMods = ParseModsForDeepLC(peptide.FullSequence, peptide.BaseSeq);
                    var modType = Ptm_tmt.ModType(peptide.FullSequence);
                    writer.WriteLine(string.Join(",", new List<string> { peptide.BaseSeq, parsedMods, Math.Round(peptide.RetentionTime, 2).ToString(), modType }));
                }
            }
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
                writer.WriteLine(string.Join("\t", new List<string> { "PredictedRt_TMT", "PredictedRt_noTMT", "Diff", "ObservedRt" }));
                for (int i = 0; i < predictions.Length; i++)
                {
                    if (predictions[i] == null || predictions_noTMT[i] == null) continue;
                    writer.WriteLine(string.Join("\t", new List<string> { predictions[i].ToString(), predictions_noTMT[i].ToString(), (predictions[i] - predictions_noTMT[i]).ToString(), observedRts[i].ToString() }));
                }
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
                writer.WriteLine(string.Join("\t", new List<string> { "LF_RT", "TMT_RT", "NumberOfTags", "Mod" }));
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

        public static string ApplyAaSubstitutions(string fullSequence)
        {
            Regex SubstitutionRx = new( @"([A-Z])\[\s*\d+\s+nucleotide\s+substitutions?\s*:\s*([A-Z])\s*->\s*([A-Z])\s+on\s+([A-Z])\s*\]", RegexOptions.Compiled);
            if (string.IsNullOrEmpty(fullSequence)) return fullSequence;

            return SubstitutionRx.Replace(fullSequence, m =>
            {
                char preceding = m.Groups[1].Value[0];
                char from = m.Groups[2].Value[0];
                char to = m.Groups[3].Value[0];
                char onWhich = m.Groups[4].Value[0];

                if (preceding != from || from != onWhich)
                    return m.Value;   // malformed -> leave as-is

                return to.ToString();
            });
        }

        [Test]
        public static void AAsubInput()
        {
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var peptideFile = new PsmFromTsvFile(peptidePath);
            var allPeptides = peptideFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(allPeptides);

            var outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\unmodPeptides_aasub_deepLC_input.csv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "seq", "modifications", "ObservedRt", "AAsub" }));
                foreach (var peptide in unmodPeptides)
                {
                    var aasub = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.Where(v => v.Contains("nucleotide substitution"));
                    var aasubString = aasub.Any() ? string.Join(";", aasub) : "";
                    var baseSeq = peptide.BaseSeq;
                    var fullSeq = peptide.FullSequence;
                    if (peptide.FullSequence.Contains("sub")) 
                    {
                        fullSeq = ApplyAaSubstitutions(peptide.FullSequence);
                        baseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSeq);
                    }
                    var parsedMods = ParseModsForDeepLC(fullSeq, baseSeq);
                    writer.WriteLine(string.Join(",", new List<string> { baseSeq, parsedMods, Math.Round(peptide.RetentionTime, 2).ToString(), aasubString }));
                }
            }
        }

        // Monoisotopic mass shifts of common PTMs that an AA substitution can mimic.
        // Used together with the < 2 Da blanket cutoff in HasPtmConfoundedAaSubstitution.
        private static readonly double[] PtmConfounderMasses =
        {
            14.0157,  // methylation
            15.9949,  // oxidation / hydroxylation
            42.0106,  // acetylation
            79.9663,  // phosphorylation
        };
        // ±tolerance around each PTM mass. 0.05 Da is conservative for high-res Orbitrap data.
        private const double PtmConfounderTolerance = 1;

        // Extracts the from/to amino acids out of e.g. "[1 nucleotide substitution:A->S on A]".
        private static readonly Regex AaSubExtractRx = new(
            @"\[\s*\d+\s+nucleotide\s+substitutions?\s*:\s*([A-Z])\s*->\s*([A-Z])\s+on\s+([A-Z])\s*\]",
            RegexOptions.Compiled);

        /// <summary>
        /// True if the peptide carries any nucleotide-substitution modification whose
        /// residue-mass delta is &lt; 2 Da (catches I/L isobars, Q/K ~0.04 Da, N/D ~0.98 Da)
        /// or within ±PtmConfounderTolerance of a common PTM shift (oxidation, acetylation,
        /// methylation/di/tri, phospho, etc.). Such substitutions are excluded because
        /// they cannot be confidently distinguished from a PTM at the precursor-mass level.
        /// </summary>
        /// <summary>
        /// Returns the per-substitution residue-mass shift (mass(to) - mass(from), signed)
        /// for every "[N nucleotide substitution:X->Y on X]" annotation in the FullSequence,
        /// in the order they appear. Empty if the peptide has no AA-sub mods.
        /// </summary>
        public static List<double> GetAaSubstitutionMassShifts(string fullSequence)
        {
            var shifts = new List<double>();
            if (string.IsNullOrEmpty(fullSequence)) return shifts;
            foreach (Match m in AaSubExtractRx.Matches(fullSequence))
            {
                char from = m.Groups[1].Value[0];
                char to   = m.Groups[2].Value[0];
                double mFrom = Residue.ResidueMonoisotopicMass[(int)from];
                double mTo   = Residue.ResidueMonoisotopicMass[(int)to];
                if (double.IsNaN(mFrom) || double.IsNaN(mTo)) continue;
                shifts.Add(mTo - mFrom);
            }
            return shifts;
        }

        /// <summary>
        /// Returns one annotation string per AA-substitution mod in the FullSequence, in the
        /// form "X -> Y: ±N.NNNN Da" where the mass shift is mass(Y) - mass(X). Empty if no
        /// AA-sub mods are present.
        /// </summary>
        public static List<string> GetAaSubstitutionAnnotations(string fullSequence)
        {
            var annotations = new List<string>();
            if (string.IsNullOrEmpty(fullSequence)) return annotations;
            foreach (Match m in AaSubExtractRx.Matches(fullSequence))
            {
                char from = m.Groups[1].Value[0];
                char to   = m.Groups[2].Value[0];
                double mFrom = Residue.ResidueMonoisotopicMass[(int)from];
                double mTo   = Residue.ResidueMonoisotopicMass[(int)to];
                if (double.IsNaN(mFrom) || double.IsNaN(mTo)) continue;
                double delta = mTo - mFrom;
                annotations.Add(string.Format(CultureInfo.InvariantCulture,
                    "{0} -> {1}: {2:F4} Da", from, to, delta));
            }
            return annotations;
        }

        public static bool HasPtmConfoundedAaSubstitution(string fullSequence)
        {
            if (string.IsNullOrEmpty(fullSequence)) return false;
            foreach (Match m in AaSubExtractRx.Matches(fullSequence))
            {
                char from = m.Groups[1].Value[0];
                char to   = m.Groups[2].Value[0];
                double mFrom = Residue.ResidueMonoisotopicMass[(int)from];
                double mTo   = Residue.ResidueMonoisotopicMass[(int)to];
                if (double.IsNaN(mFrom) || double.IsNaN(mTo)) continue;
                double delta = Math.Abs(mTo - mFrom);
                if (delta < 2.0) return true;
                for (int i = 0; i < PtmConfounderMasses.Length; i++)
                    if (Math.Abs(delta - PtmConfounderMasses[i]) <= PtmConfounderTolerance) return true;
            }
            return false;
        }

        // -----------------------------------------------------------------------
        // FragPipe Observed Modifications -> MetaMorpheus FullSequence
        //
        // We trust FragPipe's own classification in the "Observed Modifications"
        // column. Each entry is "Name(mass)" e.g.:
        //   - Oxidation(15.994900)
        //   - Carbamidomethyl(57.021400)
        //   - Val->Trp(87.010900)         <- AA substitution, 3-letter codes
        //   - Carbamyl(43.005800), Formyl(27.994900), Dioxidation(...), ...
        //
        // Rules for this test:
        //   - Acceptable mods: Oxidation, Carbamidomethyl, AA->AA substitution
        //   - Substitutions are dropped when the residue-mass delta is < 2 Da
        //     or within ±0.05 Da of a known PTM (PtmConfounderMasses).
        //   - Any other observed mod (Carbamyl/Formyl/Dioxidation/Carboxymethyl/etc.)
        //     causes the peptide to be dropped.
        //   - CAM is added to every Cys regardless of whether the column lists it
        //     (it's a universal fixed mod here).
        //   - If Oxidation is listed, every Met in the peptide is annotated as
        //     oxidised. Single-Met peptides are unambiguous; for multi-Met
        //     peptides this slightly over-annotates but is simple and conservative.
        // -----------------------------------------------------------------------

        private const string CamAnnotation = "[Common Fixed:Carbamidomethyl on C]";
        private const string OxAnnotation  = "[Common Variable:Oxidation on M]";

        // Three-letter -> one-letter AA code (for parsing "Val->Trp"-style names).
        private static readonly Dictionary<string, char> ThreeToOne = new(StringComparer.OrdinalIgnoreCase)
        {
            {"Ala",'A'},{"Arg",'R'},{"Asn",'N'},{"Asp",'D'},{"Cys",'C'},
            {"Glu",'E'},{"Gln",'Q'},{"Gly",'G'},{"His",'H'},{"Ile",'I'},
            {"Leu",'L'},{"Lys",'K'},{"Met",'M'},{"Phe",'F'},{"Pro",'P'},
            {"Ser",'S'},{"Thr",'T'},{"Trp",'W'},{"Tyr",'Y'},{"Val",'V'}
        };

        // Matches a single FragPipe Observed Modification entry: "Name(mass)".
        // Name can include `-`, `>` (for substitutions) and word chars.
        private static readonly Regex ObservedModRx = new(
            @"([A-Za-z][\w>\-]*)\s*\(\s*([-+]?\d+(?:\.\d+)?)\s*\)",
            RegexOptions.Compiled);

        // "Val->Trp" -> ('V', 'W'). Null if not a substitution.
        private static readonly Regex SubNameRx = new(
            @"^([A-Z][a-z]{2})->([A-Z][a-z]{2})$",
            RegexOptions.Compiled);

        public static List<(string Name, double Mass)> ParseObservedMods(string field)
        {
            var result = new List<(string, double)>();
            if (string.IsNullOrWhiteSpace(field)) return result;
            foreach (Match m in ObservedModRx.Matches(field))
            {
                if (!double.TryParse(m.Groups[2].Value, NumberStyles.Float, CultureInfo.InvariantCulture, out double mass))
                    continue;
                result.Add((m.Groups[1].Value.Trim(), mass));
            }
            return result;
        }

        public static (char From, char To)? TryParseAaSubstitutionName(string name)
        {
            var m = SubNameRx.Match(name);
            if (!m.Success) return null;
            if (!ThreeToOne.TryGetValue(m.Groups[1].Value, out var from)) return null;
            if (!ThreeToOne.TryGetValue(m.Groups[2].Value, out var to))   return null;
            return (from, to);
        }

        // Build {Peptide -> median RT} from FragPipe psm.tsv (target-only).
        private static Dictionary<string, double> LoadFragPipePsmMedianRt(string psmTsvPath)
        {
            var perPep = new Dictionary<string, List<double>>(capacity: 50_000);
            using (var sr = new StreamReader(psmTsvPath))
            {
                var firstLine = sr.ReadLine();
                if (firstLine == null) return new Dictionary<string, double>();
                var header = firstLine.Split('\t');
                int idxPep    = Array.IndexOf(header, "Peptide");
                int idxRt     = Array.IndexOf(header, "Retention");
                int idxDecoy  = Array.IndexOf(header, "Is Decoy");
                int idxContam = Array.IndexOf(header, "Is Contaminant");
                if (idxPep < 0 || idxRt < 0)
                    throw new Exception($"psm.tsv missing required columns (Peptide / Retention): {psmTsvPath}");

                int required = new[] { idxPep, idxRt, idxDecoy, idxContam }.Where(i => i >= 0).Max() + 1;

                string line;
                while ((line = sr.ReadLine()) != null)
                {
                    var parts = line.Split('\t');
                    if (parts.Length < required) continue;
                    if (idxDecoy  >= 0 && parts[idxDecoy].Equals("true",  StringComparison.OrdinalIgnoreCase)) continue;
                    if (idxContam >= 0 && parts[idxContam].Equals("true", StringComparison.OrdinalIgnoreCase)) continue;
                    if (!double.TryParse(parts[idxRt], NumberStyles.Float, CultureInfo.InvariantCulture, out double rt)) continue;
                    string pep = parts[idxPep];
                    if (!perPep.TryGetValue(pep, out var list))
                        perPep[pep] = list = new List<double>();
                    list.Add(rt);
                }
            }
            var result = new Dictionary<string, double>(perPep.Count);
            foreach (var kv in perPep)
            {
                var arr = kv.Value.ToArray();
                Array.Sort(arr);
                int n = arr.Length;
                result[kv.Key] = (n % 2 == 1) ? arr[n / 2] : 0.5 * (arr[n / 2 - 1] + arr[n / 2]);
            }
            return result;
        }

        [Test]
        public static void FragPipeAAsubInput_Prosit2019Predict()
        {
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\FragPipe\AAsub\peptide.tsv";
            var psmPath     = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\FragPipe\AAsub\psm.tsv";
            var outPath     = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\FragPipe\AAsub\fragpipe_unmodPeptides_aasub_Prosit2019_predictions2.csv";

            // Median RT per peptide from psm.tsv (target-only).
            var rtByPeptide = LoadFragPipePsmMedianRt(psmPath);

            var lines = File.ReadAllLines(peptidePath);
            var header = lines[0].Split('\t');
            int idxPep      = Array.IndexOf(header, "Peptide");
            int idxObserved = Array.IndexOf(header, "Observed Modifications");
            int idxQ        = Array.IndexOf(header, "Qvalue");
            int idxDecoy    = Array.IndexOf(header, "Is Decoy");
            int idxContam   = Array.IndexOf(header, "Is Contaminant");
            if (idxPep < 0 || idxObserved < 0 || idxQ < 0 || idxDecoy < 0 || idxContam < 0)
                throw new Exception($"peptide.tsv missing required columns: {peptidePath}");

            int required = new[] { idxPep, idxObserved, idxQ, idxDecoy, idxContam }.Max() + 1;

            var keepers = new List<(string FullSeqUsed, string AasubAnnot, int NumSubs, double Rt)>();

            for (int li = 1; li < lines.Length; li++)
            {
                if (string.IsNullOrWhiteSpace(lines[li])) continue;
                var parts = lines[li].Split('\t');
                if (parts.Length < required) continue;
                if (parts[idxDecoy].Equals("true",  StringComparison.OrdinalIgnoreCase)) continue;
                if (parts[idxContam].Equals("true", StringComparison.OrdinalIgnoreCase)) continue;
                if (!double.TryParse(parts[idxQ], NumberStyles.Float, CultureInfo.InvariantCulture, out double q) || q > 0.01) continue;

                // 1. Parse Observed Modifications and classify by FragPipe's labels.
                bool drop = false;
                bool hasOxidation = false;
                var subsForLabel = new List<(char From, char To)>();
                var seenSubKeys = new HashSet<string>();

                foreach (var (name, mass) in ParseObservedMods(parts[idxObserved]))
                {
                    if (name.Equals("Oxidation",       StringComparison.OrdinalIgnoreCase)) { hasOxidation = true; continue; }
                    if (name.Equals("Carbamidomethyl", StringComparison.OrdinalIgnoreCase)) { continue; } // fixed; always added to every C

                    var sub = TryParseAaSubstitutionName(name);
                    if (sub.HasValue)
                    {
                        // Same PTM-confound filter as the MetaMorpheus test: drop the peptide if
                        // the residue-mass delta is < 2 Da or within ±0.05 Da of a common PTM.
                        var (from, to) = sub.Value;
                        double dm = Math.Abs(Residue.ResidueMonoisotopicMass[(int)to] - Residue.ResidueMonoisotopicMass[(int)from]);
                        if (dm < 2.0) { drop = true; break; }
                        if (PtmConfounderMasses.Any(p => Math.Abs(dm - p) <= PtmConfounderTolerance)) { drop = true; break; }

                        if (seenSubKeys.Add($"{from}->{to}"))
                            subsForLabel.Add((from, to));
                        continue;
                    }

                    // Any other observed mod name (Carbamyl/Formyl/Dioxidation/Carboxymethyl/…) -> drop
                    drop = true;
                    break;
                }
                if (drop) continue;

                // 2. Build the MM FullSequence from FragPipe's Peptide column.
                //    CAM annotated on every Cys; Oxidation annotated on every Met when
                //    "Oxidation" appears in the Observed column. We do NOT try to apply
                //    the AA substitution to the sequence — the Peptide column is taken
                //    at face value, and the substitution name is reported as-is in the
                //    AAsub output column.
                string baseSeq = parts[idxPep];
                if (baseSeq.Length < 1 || baseSeq.Length > 30) continue;
                if (!Regex.IsMatch(baseSeq, "^[ACDEFGHIKLMNPQRSTVWY]+$")) continue;

                var sb = new StringBuilder(baseSeq.Length + 64);
                foreach (char c in baseSeq)
                {
                    sb.Append(c);
                    if (c == 'C')                     sb.Append(CamAnnotation);
                    else if (c == 'M' && hasOxidation) sb.Append(OxAnnotation);
                }
                string mmFullSeq = sb.ToString();

                if (!rtByPeptide.TryGetValue(baseSeq, out double rt)) continue;

                string aasubAnnot = subsForLabel.Count == 0
                    ? ""
                    : string.Join("; ", subsForLabel.Select(s => $"{s.From}->{s.To}"));

                keepers.Add((mmFullSeq, aasubAnnot, subsForLabel.Count, rt));
            }

            // 4. Predict with Prosit2019iRT in small chunks so a single bad peptide doesn't
            //    poison a whole batch. We keep the full PeptideRTPrediction (not just the
            //    RT) so we can also surface the Warning string in the CSV.
            var model = new Prosit2019iRT();
            const int predictChunkSize = 100;
            var predictionsList = new List<PredictionClients.Koina.AbstractClasses.PeptideRTPrediction>(keepers.Count);
            for (int start = 0; start < keepers.Count; start += predictChunkSize)
            {
                int take = Math.Min(predictChunkSize, keepers.Count - start);
                var chunkInputs = new List<RetentionTimePredictionInput>(take);
                for (int i = 0; i < take; i++)
                    chunkInputs.Add(new RetentionTimePredictionInput(keepers[start + i].FullSeqUsed));
                try
                {
                    predictionsList.AddRange(model.Predict(chunkInputs));
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"[Prosit2019iRT] chunk [{start}, {start + take}) THREW: {ex.Message}");
                    Console.WriteLine($"  first peptide in chunk: {keepers[start].FullSeqUsed}");
                    for (int i = 0; i < take; i++)
                        predictionsList.Add(new PredictionClients.Koina.AbstractClasses.PeptideRTPrediction(
                            FullSequence: keepers[start + i].FullSeqUsed,
                            PredictedRetentionTime: null,
                            IsIndexed: null) );
                }
            }

            // Summary so it's obvious what happened
            int nWithRt = predictionsList.Count(p => p.PredictedRetentionTime.HasValue);
            int nNull   = predictionsList.Count(p => !p.PredictedRetentionTime.HasValue);
            Console.WriteLine($"[FragPipe] predictions returned: {predictionsList.Count}  (with RT: {nWithRt}, null: {nNull})");
            if (nWithRt == 0 && predictionsList.Count > 0)
            {
                // Show the first few warnings so we can see WHY nothing predicted.
                var sampleWarnings = predictionsList
                    .Select(p => p.Warning?.Message ?? "(no warning)")
                    .Distinct().Take(5);
                Console.WriteLine("[FragPipe] every prediction is null. Distinct warnings (up to 5):");
                foreach (var w in sampleWarnings) Console.WriteLine($"  - {w}");
                Console.WriteLine("[FragPipe] first 3 sequences sent to Prosit:");
                foreach (var k in keepers.Take(3)) Console.WriteLine($"  - {k.FullSeqUsed}");
            }

            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "BaseSeq", "FullSequenceUsed", "PredictedRt", "ObservedRt", "AAsub", "n_subs", "PredictWarning" }));
                int n = Math.Min(keepers.Count, predictionsList.Count);
                for (int i = 0; i < n; i++)
                {
                    var k = keepers[i];
                    var p = predictionsList[i];
                    var baseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(k.FullSeqUsed);
                    var predStr = p.PredictedRetentionTime.HasValue
                        ? p.PredictedRetentionTime.Value.ToString(CultureInfo.InvariantCulture)
                        : "";
                    var warnStr = p.Warning?.Message?.Replace(',', ';') ?? "";
                    writer.WriteLine(string.Join(",", new List<string>
                    {
                        baseSeq,
                        k.FullSeqUsed,
                        predStr,
                        Math.Round(k.Rt, 4).ToString(CultureInfo.InvariantCulture),
                        k.AasubAnnot,
                        k.NumSubs.ToString(CultureInfo.InvariantCulture),
                        warnStr
                    }));
                }
            }
        }

        [Test]
        public static void AAsubInput_Prosit2019Predict()
        {
            // Mirrors AAsubInput: same source file, same filtering / unmod selection.
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var peptideFile = new PsmFromTsvFile(peptidePath);
            var allPeptides = peptideFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            // Drop peptides whose AA substitutions could be confused with a PTM
            // (Δmass < 2 Da, or within ±0.05 Da of oxidation/acetylation/methylation/phospho/etc.).
            var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(allPeptides)
                .Where(p => !HasPtmConfoundedAaSubstitution(p.FullSequence))
                .ToList();

            // Build the sequence we feed to Prosit:
            // - if the peptide carries a "nucleotide substitution" mod, apply the
            //   substitution so the model sees the actual translated peptide;
            // - otherwise pass the FullSequence through unchanged.
            var sequencesForPrediction = unmodPeptides
                .Select(p => p.FullSequence.Contains("sub") ? ApplyAaSubstitutions(p.FullSequence) : p.FullSequence)
                .ToList();

            var model = new Prosit2019iRT();
            var inputs = sequencesForPrediction.Select(s => new RetentionTimePredictionInput(s)).ToList();
            var predictions = model.Predict(inputs).Select(p => p.PredictedRetentionTime).ToArray();

            var outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\unmodPeptides_aasub_Prosit2019_predictions_filtered.csv";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join(",", new List<string> { "BaseSeq", "FullSequenceUsed", "PredictedRt", "ObservedRt", "AAsub" }));
                for (int i = 0; i < unmodPeptides.Count; i++)
                {
                    var peptide = unmodPeptides[i];
                    var aasub = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.Where(v => v.Contains("nucleotide substitution"));
                    var aasubString = aasub.Any() ? string.Join(";", aasub) : "";
                    var fullSeqUsed = sequencesForPrediction[i];
                    var baseSeq = peptide.FullSequence.Contains("sub")
                        ? IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSeqUsed)
                        : peptide.BaseSeq;
                    var predRtStr = predictions[i] == null ? "" : predictions[i].ToString();
                    writer.WriteLine(string.Join(",", new List<string> { baseSeq, fullSeqUsed, predRtStr, Math.Round(peptide.RetentionTime, 2).ToString(), aasubString }));
                }
            }
        }

    }
}
