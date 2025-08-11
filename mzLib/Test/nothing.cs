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
using CsvHelper.Configuration;
using CsvHelper;
using Omics.Modifications;
using Proteomics;
using System.Globalization;
using UsefulProteomicsDatabases;
using MzLibUtil;
using TopDownProteomics;
using Proteomics.ProteolyticDigestion;
using Easy.Common.Extensions;


namespace Test
{
    public class nothing
    {
        public static List<string> ModifiedSequences { get; set; }
        public static List<string> OriginalSequences { get; set; }
        public static List<Protein> OriginalProteins { get; set; }
        public static List<Protein> ModifiedProteins { get; set; } 
        public static List<PsmFromTsv> SelectedPeptides { get; set; }
        public static List<Protein> AllProteinsAfterMutation { get; set; }

        [OneTimeSetUp]
        public void SetUp()
        {
            var xmlDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var fastaDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            //var allProteins = ProteinDbLoader.LoadProteinXML(xmlDb, true, DecoyType.None, null, false, new List<string>(), out Dictionary<string, Modification> un);
            OriginalProteins = ProteinDbLoader.LoadProteinFasta(fastaDb, true, DecoyType.None, false, out List<string> errors);
            var psmPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\fasta\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            //var allPeptidesNoModFull = allPeptides.Where(p => !p.FullSequence.Contains("[") && p.Description == "full").ToList();
            var allPeptidesNoModFull = allPeptides.Where(p => !p.FullSequence.Contains("[") && p.Description == "full" && !p.FullSequence.Contains("|")).ToList();

            Random rand = new Random(42);
            int numberOfPeptideSeqToMutate = 100;
            int[] randomIndices = Enumerable.Range(0, allPeptidesNoModFull.Count).OrderBy(_ => rand.Next()).Take(numberOfPeptideSeqToMutate).ToArray();
            SelectedPeptides = randomIndices.Select(i => allPeptidesNoModFull[i]).ToList();

            ModifiedProteins = new List<Protein>();
            ModifiedSequences = new List<string>();
            OriginalSequences = new List<string>();
            AllProteinsAfterMutation = new List<Protein>(OriginalProteins);
            foreach (var peptide in SelectedPeptides)
            {
                var mutatedSequence = RandomlyMutateAminoAcid(peptide.BaseSeq);
                var proteins = AllProteinsAfterMutation.Where(p => p.BaseSequence.Contains(peptide.BaseSeq)).ToList();
                if (proteins.IsNullOrEmpty()) continue;
                foreach(var protein in proteins)
                {
                    int index = protein.BaseSequence.IndexOf(peptide.BaseSeq);
                    string newProteinSequence = null;
                    if (index >= 0)
                    {
                        newProteinSequence = protein.BaseSequence.Substring(0, index) + mutatedSequence + protein.BaseSequence.Substring(index + mutatedSequence.Length);
                    }
                    var modifiedProtein = new Protein(protein, newProteinSequence);
                    AllProteinsAfterMutation.Remove(protein);
                    AllProteinsAfterMutation.Add(modifiedProtein);
                    ModifiedProteins.Add(modifiedProtein);
                    ModifiedSequences.Add(mutatedSequence);
                    OriginalSequences.Add(peptide.BaseSeq);
                    if (newProteinSequence.Length != protein.BaseSequence.Length || newProteinSequence == protein.BaseSequence)
                    {
                        throw new MzLibException("Error in modifying protein sequence");
                    }
                }
            }
        }

        [Test]
        public void WriteNewFasta()
        {
            var outPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\modifiedYeastFasta100_noAmb.fasta";
            ProteinDbWriter.WriteFastaDatabase(AllProteinsAfterMutation, outPath, " ");
        }

        [Test]
        public static void WriteNewFastaOnePep()
        {
            var xmlDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var fastaDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            //var allProteins = ProteinDbLoader.LoadProteinXML(fastaDb, true, DecoyType.None, null, false,
            //    new List<string>(), out Dictionary<string, Modification> un);
            var allProteins = ProteinDbLoader.LoadProteinFasta(fastaDb, true, DecoyType.None, false, out List<string> errors);
            var psmPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\fasta\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var peptide = allPeptides[0];

            var mutatedSequence = RandomlyMutateAminoAcid(peptide.FullSequence);
            var proteins = allProteins.Where(p => p.BaseSequence.Contains(peptide.BaseSeq)).ToList();
            foreach (var protein in proteins)
            {
                int index = protein.BaseSequence.IndexOf(peptide.FullSequence);
                string newProteinSequence = null;
                if (index >= 0)
                {
                    newProteinSequence = protein.BaseSequence.Substring(0, index) + mutatedSequence + protein.BaseSequence.Substring(index + mutatedSequence.Length);
                }
                var modifiedProtein = new Protein(protein, newProteinSequence);
                allProteins.Remove(protein);
                allProteins.Add(modifiedProtein);
                if (newProteinSequence.Length != protein.BaseSequence.Length || newProteinSequence == protein.BaseSequence)
                {
                    int stop = 0;
                }
            }
            var outPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\modifiedYeastFastaOnePeptide.fasta";
            ProteinDbWriter.WriteFastaDatabase(allProteins, outPath, " ");
        }

        [Test]
        public static void CheckModifiedDb()
        {
            var xml = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var xmlDb = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\2025-08-10-10-08-46\Task1-GPTMDTask\modifiedYeastFastaOnePeptideGPTMD.xml";
            var fastaDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            var allProteins = ProteinDbLoader.LoadProteinXML(xmlDb, true, DecoyType.None, null, false,
                new List<string>(), out Dictionary<string, Modification> un);
            //var allProteins = ProteinDbLoader.LoadProteinFasta(fastaDb, true, DecoyType.None, false, out List<string> errors);

            //ATDGGAHGVINVSVSEAAIEASTR
            //ATDGGAHGVINVFVSEAAIEASTR
            var proteins = allProteins.Where(p => p.BaseSequence.Contains("ATDGGAHGVINVFVSEAAIEASTR")).ToList();
            var mods = proteins[0].OriginalNonVariantModifications.Values.SelectMany(v => v).ToList();
            var psmPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\2025-08-10-10-08-46\Task2-SearchTask\AllPSMs.psmtsv";
            var allPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var psms = allPsms.Where(p => p.BaseSeq == "ATDGGAHGVINVFVSEAAIEASTR").ToList();
            var outPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\oneProtein.xml";
            var dic = ProteinDbWriter.WriteXmlDatabase(null, proteins, outPath);
            var peptides = proteins[0].Digest(new DigestionParams(), null, null);
            var modPeptides = peptides.Where(p => p.AllModsOneIsNterminus.Count > 0).ToList();

            ////Q08581
            //var testXml = @"E:\Aneuploidy\DDA\071525\1611_E1-8_cali-generalGPTMD_noTrunc\Task2-GPTMDTask\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16GPTMD.xml";
            //var testProteins = ProteinDbLoader.LoadProteinXML(testXml, true, DecoyType.None, null, false,
            //                   new List<string>(), out Dictionary<string, Modification> un2);
            //var testProtein = testProteins.Where(p => p.Accession == "Q08581").ToList();
            //var peptides2 = testProtein[0].Digest(new DigestionParams(), null, null);

            //test
            var tryxml = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\xml_test.xml";
            var tryProteins = ProteinDbLoader.LoadProteinXML(tryxml, true, DecoyType.None, null, false,
                                              new List<string>(), out Dictionary<string, Modification> un3);
        }

        public static string RandomlyMutateAminoAcid(string sequence)
        {
            var subModsPath = @"E:\GitClones\mzLib\mzLib\Omics\Modifications\substitutions.txt";
            var subMods = PtmListLoader.ReadModsFromFile(subModsPath, out var errorMods).Where(m => m.ModificationType == "1 nucleotide substitution").ToList();
            var subDic = new Dictionary<char, List<char>>();

            foreach (var sub in subMods)
            {
                var aa = sub.OriginalId.Split("->")[0][0];
                var subAA = sub.OriginalId.Split("->")[1];
                subDic[aa] = subDic.ContainsKey(aa) ? subDic[aa] : new List<char>();
                subDic[aa].Add(subAA[0]);
            }

            var rand = new Random(42);
            int pos = rand.Next(sequence.Length);
            char originalAA = sequence[pos];
            char newAA = originalAA;
            var doNotTouchAAs = new List<char> { 'C', 'M', 'K', 'R' };
            while (originalAA == newAA || doNotTouchAAs.Contains(newAA) || doNotTouchAAs.Contains(newAA))
            {
                pos = rand.Next(sequence.Length);
                originalAA = sequence[pos];
                var allMutants = subDic[originalAA];
                newAA = allMutants[rand.Next(allMutants.Count)];
            }
            char[] arr = sequence.ToCharArray();
            arr[pos] = newAA;
            return new string(arr);
        }

        [Test]
        public static void TestMakeNewProteins()
        {
            var psmPath = @"E:\DIA\Data\DIA_bu_250114\fasta\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var peptide = allPeptides.First(p => !p.FullSequence.Contains("["));

            Random rand = new Random();
            var protein = OriginalProteins.FirstOrDefault(p => p.Accession == peptide.ProteinAccession);
            var mutatedSequence = RandomlyMutateAminoAcid(peptide.FullSequence);
            var pos = peptide.StartAndEndResiduesInProtein.Trim('[', ']').Split(new[] { " to " }, StringSplitOptions.None);
            var newProteinSequence = protein.BaseSequence.Substring(0, int.Parse(pos[0]) - 1) + mutatedSequence + protein.BaseSequence.Substring(int.Parse(pos[1]));
            var modifiedProtein = new Protein(protein, newProteinSequence);
            Assert.That(newProteinSequence.Length, Is.EqualTo(protein.BaseSequence.Length));
        }

        [Test]
        public static void ModifiedSearch()
        {
            var allGptmdPsmPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\modifiedYeastFasta100_noAmb_NicGPTMD\Task2-SearchTask\AllPSMs.psmtsv";
            var allGptmdPsms = SpectrumMatchTsvReader.ReadPsmTsv(allGptmdPsmPath, out List<string> _).ToList();//.Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01)
            var subResults = new List<SubResult>();
            foreach (var peptide in SelectedPeptides)
            {
                var matchedGptmdPsm = allGptmdPsms.FirstOrDefault(p => p.Ms2ScanNumber == peptide.Ms2ScanNumber && p.PrecursorCharge == peptide.PrecursorCharge
                && Math.Round(p.PrecursorMass, 0) == Math.Round(peptide.PrecursorMass, 0)); //&& Math.Round(p.PrecursorMass, 0) == Math.Round(peptide.PrecursorMass, 0)
                if (matchedGptmdPsm != null)
                {
                    var subResult = new SubResult(peptide, matchedGptmdPsm, null, null);
                    subResults.Add(subResult);
                }
            }
            var subResultFile = new SubResultFile { Results = subResults };
            var path = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\subResultsFasta100subOnly_noOtherMods_noAmb.tsv";
            subResultFile.WriteResults(path);
        }

        [Test]
        public static void CheckModifiedGptmdDb()
        {
            var gptmdDbPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\modifiedYeastFasta100_noAmb\Task1-GPTMDTask\modifiedYeastFasta100_noAmbGPTMD.xml";
            var gptmdProteins = ProteinDbLoader.LoadProteinXML(gptmdDbPath, true, DecoyType.None, ModificationConverter.AllKnownMods, false, null, out Dictionary<string, Modification> um);
            var foundProteins = ModifiedProteins.Where(p => gptmdProteins.Any(pro => pro.BaseSequence == p.BaseSequence)).ToList();
            var allGptmdPeptides = new List<PeptideWithSetModifications>();
            foreach (var protein in gptmdProteins)
            {
                var peptides = protein.Digest(new DigestionParams(), null, null);
                allGptmdPeptides.AddRange(peptides);
            }
            var subSequences = allGptmdPeptides.Where(p => p.FullSequence.Contains("substitution")).ToList();
            var found = new List<string>();
            foreach (var seq in ModifiedSequences)
            {
                if (allGptmdPeptides.Any(p => p.BaseSequence == seq))
                {
                    found.Add(seq);
                }
            }

            var modifiedDbPath = @"E:\Aneuploidy\DDA\071525\RtPredictionResults\modifiedYeastFasta100_noAmb.fasta";
            var allProteinsAfterMutation = ProteinDbLoader.LoadProteinFasta(modifiedDbPath, true, DecoyType.None, false, out List<string> errors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex);
            var allPep = new List<PeptideWithSetModifications>();
            foreach (var protein in allProteinsAfterMutation)
            {
                var peptides = protein.Digest(new DigestionParams(), null, null);
                allPep.AddRange(peptides);
            }
            var found2 = new List<string>();
            foreach (var seq in ModifiedSequences)
            {
                if (allPep.Any(p => p.BaseSequence == seq))
                {
                    found2.Add(seq);
                }
            }
        }

        public static string ParseSubstitutedFullSequence(string input)
        {
            string pattern = @"\[1 nucleotide substitution:([A-Z])->([A-Z]) on ([A-Z])\]";
            var match = Regex.Match(input, pattern);
            if (match.Success)
            {
                string from = match.Groups[3].Value; // Z (original)
                string to = match.Groups[2].Value;   // Y (substitute)
                int patternIndex = match.Index;
                int replaceIndex = input.LastIndexOf(from, patternIndex);
                if (replaceIndex != -1)
                {
                    // Replace the first occurrence of Z before the pattern with Y
                    input = input.Remove(replaceIndex, 1).Insert(replaceIndex, to);
                }
                // Remove the pattern
                input = Regex.Replace(input, pattern, "");
            }
            return input;
        }

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

    public class SubResult
    {
        public int? Index { get; set; }
        public int Ms2ScanNumber { get; set; }
        public double PrecursorMass { get; set; }
        public int PrecursorCharge { get; set; }
        public string OriginalFullSequence { get; set; }

        public string? MutatedSequence { get; set; }
        public string ReIdentifiedFullSequence { get; set; }
        public double OriginalScore { get; set; }
        public double ReIdentifiedScore { get; set; }

        public SubResult(PsmFromTsv origPsm, PsmFromTsv newPsm, int? index, string? mutatedSequence)
        {
            Index = index;
            Ms2ScanNumber = origPsm.Ms2ScanNumber;
            PrecursorMass = origPsm.PrecursorMass;
            PrecursorCharge = origPsm.PrecursorCharge;
            OriginalFullSequence = origPsm.FullSequence;
            MutatedSequence = mutatedSequence;
            ReIdentifiedFullSequence = newPsm?.FullSequence;
            OriginalScore = origPsm.Score;
            ReIdentifiedScore = newPsm?.Score ?? 0;
        }
        public SubResult() { }
    }

    public class SubResultFile : ResultFile<SubResult>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public SubResultFile() : base() { }
        public SubResultFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<SubResult>().ToList();
        }

        public string FullFileName { get; set; }
        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<SubResult>();
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
