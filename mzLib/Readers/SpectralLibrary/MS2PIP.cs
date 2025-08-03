using CsvHelper.Configuration;
using CsvHelper;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Readers.SpectralLibrary
{
    public class Ms2PipInput
    {
        public string peptidoform { get; set; }
        public string spectrum_id { get; set; }

        public Ms2PipInput(PsmFromTsv psmTsv)
        {
            string updatedFullSeq = psmTsv.FullSequence;
            if (psmTsv.FullSequence.Contains("substitution"))
            {
                updatedFullSeq = ParseSubstitutedFullSequence(psmTsv.FullSequence);
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

        public static string ParseSubstitutedFullSequence(string fullSequence)
        {
            var subMatch = Regex.Match(fullSequence, @"\[(\d+)[^\:]*:([A-Z])->([A-Z])");
            if (!subMatch.Success)
                return fullSequence; // No substitution found

            int position = int.Parse(subMatch.Groups[1].Value) - 1; // 0-based
            char newAminoAcid = subMatch.Groups[3].Value[0];

            // Remove only the substitution annotation
            string cleaned = Regex.Replace(fullSequence, @"\[\d+[^\]]*substitution:[A-Z]->[A-Z][^\]]*\]", "");

            // Find the index after the last closing bracket (end of all annotations)
            int seqStart = cleaned.LastIndexOf(']') + 1;
            if (seqStart < 0 || seqStart >= cleaned.Length)
                return cleaned; // No sequence found

            // The sequence is everything after the last bracket
            string prefix = cleaned.Substring(0, seqStart);
            string sequence = cleaned.Substring(seqStart);

            // Only substitute if the sequence is long enough
            if (sequence.Length > position)
            {
                char[] seqArray = sequence.ToCharArray();
                seqArray[position] = newAminoAcid;
                string modifiedSequence = new string(seqArray);
                return prefix + modifiedSequence;
            }

            // If not long enough, return cleaned string
            return cleaned;
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
