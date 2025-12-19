using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Readers
{
    public class CometPsm
    {
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        [Name("start")]
        public int StartScan { get; set; }

        [Name("end")]
        public int EndScan { get; set; }

        [Name("Z")]
        public int Charge { get; set; }

        [Name("expM")]
        public double ExperimentalMass { get; set; }

        [Name("SpRank")]
        public int SpRank { get; set; }

        [Name("theoM")]
        public double TheoreticalMass { get; set; }

        [Name("deltaCN")]
        public double DeltaCN { get; set; }

        [Name("Xcorr")]
        public double Xcorr { get; set; }

        [Name("Sequence")]
        public string Sequence { get; set; }

        [Name("Loci")]
        public string Loci { get; set; }

        [Name("NewDeltaCN")]
        public double NewDeltaCN { get; set; }

        [Name("ISBDisc")]
        public double? ISBDisc { get; set; }

        [Name("NewDisc")]
        public double NewDisc { get; set; }

        [Name("ntt")]
        public int Ntt { get; set; }

        [Name("ForR")]
        public string ForwardOrReverse { get; set; }  // F or R

        [Name("Sp")]
        public double SpScore { get; set; }

        public CometPsm() { }

        public static string ConvertMMSeqToCometSeq(string mmSeq)
        {
            string removeOxidation = Regex.Replace(mmSeq, @"\[[^\]]*oxidation[^\]]*\]", "", RegexOptions.IgnoreCase);
            string cometSeq = $"{removeOxidation[0]}.{removeOxidation.Substring(1, removeOxidation.Length - 2)}.{removeOxidation[^1]}";
            return cometSeq;
        }

        public static string TrimAccession(string accession)
        {
            return accession.Split(' ')[0];
        }

        public static CometPsm ConvertPsmTsvToCometPsm(PsmFromTsv psmTsv)
        {
            var cometPsm = new CometPsm();
            cometPsm.StartScan = psmTsv.Ms2ScanNumber;
            cometPsm.EndScan = psmTsv.Ms2ScanNumber;
            cometPsm.Charge = psmTsv.PrecursorCharge;
            cometPsm.ExperimentalMass = psmTsv.PrecursorMass;
            cometPsm.TheoreticalMass = double.Parse(psmTsv.PeptideMonoMass);
            cometPsm.DeltaCN = psmTsv.DeltaScore.Value;
            cometPsm.Xcorr = psmTsv.Score;
            cometPsm.Sequence = ConvertMMSeqToCometSeq(psmTsv.BaseSeq);
            cometPsm.SpRank = 1; 
            cometPsm.Loci = TrimAccession(psmTsv.ProteinAccession);
            cometPsm.NewDeltaCN = psmTsv.DeltaScore.Value;
            cometPsm.ISBDisc = null; // Not available
            cometPsm.NewDisc = psmTsv.Score;
            cometPsm.Ntt = 2;
            cometPsm.ForwardOrReverse = psmTsv.IsDecoy ? "R" : "F";
            cometPsm.SpScore = psmTsv.Score;

            return cometPsm;
        }
    }

    public class CometPsmFile : ResultFile<CometPsm>, IResultFile
    {
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CometPsm.CsvConfiguration);
            Results = csv.GetRecords<CometPsm>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CometPsm.CsvConfiguration);

            csv.WriteHeader<CometPsm>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public CometPsmFile() : base() { }
        public CometPsmFile(string filePath) : base(filePath) { }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }
}
