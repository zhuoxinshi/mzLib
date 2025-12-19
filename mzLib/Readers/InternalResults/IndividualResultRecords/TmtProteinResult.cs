using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{

    public class TmtPsmResult
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = CsvHelper.Configuration.TrimOptions.Trim,
        };

        [Name("Scan Number")]
        public int ScanNumber { get; set; }
        [Name("Base Sequence")]
        public string BaseSequence { get; set; }

        [Name("126")]
        public double Chann126 { get; set; }

        [Name("127N")]
        public double Chann127N { get; set; }

        [Name("127C")]
        public double Chann127C { get; set; }

        [Name("128N")]
        public double Chann128N { get; set; }
        [Name("128C")]
        public double Chann128C { get; set; }

        [Name("129N")]
        public double Chann129N { get; set; }

        [Name("129C")]
        public double Chann129C { get; set; }

        [Name("130N")]
        public double Chann130N { get; set; }
        [Name("130C")]
        public double Chann130C{ get; set; }
        [Name("131N")]
        public double Chann131N { get; set; }
        [Name("131C")]
        public double Chann131C { get; set; }

        public TmtPsmResult() { }
    }

    public class TmtPsmResultFile : ResultFile<TmtPsmResult>, IResultFile
    {
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), TmtProteinResult.CsvConfiguration);
            Results = csv.GetRecords<TmtPsmResult>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), TmtProteinResult.CsvConfiguration);

            csv.WriteHeader<TmtPsmResult>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public TmtPsmResultFile() : base() { }
        public TmtPsmResultFile(string filePath) : base(filePath) { }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }

    public class TmtProteinResult
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = CsvHelper.Configuration.TrimOptions.Trim,
        };

        [Name("Protein Accession")]
        public string Accession { get; set; }

        [Name("Number of Proteins in Group")]
        public string NumProteinsInGroup { get; set; }

        [Name("Unique Peptides")]
        public string UniquePeptides { get; set; }

        [Name("Shared Peptides")]
        public string SharedPeptides { get; set; }
        [Name("Number of Unique Peptides")]
        public int NumUniquePeptides { get; set; }
        [Name("Protein Decoy/Contaminant/Target")]
        public string TargetDecoy { get; set; }

        [Name("Protein QValue")]
        public double ProteinQValue { get; set; }
        [Optional]
        [Name("126")]
        public double Chann126 { get; set; }
        [Optional]
        [Name("127N")]
        public double Chann127N { get; set; }
        [Optional]
        [Name("127C")]
        public double Chann127C { get; set; }
        [Optional]
        [Name("128N")]
        public double Chann128N { get; set; }
        [Optional]
        [Name("128C")]
        public double Chann128C { get; set; }
        [Optional]
        [Name("129N")]
        public double Chann129N { get; set; }
        [Optional]
        [Name("129C")]
        public double Chann129C { get; set; }
        [Optional]
        [Name("130N")]
        public double Chann130N { get; set; }
        [Optional]
        [Name("130C")]
        public double Chann130C { get; set; }
        [Optional]
        [Name("131N")]
        public double Chann131N { get; set; }
        [Optional]
        [Name("131C")]
        public double Chann131C { get; set; }
        [Ignore]
        public List<string> UniquePeptideList => UniquePeptides.Split('|').ToList();
        public TmtProteinResult() { }

        public void UpdateReporterIntensities(SpectrumMatchFromTsv psm)
        {
            Chann126 += psm.Chann126;
            Chann127N += psm.Chann127N;
            Chann127C += psm.Chann127C;
            Chann128N += psm.Chann128N;
            Chann128C += psm.Chann128C;
            Chann129N += psm.Chann129N;
            Chann129C += psm.Chann129C;
            Chann130N += psm.Chann130N;
            Chann130C += psm.Chann130C;
            Chann131N += psm.Chann131N;
            Chann131C += psm.Chann131C;
        }

        public static void FdrReanalysis(List<PsmFromTsv> psmTsvs)
        {
            var psmGroups = psmTsvs.Where(p => p.ChargeState <= 4).GroupBy(p => new { p.ChargeState, p.Notch, p.Oxidation });

            foreach (var group in psmGroups)
            {
                var sortedPsms = group.OrderByDescending(p => p.Score).ToList();
                double cumulativeTarget = 0;
                double cumulativeDecoy = 0;

                foreach (var psm in sortedPsms)
                {
                    if (psm.IsDecoy)
                    {
                        cumulativeDecoy++;
                    }
                    else
                    {
                        cumulativeTarget++;
                    }

                    psm.PAW_qvalue = cumulativeDecoy / (double)cumulativeTarget;
                }
            }
        }
    }

    public class TmtProteinResultFile : ResultFile<TmtProteinResult>, IResultFile
    {
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), TmtProteinResult.CsvConfiguration);
            Results = csv.GetRecords<TmtProteinResult>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), TmtProteinResult.CsvConfiguration);

            csv.WriteHeader<TmtProteinResult>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public TmtProteinResultFile() : base() { }
        public TmtProteinResultFile(string filePath) : base(filePath) { }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }
}
