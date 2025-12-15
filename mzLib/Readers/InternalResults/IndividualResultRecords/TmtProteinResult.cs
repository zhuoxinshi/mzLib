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

        [Name("Protein QValue")]
        public double ProteinQValue { get; set; }
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
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

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

        public override SupportedFileType FileType => SupportedFileType.Tsv_FlashDeconv;
        public override Software Software { get; set; }
    }
}
