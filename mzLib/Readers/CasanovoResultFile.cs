using CsvHelper.Configuration;
using CsvHelper;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class CasanovoResult
    {
        public string DenovoSequence { get; set; }
        public string FullSequenceFromMM { get; set; }
        public string Scores { get; set; }

        public CasanovoResult(string line)
        {
            
        }
        public CasanovoResult() { }
    }

    public class CasanovoResultFile : ResultFile<CasanovoResult>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public static void ReadInCasanovoResults(string path)
        {
            var dataLines = new List<string>();

            foreach (var line in File.ReadLines(path))
            {
                // Skip metadata lines (those starting with 'MTD', 'COM', etc.)
                if (line.StartsWith("MTD") || line.StartsWith("COM") || string.IsNullOrWhiteSpace(line))
                    continue;

                // Only keep data table rows (e.g., PSM section starts with 'PSH' or 'PSM')
                if (line.StartsWith("PSH") || line.StartsWith("PSM")) ;
                    //var columns = line.Split('\t');
            }
        }

        public CasanovoResultFile() : base() { }
        public CasanovoResultFile(string filePath) : base(filePath, Software.Unspecified) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), CsvConfiguration);
            Results = csv.GetRecords<CasanovoResult>().ToList();
        }

        public string FullFileName { get; set; }
        public override void WriteResults(string outputPath)
        {
            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), CsvConfiguration);

            csv.WriteHeader<CasanovoResult>();
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
