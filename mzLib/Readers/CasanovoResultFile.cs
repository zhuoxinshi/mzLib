using CsvHelper.Configuration;
using CsvHelper;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static TorchSharp.torch.nn;

namespace Readers
{
    public class CasanovoResult
    {
        public string DenovoSequence { get; set; }
        public string FullSequenceFromMM { get; set; }
        public string ScanNumber { get; set; }
        public string Scores { get; set; }

        public CasanovoResult(string line)
        {
            var columns = line.Split('\t');
            DenovoSequence = columns[1];
            string[] parts = columns[14].Split(new[] { ": ", " " }, 3, StringSplitOptions.None);
            FullSequenceFromMM = parts[2];
            ScanNumber = parts[1];
            Scores = columns[columns.Length - 1];
        }
        public CasanovoResult() { }
    }

    public class CasanovoResultFile : ResultFile<CasanovoResult>, IResultFile
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        public static List<CasanovoResult> ReadInCasanovoResults(string path)
        {
            var allResults = new List<CasanovoResult>();
            foreach (var line in File.ReadLines(path))
            {
                // Only keep data table rows (e.g., PSM section starts with 'PSH' or 'PSM')
                if (line.StartsWith("PSM"))
                {
                    var result = new CasanovoResult(line);
                    allResults.Add(result);
                }
            }
            return allResults;
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
