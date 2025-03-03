namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// Fastqファイルペア情報
    /// </summary>
    internal class FastqFilePair
    {
        private static readonly string _delimiter = ".";
        private static readonly char[] _splitter = ['.', '_'];

        /// <summary>
        /// Fastqファイルペアを作成する。
        /// </summary>
        /// <param name="fastq1Path">Fastq1ファイルのPath</param>
        /// <param name="fastq2Path">Fastq2ファイルのPath</param>
        public FastqFilePair(string fastq1Path, string fastq2Path)
        {
            if (string.IsNullOrWhiteSpace(fastq1Path)) throw new ArgumentException(null, nameof(fastq1Path));
            if (string.IsNullOrWhiteSpace(fastq2Path)) throw new ArgumentException(null, nameof(fastq2Path));

            Fastq1Path = fastq1Path;
            Fastq2Path = fastq2Path;
            BaseName = GetBaseName();
        }

        /// <summary>
        /// Fastq1ファイルのPathを取得する。
        /// </summary>
        public string Fastq1Path { get; }

        /// <summary>
        /// Fastq2ファイルのPathを取得する。
        /// </summary>
        public string Fastq2Path { get; }

        /// <summary>
        /// Fastqのベース名を取得する。
        /// </summary>
        public string BaseName { get; }

        private string GetBaseName()
        {
            var fileName1 = Path.GetFileName(Fastq1Path);
            var fileName2 = Path.GetFileName(Fastq2Path);

            var fileName1Parts = fileName1.Split(_splitter);
            var fileName2Parts = fileName2.Split(_splitter);

            var samePartQuery = fileName1Parts.Zip(fileName2Parts, (a, b) => new { a, b })
                .TakeWhile(x => x.a == x.b)
                .Select(x => x.a);
            return string.Join(_delimiter, samePartQuery);
        }
    }
}
