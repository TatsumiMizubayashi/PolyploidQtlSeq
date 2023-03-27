namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// Fastqファイルコレクション
    /// </summary>
    internal class FastqFiles
    {
        private readonly string[] _sortedFastqFilePaths;

        /// <summary>
        /// Fastqファイルコレクションを作成する。
        /// </summary>
        /// <param name="fastqFilePaths">Fastq File Path配列</param>
        public FastqFiles(string[] fastqFilePaths)
        {
            if (fastqFilePaths.Length == 0) throw new ArgumentException($"The number of files in {nameof(fastqFilePaths)} is 0.");
            if (fastqFilePaths.Length % 2 == 1) throw new ArgumentException($"The number of {nameof(fastqFilePaths)} files is odd.");

            _sortedFastqFilePaths = fastqFilePaths;
        }

        /// <summary>
        /// Fastqファイルペアに変換する。
        /// </summary>
        /// <returns>Fastqファイルペア配列</returns>
        public FastqFilePair[] ToFastqFilePairs()
        {
            return _sortedFastqFilePaths
                .Chunk(2)
                .Select(x => new FastqFilePair(x[0], x[1]))
                .ToArray();
        }
    }
}
