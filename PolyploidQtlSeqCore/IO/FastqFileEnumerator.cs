namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// Fastqファイルペアを列挙する
    /// </summary>
    internal static class FastqFilePairEnumerator
    {
        private static readonly string[] _fastqFilePatterns =
        [
            "*.fq.gz",
            "*.fastq.gz"
        ];

        /// <summary>
        /// FastqFilePairを列挙する。
        /// </summary>
        /// <param name="dirPath">調査ディレクトリのPath</param>
        /// <returns>Fastqファイルペア配列</returns>
        public static FastqFilePair[] Enumerate(string dirPath)
        {
            var fastqPairList = new List<FastqFilePair>();

            foreach (var fastqFilePattern in _fastqFilePatterns)
            {
                var sortedFastqFilePaths = FileEnumerator.Enumerate(dirPath, fastqFilePattern);
                if (sortedFastqFilePaths.Length == 0) continue;

                var fastqFiles = new FastqFiles(sortedFastqFilePaths);
                fastqPairList.AddRange(fastqFiles.ToFastqFilePairs());
            }

            return [.. fastqPairList];
        }
    }
}
